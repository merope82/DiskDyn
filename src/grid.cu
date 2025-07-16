#include <stdio.h>
#include "ddyn.cuh"

void Grid::null_gpu(){
    if ( n_all >0 ){
	int block = int(( n_all - 1.0 ) / 128.0 + 1.0);
	null_coord_driver<<<block,128>>>(n_all,gpu);
	synced = false;
    }
}

void Grid::null_flux(){
    configdata *cfg  = &cconfig;
    
    if ( cfg->xsize>0 && cfg->ysize>0 ){
	int block = int(( cfg->xsize*cfg->ysize - 1.0 ) / 128.0 + 1.0);
	null_flux_driver<<<block,128>>>(cfg->xsize*cfg->ysize,Fd);
    }
}

void Grid::null_flux_part(){
    if ( n_dust > 0 ){
	int block = int(( n_dust - 1.0 ) / 128.0 + 1.0);
	null_flux_part_driver<<<block,128>>>(n_dust,Fdust);
    }
}

void Grid::copy_sed(int j){
    configdata *cfg  = &cconfig;

    if ( n_dust > 0 ){
	int block = int(( n_dust - 1.0 ) / 128.0 + 1.0);

	sum_flux_driver<<<block,128,128*sizeof(double)>>>(n_dust,Fdust);

	// Download Fluxsum
	double flux=0.0;
	cudaMemcpyFromSymbol(&flux,doubletmp,sizeof(double));
	cfg->inexc[j] = flux;
    }
}

void Grid::calc_flux(int wj){
    configdata	*cfg  = &cconfig;

    printcalcflux(wj);
    if ( n_dust > 0 ){
	int n_grav = n_mass + cfg->n_planets;
	int block = int(( n_dust - 1.0 ) / 128.0 + 1.0);	// Only massless dust particles emit!

	calc_flux_driver<<<block,128>>>(n_grav,n_dust,gpu,cfg->dist,Rmin,Rmax,Qabs,Qsca,Qpfc,Td,Fd,w,in,\
                cfg->n_wav,cfg->n_dust,cfg->n_theta,wj,cfg->longitude,cfg->inc,cfg->PA,cfg->xsize,cfg->ysize,cfg->fov);
    }
}

void Grid::calc_flux_sed(int wj){
    configdata	*cfg  = &cconfig;

    if ( n_dust > 0 ){
	int n_grav = n_mass + cfg->n_planets;
	int block = int(( n_dust - 1.0) / 64.0 + 1.0);		// Only massless dust particles emit!

	calc_flux_sed_driver<<<block,64>>>(n_grav,n_dust,gpu,cfg->dist,Rmin,Rmax,Qabs,Qsca,Qpfc,Td,Fdust,w,in,\
                cfg->n_wav,cfg->n_dust,cfg->n_theta,wj,cfg->longitude,cfg->inc,cfg->PA);
    }
}

void Grid::setup(){
    configdata	*cfg  = &cconfig;
    double	x,y,z,vx,vy,vz,N,m;
    int		totest;

    if ( cfg->verb>1 ) printf("Interpolating optical properties to grid\n");
    interpolate_dust_to_grid();
    interpolate_stellar_to_grid();
    calculate_beta();
    calculate_dust_temp();
    printsystemparams();

    seed   = -2;
    time   = 0;
    dt     = 0;
    tlastd = 0;
    tlastf = 0;
    tlasts = 0;

    n_dust = 0;
    n_mass = 0;

    cfg->seed = seed;
    for ( int j=0 ; j< cfg->n_disks ; j++ ){ n_mass += cfg->disks[j].Ngrav; }
    for ( int j=0 ; j< cfg->n_blobs ; j++ ){ n_mass += cfg->blobs[j].Ngrav; }

    n_all  = cfg->n_planets + n_mass;
    totest = n_all;

    for ( int j=0 ; j< cfg->n_disks ; j++ ){ totest += cfg->disks[j].Nsmall; }
    for ( int j=0 ; j< cfg->n_blobs ; j++ ){ totest += cfg->blobs[j].Nsmall; }

    // allochost needs n_all
    allochost();

    // Add planets to array first    
    n_all = 0;

    for ( int j=0 ; j < cfg->n_planets ; j++ ){
	add_to_host(cfg->planets[j].x,cfg->planets[j].y,cfg->planets[j].z,\
		    cfg->planets[j].vx,cfg->planets[j].vy,cfg->planets[j].vz,\
		    1,cfg->planets[j].m,-1,totest);
    }

    for ( int j=0 ; j < cfg->n_disks ; j++ ){
	long double mmax = cfg->disks[j].mmax;
	for ( int l=0 ; l<cfg->disks[j].Ngrav ; l++ ){
	    getbod(&mmax,&m,cfg->disks[j].C,1.0-cfg->disks[j].eta);
    	    genorb_disk(j,m,&seed,&x,&y,&z,&vx,&vy,&vz,-1);
	    add_to_host(x,y,z,vx,vy,vz,1,m,-1,totest);
	    if ( l==cfg->disks[j].Ngrav-1 ){ 
		printmin((double)mmax,j,0); 
	    }
	}
    }

    for ( int j=0 ; j < cfg->n_blobs ; j++ ){
	long double mmax = cfg->blobs[j].mmax;
	for ( int l=0 ; l<cfg->blobs[j].Ngrav ; l++ ){
	    getbod(&mmax,&m,cfg->blobs[j].C,1.0-cfg->blobs[j].eta);
	    genorb_blob(j,m,&seed,&x,&y,&z,&vx,&vy,&vz,-1);
	    add_to_host(x,y,z,vx,vy,vz,1,m,-1,totest);
	    if ( l==cfg->blobs[j].Ngrav-1 ){ 
		printmin((double)mmax,j,1); 
	    }
	}
    }

    // Add "massless" particles
    totest = n_all;

    for ( int j=0 ; j < cfg->n_disks ; j++ ){
	// Redefine limit to grid
	double mmax = cfg->mmin * exp(round(log(cfg->disks[j].dmmax/cfg->mmin)/cfg->f)*cfg->f);
	double mmin = cfg->mmin * exp(round(log(cfg->disks[j].mmin /cfg->mmin)/cfg->f)*cfg->f);
	int Ngrid = round(log(mmax/mmin)/log(mgrid)+1);
	int Npg   = round((double)(cfg->disks[j].Nsmall)/((double)Ngrid));
        reallochost(Npg*Ngrid);
	n_dust += Npg*Ngrid;
	totest += Npg*Ngrid;
	for ( int l=0 ; l < Ngrid ; l++ ){
	    m = mmin * exp((double)l*cfg->f);
	    double size = pow(m*3.0/4.0/PI/cfg->rho,1.0/3.0);
	    int id = round(log(m/cfg->mmin)/cfg->f);
	    // Safety switches
		if ( size < cfg->dmin ) id = -1;
		if ( size - cfg->dmax > 1e-9 ) id = -1;
	    double in = mmin * exp(((double)l-0.5)*cfg->f);
	    double ax = mmin * exp(((double)l+0.5)*cfg->f);
	    m *= Grav/MSun;
    	    if      ( cfg->disks[j].eta==0 ) N = cfg->disks[j].C * (ax-in);
    	    else if ( cfg->disks[j].eta==1 ) N = cfg->disks[j].C * log(ax/in);
    	    else		 	     N = cfg->disks[j].C * (pow(ax,1.0-cfg->disks[j].eta)-pow(in,1.0-cfg->disks[j].eta))/(1.0-cfg->disks[j].eta);
	    for ( int k=0 ; k<Npg ; k++ ){
		genorb_disk(j,m,&seed,&x,&y,&z,&vx,&vy,&vz,id);
		add_to_host(x,y,z,vx,vy,vz,N/Npg,m,id,totest);
    }}}

    for ( int j=0 ; j < cfg->n_blobs ; j++ ){
	double mmax = cfg->mmin * exp(round(log(cfg->blobs[j].dmmax/cfg->mmin)/cfg->f)*cfg->f);
	double mmin = cfg->mmin * exp(round(log(cfg->blobs[j].mmin /cfg->mmin)/cfg->f)*cfg->f);
	int Ngrid = round(log(mmax/mmin)/log(mgrid)+1);
	int Npg   = round((double)(cfg->blobs[j].Nsmall)/((double)Ngrid));
        reallochost(Npg*Ngrid);
	n_dust += Npg*Ngrid;
	totest += Npg*Ngrid;
	for ( int l=0 ; l < Ngrid ; l++ ){
	    m = mmin * exp((double)l*cfg->f);
	    double size = pow(m*3.0/4.0/PI/cfg->rho,1.0/3.0);
	    int id = round(log(m/cfg->mmin)/cfg->f);
		if ( size < cfg->dmin ) id = -1;
		if ( size - cfg->dmax > 1e-9 ) id = -1;
	    double in = mmin * exp(((double)l-0.5)*cfg->f);
	    double ax = mmin * exp(((double)l+0.5)*cfg->f);
	    m *= Grav/MSun;
    	    if      ( cfg->blobs[j].eta==0 ) N = cfg->blobs[j].C * (ax-in);
    	    else if ( cfg->blobs[j].eta==1 ) N = cfg->blobs[j].C * log(ax/in);
    	    else			     N = cfg->blobs[j].C * (pow(ax,1.0-cfg->blobs[j].eta)-pow(in,1.0-cfg->blobs[j].eta))/(1.0-cfg->blobs[j].eta);
	    for ( int k=0 ; k<Npg ; k++ ){
		genorb_blob(j,m,&seed,&x,&y,&z,&vx,&vy,&vz,id);
		add_to_host(x,y,z,vx,vy,vz,N/Npg,m,id,totest);
    }}}

    if ( cfg->verb>1 ) printf("Particles set up and allocated on RAM\n");

    checkmemusage(n_all,n_dust);
    printmemall();
    allocgpu();
    null_gpu();

    if ( cfg->verb>1 ) printf("Syncing data up to GPU\n");
    sync_to_gpu();
    sync_dust_to_gpu();

    if ( cfg->verb>1 ) printf("Starting evolution\n\n");
}

void Grid::read_in(char *model){
    FILE	*fr;
    int		i,t,itemp,tot,offset;
    double	dtemp;
    double	x,y,z,vx,vy,vz,N,m;
    char	buff[2048],*dat[64];
    configdata	*cfg=&cconfig;

    fr=fopen(model,"rb");
    if ( fr==NULL )        exit_with_usage(3);

    i=0;
    while ( ! feof(fr) ){
	if ( fgets(buff,2048,fr)==NULL )	break;
	t=tokenize_spaces(buff,dat,32);
	if ( t==0 )				continue;
	if ( i==48 ){
	    sscanf(dat[1],"%d",&itemp);
	    cfg->seed = itemp;
	}
	if ( i==98 ){
	    sscanf(dat[1],"%lg",&dtemp);
	    time = dtemp;
	}
	if ( i==99 ){
	    sscanf(dat[1],"%lg",&dtemp);
	    dt = dtemp;
	}
	if ( i==100 ){
	    sscanf(dat[1],"%lg",&dtemp);
	    tlastd = dtemp;
	}
	if ( i==101 ){
	    sscanf(dat[1],"%lg",&dtemp);
	    tlastf = dtemp;
	}
	if ( i==102 ){
	    sscanf(dat[1],"%lg",&dtemp);
	    tlasts = dtemp;
	}
	if ( i==103 ){
	    sscanf(dat[1],"%d",&itemp);
	    n_mass = itemp;
	}
	if ( i==104 ){
	    sscanf(dat[1],"%d",&itemp);
	    n_dust = itemp;
	}
	if ( i==105 ){
	    sscanf(dat[1],"%d",&itemp);
	    n_all = itemp;
	    allochost();
	    tot = n_all;
	    n_all = 0;
	}
	if ( i>=108 && i<108+cfg->n_planets ){
	    if ( strcmp(dat[0],"##")==0 ) offset=1;
	    else			  offset=0;
	    sscanf(dat[0+offset],"%lg",&x);
	    sscanf(dat[1+offset],"%lg",&y);
	    sscanf(dat[2+offset],"%lg",&z);
	    sscanf(dat[3+offset],"%lg",&vx);
	    sscanf(dat[4+offset],"%lg",&vy);
	    sscanf(dat[5+offset],"%lg",&vz);
	    sscanf(dat[6+offset],"%lg",&N);
	    sscanf(dat[7+offset],"%lg",&m);
	    add_to_host(x,y,z,vx,vy,vz,N,m*Grav/MSun,-1,tot);
	}
	if ( i>=108+cfg->n_planets+1 ){
	    if ( strcmp(dat[0],"##")==0 ) offset=1;
	    else			  offset=0;
	    sscanf(dat[0+offset],"%lg",&x);
	    sscanf(dat[1+offset],"%lg",&y);
	    sscanf(dat[2+offset],"%lg",&z);
	    sscanf(dat[3+offset],"%lg",&vx);
	    sscanf(dat[4+offset],"%lg",&vy);
	    sscanf(dat[5+offset],"%lg",&vz);
	    sscanf(dat[6+offset],"%lg",&N);
	    sscanf(dat[7+offset],"%lg",&m);
	    int id = round(log(m/cfg->mmin)/cfg->f);
	    if ( id > cfg->n_grid-1 ) id=-1;
	    add_to_host(x,y,z,vx,vy,vz,N,m*Grav/MSun,id,tot);
	}
	i++;
    }
    fclose(fr);

    if ( cfg->verb>1 ) printf("Particles set up and allocated on RAM\n");
    if ( cfg->verb>1 ) printf("Interpolating optical properties to grid\n");
    interpolate_dust_to_grid();
    interpolate_stellar_to_grid();
    calculate_beta();
    calculate_dust_temp();
    printsystemparams();

    checkmemusage(n_all,n_dust);
    printmemall();
    allocgpu();
    null_gpu();

    if ( cfg->verb>1 ) printf("Syncing data up to GPU\n");
    sync_to_gpu();
    sync_dust_to_gpu();

    if ( cfg->verb>1 ) printf("Starting evolution\n\n");
}

double Grid::det_step(){
    configdata	*cfg = &cconfig;
    double dtdynmin=1e137;
    
    // I could do it on the GPU, but I am lazy
    sync_from_gpu();

    if ( cfg->n_planets + n_mass > 1 ){
        double R =  pow((host.x[1]-host.x[0])*(host.x[1]-host.x[0])+\
	    	        (host.y[1]-host.y[0])*(host.y[1]-host.y[0])+\
		        (host.z[1]-host.z[0])*(host.z[1]-host.z[0]),3.0/2.0);
    	dtdynmin = PI2 * sqrt(R/(host.m[0]+host.m[1])) / cfg->n_step;
        for ( int i=2 ; i < cfg->n_planets + n_mass ; i++ ){
               R =  pow((host.x[i]-host.x[0])*(host.x[i]-host.x[0])+\
	    	        (host.y[i]-host.y[0])*(host.y[i]-host.y[0])+\
		        (host.z[i]-host.z[0])*(host.z[i]-host.z[0]),3.0/2.0);
	    double dttmp = PI2 * sqrt(R/(host.m[0]+host.m[i])) / cfg->n_step;
	    if ( dttmp < dtdynmin ) dtdynmin = dttmp;
	}
    }
    for ( int i=0 ; i < cfg->n_disks ; i++ ){
	    double dttmp = PI2 * sqrt(cfg->disks[i].Rin/(host.m[0])) / cfg->n_step;
	    if ( dttmp < dtdynmin ) dtdynmin = dttmp;
    }
    for ( int i=0 ; i < cfg->n_blobs ; i++ ){
	    double dttmp = PI2 * sqrt((cfg->blobs[i].a*(1.0-cfg->blobs[i].ecc)-cfg->blobs[i].sigma)/(host.m[0])) / cfg->n_step;
	    if ( dttmp < dtdynmin ) dtdynmin = dttmp;
    }

    return(dtdynmin);
}

void Grid::evolve(const char * Version){
    configdata *cfg = &cconfig;
    double dtmp=0.0;

    // Determine smallest dt; it will define the timestep.
    double dt = det_step();
    printtimestep(dt);
    double dttmp = dt;

    // Why twice? I don't think this is needed, whether param file or model input.
//    cfg->betaswconst = 3.0 * cfg->mloss * Cd * MSun * year * cfg->vsw / 32.0 / PI / cfg->planets[0].m / cfg->rho / AU / AU / AU;

    // Print initial array
    switch( cfg->model ){ 
	case 1: 
	    for ( int j=0 ; j<cfg->nwavs ; j++ ){	// Observing wavelengths
		null_flux();
		calc_flux(cfg->wj[j]);
		write_image(cfg->wj[j],Version);
		cudaDeviceSynchronize();
	    }
	    // SED
	    printsedcalc();
	    for ( int j=0 ; j<cfg->n_wav ; j++ ){	// All wavelengths
		cudaMemcpyToSymbol(doubletmp,&dtmp,sizeof(double));
		null_flux_part();
		calc_flux_sed(j);
		copy_sed(j);
		cudaDeviceSynchronize();
	    }
	    write_sed(Version);
	    writegrid(Version); 
    }

    // Getting the shortest possible time necessary
    if ( time + dttmp > tlastd + cfg->t_write ) dttmp = cfg->t_write + tlastd - time;
    if ( time + dttmp > tlastf + cfg->t_fits  ) dttmp = cfg->t_fits  + tlastf - time;
    if ( time + dttmp > tlasts + cfg->t_sed   ) dttmp = cfg->t_sed   + tlasts - time;
    if ( time + dttmp > cfg->t_end )            dttmp = cfg->t_end   - time;

    // Now, let's evolve!
    do{
	// Evolve the massive bodies dynamically first
	massive_RK4(dttmp);
	small_RK4(dttmp);
	massive_final(dttmp);
	synced = false;

	time += dttmp;

	printstatus(time);

	if ( time >= tlasts + cfg->t_sed || time >= cfg->t_end ){ 
	    printnewline();
	    printsedcalc();
	    for ( int j=0 ; j<cfg->n_wav ; j++ ){	// All wavelengths
	        cudaMemcpyToSymbol(doubletmp,&dtmp,sizeof(double));
	        null_flux_part();
	        calc_flux_sed(j);
	        copy_sed(j);
	        cudaDeviceSynchronize();
	    }
	    write_sed(Version);
	    tlasts = time;
	    dttmp = dt;
	}

	if ( time >= tlastf + cfg->t_fits || time >= cfg->t_end ){ 
	    printnewline();
	    for ( int j=0 ; j<cfg->nwavs ; j++ ){
//		cudaDeviceSynchronize();
		null_flux();
//		cudaDeviceSynchronize();
		calc_flux(cfg->wj[j]);
//		cudaDeviceSynchronize();
		write_image(cfg->wj[j],Version);
		cudaDeviceSynchronize();
	    }
	    tlastf = time;
	    dttmp = dt;
	}

	if ( time >= tlastd + cfg->t_write || time >= cfg->t_end ){ 
	    printnewline();
	    tlastd = time;
	    writegrid(Version); 
	    dttmp  = dt;
	}

	if ( time + dttmp > tlastd + cfg->t_write ) dttmp = cfg->t_write + tlastd - time;
	if ( time + dttmp > tlastf + cfg->t_fits  ) dttmp = cfg->t_fits  + tlastf - time;
	if ( time + dttmp > tlasts + cfg->t_sed   ) dttmp = cfg->t_sed   + tlasts - time;
	if ( time + dttmp > cfg->t_end )            dttmp = cfg->t_end   - time;

    }while(time < cfg->t_end);
}

void Grid::writegrid(const char *ver){
    configdata *cfg = &cconfig;
    sync_from_gpu();

    FILE   *fw;
    char	out[512];

    printdata(time);
    sprintf(out,"%s_%.2f.dat",cfg->name_stub,time);
    fw=fopen(out,"wb");

    // Print header for read-in
    fprintf(fw,"#\n");
    fprintf(fw,"#Calculated_with_DiskDyn_Version %s\n",ver);
    fprintf(fw,"#\n");
    fprintf(fw,"#Config file parameters\n");
    fprintf(fw,"#Verbose\t\t\t%d\n",cfg->verb);
    fprintf(fw,"#GPUID\t\t\t\t%d\n",cfg->gpuid);
    fprintf(fw,"#Graphical_Output\t\t%d\n",cfg->gl);
    fprintf(fw,"#Time_Out_for_Graphics\t\t%.3f\n",cfg->t_gl);
    fprintf(fw,"#Time_Out_for_Data\t\t%.3f\n",cfg->t_write);
    fprintf(fw,"#Time_Out_for_SED\t\t%.3f\n",cfg->t_sed);
    fprintf(fw,"#Time_Out_for_FITS\t\t%.3f\n",cfg->t_fits);
    fprintf(fw,"#Time_End_for_model\t\t%.3f\n",cfg->t_end);
    fprintf(fw,"#N_step_per_orbit\t\t%d\n",cfg->n_step);
    fprintf(fw,"#Spectral-type_of_star\t\t%s\n",cfg->sptype);
    fprintf(fw,"#Temperature_of_star\t\t%.2f\n",cfg->teff);
    fprintf(fw,"#Radius_of_star\t\t\t%.5f\n",cfg->rstar*AU/RSun);
    fprintf(fw,"#Magnetic_field_of_star\t\t%.5e\n",cfg->Bstar*MSun*1.0e4/year/year);
    fprintf(fw,"#Stellar_rotation_of_star\t%.5f\n",365.2422/cfg->Wz);
    fprintf(fw,"#Stellar_mass_loss\t\t%.5e\n",cfg->mloss);
    fprintf(fw,"#Stellar_wind_speed\t\t%.5e\n",cfg->vsw);
    fprintf(fw,"#Include_rad_press\t\t%d\n",cfg->radQ);
    fprintf(fw,"#Include_SWD\t\t\t%d\n",cfg->SWQ);
    fprintf(fw,"#Include_Lorentz\t\t%d\n",cfg->FLQ);
    fprintf(fw,"#Bulk_density\t\t\t%.3f\n",cfg->rho);
    fprintf(fw,"#Electric_charge_coeff\t\t%.3f\n",cfg->kqe*year);
    fprintf(fw,"#Name_Stub\t\t\t%s\n",cfg->name_stub);
    fprintf(fw,"#Distance\t\t\t%.5f\n",cfg->dist/pc);
    fprintf(fw,"#Field_of_View\t\t\t%.6f\n",cfg->fov);
    fprintf(fw,"#Position_Angle\t\t\t%.5f\n",cfg->PA*180.0/PI);
    fprintf(fw,"#Inclination_Angle\t\t%.5f\n",cfg->inc*180.0/PI);
    fprintf(fw,"#Longitude_Angle\t\t%.5f\n",cfg->longitude*180.0/PI);
    fprintf(fw,"#Image_x_size\t\t\t%d\n",cfg->xsize);
    fprintf(fw,"#Image_y_size\t\t\t%d\n",cfg->ysize);
    fprintf(fw,"#Spectral-file\t\t\t%s\n",cfg->specfile);
    fprintf(fw,"#Minimum_particle_size\t\t%.8e\n",cfg->amin);
    fprintf(fw,"#Maximum_particle_size\t\t%.8e\n",cfg->amax);
    fprintf(fw,"#Minimum_particle_mass\t\t%.8e\n",cfg->mmin);
    fprintf(fw,"#Maximum_particle_mass\t\t%.8e\n",cfg->mmax);
    fprintf(fw,"#Number_of_dust_grids\t\t%d\n",cfg->n_grid);
    fprintf(fw,"#Mass_grid_ratio\t\t%.8e\n",cfg->f);
    fprintf(fw,"#Optical_composition_file\t%s\n",cfg->comp);
    fprintf(fw,"#Wavelengths_in_interpolation\t%d\n",cfg->n_wav);
    fprintf(fw,"#Number_of_thetas\t\t%d\n",cfg->n_theta);
    fprintf(fw,"#Number_of_dust_grids\t\t%d\n",cfg->n_dust);
    fprintf(fw,"#Number_of_planets\t\t%d\n",cfg->n_planets);
    fprintf(fw,"#Number_of_disks\t\t%d\n",cfg->n_disks);
    fprintf(fw,"#Number_of_blobs\t\t%d\n",cfg->n_blobs);
    fprintf(fw,"#Image_wavelengths\t");
    for ( int i=0 ; i<cfg->nwavs ; i++ ) fprintf(fw,"\t%.3f",cfg->wavs[i]*1e6);
    fprintf(fw,"\n");
    fprintf(fw,"#Seed_of_model\t\t\t%ld\n",cfg->seed);			// For Reference
    fprintf(fw,"#Disk_Rin\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t\t%.3f",cfg->disks[i].Rin);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_Rout\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t\t%.3f",cfg->disks[i].Rout);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_dhdr\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t\t%.4f",cfg->disks[i].dhdr);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_sigmae\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t\t%.4f",cfg->disks[i].sigmae);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_amax\t\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t%.6e",cfg->disks[i].amax);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_mmax\t\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t%.6e",cfg->disks[i].mmax);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_amin\t\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t%.6e",cfg->disks[i].amin);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_mmin\t\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t%.6e",cfg->disks[i].mmin);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_dmax\t\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t%.6e",cfg->disks[i].dmax);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_dmmax\t\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t%.6e",cfg->disks[i].dmmax);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_eta\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t\t%.4f",cfg->disks[i].eta);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_gammas\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t\t%.3f",cfg->disks[i].gammas);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_C\t\t\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t%.8e",cfg->disks[i].C);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_Mtot\t\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t%.8e",cfg->disks[i].Mtot);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_Ngrav\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t\t%ld",cfg->disks[i].Ngrav);
    fprintf(fw,"\n");
    fprintf(fw,"#Disk_Nsmall\t");
	for ( int i=0 ; i<cfg->n_disks ; i++ ) fprintf(fw,"\t\t%ld",cfg->disks[i].Nsmall);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_a\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.5f",cfg->blobs[i].a);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_ecc\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.5f",cfg->blobs[i].ecc);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_inc\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.5f",cfg->blobs[i].inc);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_w\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.5f",cfg->blobs[i].w);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_o\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.5f",cfg->blobs[i].o);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_nu\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.5f",cfg->blobs[i].nu);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_sigma\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.3f",cfg->blobs[i].sigma);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_amax\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t%.6e",cfg->blobs[i].amax);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_mmax\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t%.6e",cfg->blobs[i].mmax);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_amin\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t%.6e",cfg->blobs[i].amin);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_mmin\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t%.6e",cfg->blobs[i].mmin);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_dmax\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t%.6e",cfg->blobs[i].dmax);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_dmmax\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t%.6e",cfg->blobs[i].dmmax);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_eta\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.3f",cfg->blobs[i].eta);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_C\t\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t%.8e",cfg->blobs[i].C);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_Mtot\t\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t%.8e",cfg->blobs[i].Mtot);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_Ngrav\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%ld",cfg->blobs[i].Ngrav);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_Nsmall\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%ld",cfg->blobs[i].Nsmall);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_collcloud\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%d",cfg->blobs[i].collcloud);
    fprintf(fw,"\n");
    fprintf(fw,"#Blob_vkepp\t");
	for ( int i=0 ; i<cfg->n_blobs ; i++ ) fprintf(fw,"\t\t%.5f",cfg->blobs[i].vkepp);
    fprintf(fw,"\n");
    fprintf(fw,"#\n");
    fprintf(fw,"# t=0 positions and velocities of planets and Star\n");
    fprintf(fw,"#Planet_x\t");
	for ( int i=0 ; i<cfg->n_planets ; i++ ) fprintf(fw,"\t\t%.4f",cfg->planets[i].x);
    fprintf(fw,"\n");
    fprintf(fw,"#Planet_y\t");
	for ( int i=0 ; i<cfg->n_planets ; i++ ) fprintf(fw,"\t\t%.4f",cfg->planets[i].y);
    fprintf(fw,"\n");
    fprintf(fw,"#Planet_z\t");
	for ( int i=0 ; i<cfg->n_planets ; i++ ) fprintf(fw,"\t\t%.4f",cfg->planets[i].z);
    fprintf(fw,"\n");
    fprintf(fw,"#Planet_vx\t");
	for ( int i=0 ; i<cfg->n_planets ; i++ ) fprintf(fw,"\t\t%.4f",cfg->planets[i].vx);
    fprintf(fw,"\n");
    fprintf(fw,"#Planet_vy\t");
	for ( int i=0 ; i<cfg->n_planets ; i++ ) fprintf(fw,"\t\t%.4f",cfg->planets[i].vy);
    fprintf(fw,"\n");
    fprintf(fw,"#Planet_vz\t");
	for ( int i=0 ; i<cfg->n_planets ; i++ ) fprintf(fw,"\t\t%.4f",cfg->planets[i].vz);
    fprintf(fw,"\n");
    fprintf(fw,"#Planet_m\t\t");
	for ( int i=0 ; i<cfg->n_planets ; i++ ) fprintf(fw,"\t%.6e",cfg->planets[i].m);
    fprintf(fw,"\n");
    fprintf(fw,"#NR_used\t\t\t%d\n",nR);
    fprintf(fw,"#\n");
    fprintf(fw,"# Remaining parameters and data are for Grid class read-in\n");
    fprintf(fw,"#\n");
    fprintf(fw,"#Time_of_model\t\t\t%.8f\n",time);
    fprintf(fw,"#dt_of_model\t\t\t%.8f\n",dt);
    fprintf(fw,"#Last_t_of_data\t\t\t%.8f\n",tlastd);
    fprintf(fw,"#Last_t_of_fits\t\t\t%.8f\n",tlastf);
    fprintf(fw,"#Last_t_of_SED\t\t\t%.8f\n",tlasts);
    fprintf(fw,"#Number_of_grav_bodies\t\t%d\n",n_mass);
    fprintf(fw,"#Number_of_dust_particles\t%d\n",n_dust);
    fprintf(fw,"#Number_of_all_particles\t%d\n",n_all);

    // Print planets
    fprintf(fw,"#Data: x[AU] y[AU] z[AU] vx[AU/yr] vy[AU/yr] vz[AU/yr] N m[kg] r[m] beta\n");
    fprintf(fw,"#Planets:\n");
    for ( int i=0 ; i<cfg->n_planets ; i++ ){
	fprintf(fw,"% .8e % .8e % .8e % .8e % .8e % .8e %.12e %.8e %d\n",host.x[i],host.y[i],host.z[i],host.vx[i],host.vy[i],host.vz[i],1.0,host.m[i]*MSun/Grav,host.id[i]);
    }
    fprintf(fw,"#Bodies:\n");
    for ( int i=cfg->n_planets ; i<n_all ; i++ ){
	double R = host.x[i]*host.x[i]+host.y[i]*host.y[i]+host.z[i]*host.z[i];
	if ( R > Rmax*Rmax )             fprintf(fw,"## ");
	if ( R < cfg->rstar*cfg->rstar ) fprintf(fw,"## ");
	fprintf(fw,"% .8e % .8e % .8e % .8e % .8e % .8e %.12e %.8e %.8e %.8e %d\n",
	host.x[i],
        host.y[i],
	host.z[i],
	host.vx[i],
	host.vy[i],
	host.vz[i],
	host.N[i],
	host.m[i]*MSun/Grav,
	pow(3.0*host.m[i]*MSun/Grav/4.0/PI/cfg->rho,1.0/3.0),
	cfg->beta[host.id[i]],host.id[i]);
    }

 fclose(fw);
}

void Grid::write_sed(const char *ver){
    configdata *cfg = &cconfig;
    FILE   *fw;
    char   out[512];

    printsed(time);
    sprintf(out,"%s_SED_%.2f.dat",cfg->name_stub,time);
    fw=fopen(out,"wb");

    fprintf(fw,"#Calculated with DiskDyn Version %s\n",ver);
    fprintf(fw,"#Wav[mu]\tFlux_star[Jy]\tFlux_exc[Jy]\n");
    for ( int j=0 ; j<cfg->n_wav ; j++ ){ 
	fprintf(fw,"%.3f %.6e %.6e\n",cfg->w[j]*1e6,cfg->in[j]*1e26*cfg->w[j]*cfg->w[j]/cspeed,cfg->inexc[j]);
    }

    fclose(fw);
}

void Grid::massive_RK4(double dt){
    configdata *cfg = &cconfig;

    int n_grav = n_mass + cfg->n_planets;

    if ( n_grav >1 ){
	int TILE = 32;	// ( 8x8 - 16x16)
	int n_thread = 64;
	int block0 = int(( n_grav - 1.0) / n_thread + 1.0);

	int num      = int (( n_grav - 1.0 ) / TILE + 1 );	// 37500

	RK4calc_massive0<<<block0,n_thread>>>(n_grav,gpu);	// nulls Ks
	RK4calc_massive1<<<block0,n_thread>>>(n_grav,gpu);
	for ( int i=0 ; i < num ; i++ ){
	    int NLblocks = num - i;
	    RK4calc_massive2<<<NLblocks,TILE*TILE,2*TILE*(sizeof(struct array))+TILE*TILE*3*(sizeof(double))>>>(i,n_grav,TILE,gpu,cfg->n_planets);
	}
	RK4calc_massive3<<<block0,n_thread>>>(n_grav,gpu,dt*0.5);
	for ( int i=0 ; i < num ; i++ ){
	    int NLblocks = num - i;
	    RK4calc_massive4<<<NLblocks,TILE*TILE,2*TILE*(sizeof(struct array))+TILE*TILE*3*(sizeof(double))>>>(i,n_grav,TILE,gpu,cfg->n_planets);
	}
	RK4calc_massive5<<<block0,n_thread>>>(n_grav,gpu,dt*0.5);
	for ( int i=0 ; i < num ; i++ ){
	    int NLblocks = num - i;
	    RK4calc_massive6<<<NLblocks,TILE*TILE,2*TILE*(sizeof(struct array))+TILE*TILE*3*(sizeof(double))>>>(i,n_grav,TILE,gpu,cfg->n_planets);
	}
	RK4calc_massive7<<<block0,n_thread>>>(n_grav,gpu,dt);
	for ( int i=0 ; i < num ; i++ ){
	    int NLblocks = num - i;
	    RK4calc_massive8<<<NLblocks,TILE*TILE,2*TILE*(sizeof(struct array))+TILE*TILE*3*(sizeof(double))>>>(i,n_grav,TILE,gpu,cfg->n_planets);
	}
    }
}

void Grid::small_RK4(double dt){
    configdata *cfg = &cconfig;

    int n_grav = n_mass + cfg->n_planets;

    if ( n_dust >0 ){
	int block = int(( n_dust - 1.0) / 128.0 + 1.0);

        RK4calc_dust<<<block,128>>>(n_grav,n_dust,gpu,dt,cfg->FLQ,cfg->Wz,cfg->Bstar,cfg->rstar,cfg->rho,cfg->SWQ,cfg->vsw*year/AU,Rmin);
    }
}

// Overwrites original x,y,z,vx,vy,vz ... seperated off, because we need them for the small particle calcs first
void Grid::massive_final(double dt){
    configdata *cfg = &cconfig;

    int n_grav = n_mass + cfg->n_planets;

    if ( n_grav > 1 ){
	int block0 = int(( n_grav - 1.0) / 128.0 + 1.0);

	RK4calc_massive_final<<<block0,128>>>(n_grav,gpu,dt/6.0);
    }
}

void Grid::calculate_dust_temp(){
    configdata *cfg = &cconfig;

    double	star,gflux;

    Rmax = 0.0;
    for ( int i=0 ; i < cfg->n_disks   ; i++ ){ if ( cfg->disks[i].Rout > Rmax                    ) Rmax = cfg->disks[i].Rout; }
    for ( int i=0 ; i < cfg->n_blobs   ; i++ ){ if ( cfg->blobs[i].a*(1.0+cfg->blobs[i].ecc) + cfg->blobs[i].sigma > Rmax ) Rmax = cfg->blobs[i].a*(1.0+cfg->blobs[i].ecc) + cfg->blobs[i].sigma; }
    for ( int i=0 ; i < cfg->n_planets ; i++ ){ 
	double R = sqrt(cfg->planets[i].x*cfg->planets[i].x + cfg->planets[i].y*cfg->planets[i].y + cfg->planets[i].z*cfg->planets[i].z);
	if ( R > Rmax ) Rmax = R; 
    }
    Rmax *= 2*AU;	// Convert it to m;
    Rmin = Rmax;

    // Calculate Sublimation distance with smallest particle - adjust Rmin to that
    for ( int k=0 ; k<cfg->n_dust ; k++ ){
	double Rtmp = 0;

	star=0.0;
	for ( int i=1; i<cfg->n_wav ; i++ )
    	    star += 0.5*(cfg->w[i]-cfg->w[i-1])*(cfg->Qabs[k*cfg->n_wav+i-1]*cfg->in[i-1]+cfg->Qabs[k*cfg->n_wav+i]*cfg->in[i])/PI; 				/* Divide by pi, because model is B*pi */
	star *= pow(cfg->dist,2.0);

	gflux=0.0;
	for ( int i=1; i<cfg->n_wav ; i++ )
	    gflux += 0.5*(cfg->w[i]-cfg->w[i-1])*(cfg->Qabs[k*cfg->n_wav+i-1]*intens(Tsub,cfg->w[i-1])+cfg->Qabs[k*cfg->n_wav+i]*intens(Tsub,cfg->w[i])); 

	Rtmp = 0.5*sqrt(star/gflux);
	if ( Rtmp < Rmin ) Rmin = Rtmp;
    }

    // Now, calculate R-T array for all dust particles
    int nT = floor((2*Tsub-3)*0.5)+1;
    double *Rarr = (double *)malloc(sizeof(double)*nT);
    if (!Rarr) exit_with_usage(49);
    double *Tarr = (double *)malloc(sizeof(double)*nT);
    if (!Tarr) exit_with_usage(49);

    Th = (double *)malloc(sizeof(double)*nR*cfg->n_dust);
    if (!Th) exit_with_usage(49);

    for ( int k=0 ; k<cfg->n_dust ; k++ ){
	star=0.0;
	for ( int i=1; i<cfg->n_wav ; i++ )
	    star += 0.5*(cfg->w[i]-cfg->w[i-1])*(cfg->Qabs[k*cfg->n_wav+i-1]*cfg->in[i-1]+cfg->Qabs[k*cfg->n_wav+i]*cfg->in[i])/PI; 				/* Divide by pi, because model is B*pi */

	star *= pow(cfg->dist,2.0);


	int j=0;
	for ( double T = 2*Tsub ; T>=3.0 ; T-=2.0 ){
	    gflux=0.0;

	    for ( int i=1; i<cfg->n_wav ; i++ )
		gflux += 0.5*(cfg->w[i]-cfg->w[i-1])*(cfg->Qabs[k*cfg->n_wav+i-1]*intens(T,cfg->w[i-1])+cfg->Qabs[k*cfg->n_wav+i]*intens(T,cfg->w[i]));

	    Rarr[j] = 0.5*sqrt(star/gflux);
	    Tarr[j] = T;
	    j++;
	}
	// Now, interpolate

        int im=0;
        int ip=0;
	for ( int i=0 ; i<nR ; i++ ){
	    double R = Rmin * exp(log(Rmax/Rmin)*((double)i)/((double)nR-1.0));
	    while ( Rarr[ip] <= R ) ip++;
	    im = ip-1;
    	    Th[k*nR + i] = exp(log(Tarr[im])+log(Tarr[ip]/Tarr[im])*log(R/Rarr[im])/log(Rarr[ip]/Rarr[im]));
	    ip--;
	}
    }

    free(Tarr);
    free(Rarr);

    // Convert Rmin and Rmax back to AUs
    Rmin /= AU;
    Rmax /= AU;
}
