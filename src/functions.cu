#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ddyn.cuh"

void coe2rv(double M,double m,double a,double e,double i,double w,double o,\
	    double nu,double *x,double *y,double *z,double *vx,double *vy,double *vz,double beta){

    double p = a*(1.0-e*e);			// p = semi-latus rectum (semiparameter) a in AU
   
    double mu = (M+m)*(1.0-beta);		// Assuming masses are multiplied by Grav constant and are in Solar Units

    //  if circular and equatorial: 	w = 0, o = 0, nu = lambda_true
    //  if circular and inclined:   	w = 0, nu = u
    //  if elliptical and equatorial:	w = w_true, o = 0, 

    // Position Coordinates in Perifocal Coordinate System
    double xe  =  p*cos(nu)/(1.0+e*cos(nu)); 	// x-coordinate (AU)
    double ye  =  p*sin(nu)/(1.0+e*cos(nu)); 	// y-coordinate (AU)
    double ze  =  0.0;                          // z-coordinate (AU)
    double vxe = -sqrt(mu/p)*sin(nu);        	// velocity in x (AU/yr)
    double vye =  sqrt(mu/p)*(e+cos(nu));    	// velocity in y (AU/yr)
    double vze =  0;                            // velocity in z (AU/yr)

    // Rotate 

    *x  = xe*(cos(o)*cos(w)-sin(o)*sin(w)*cos(i))  + \
	  ye*(-cos(o)*sin(w)-sin(o)*cos(w)*cos(i)) + \
          ze*( sin(o)*sin(i));
    *y  = xe*(sin(o)*cos(w)+cos(o)*sin(w)*cos(i))  + \
	  ye*(-sin(o)*sin(w)+cos(o)*cos(w)*cos(i)) + \
          ze*(-cos(o)*sin(i));
    *z  = xe*sin(w)*sin(i) + \
	  ye*cos(w)*sin(i) + \
          ze*cos(i);

    *vx = vxe*(cos(o)*cos(w)-sin(o)*sin(w)*cos(i))  + \
	  vye*(-cos(o)*sin(w)-sin(o)*cos(w)*cos(i)) + \
	  vze*( sin(o)*sin(i));
    *vy = vxe*(sin(o)*cos(w)+cos(o)*sin(w)*cos(i))  + \
	  vye*(-sin(o)*sin(w)+cos(o)*cos(w)*cos(i)) + \
	  vze*(-cos(o)*sin(i));
    *vz = vxe*sin(w)*sin(i) + \
	  vye*cos(w)*sin(i) + \
	  vze*cos(i);
}

void rv2coe(double M,double m,double x,double y,double z,double vx,double vy,\
	    double vz,double *a,double *e,double *i,double *w,double *o,double *nu,double beta){
    double mu = (M+m)*(1.0-beta);		// Assuming masses are multiplied by Grav constant and are in Solar Units
						// a in AU, velocity in AU/year
    // magr and magv
    double  R = sqrt(x*x + y*y + z*z);
    double  V = sqrt(vx*vx + vy*vy + vz*vz);

    // hbar
    double hx = y*vz - z*vy;
    double hy = z*vx - x*vz;
    double hz = x*vy - y*vx;

    // hmag
    double H  = sqrt(hx*hx + hy*hy + hz*hz);

    if ( H > 1e-12){
	double nx = -hy;
	double ny =  hx;
    
	double  N = sqrt(nx*nx + ny*ny);
        double c1 = V*V - mu/R;
    
	double rv = x*vx + y*vy + z*vz;
	
	double ex = (c1*x - rv*vx)/mu;
	double ey = (c1*y - rv*vy)/mu;
	double ez = (c1*z - rv*vz)/mu;

	*e = sqrt(ex*ex+ey*ey+ez*ez);

        // ------------  find a e and semi-latus rectum   ----------
	double sme = V*V*0.5 - mu/R;
	if ( fabs(sme) > 1e-12 ) *a = -mu/2.0/sme;
	else                     *a = INFINITY;
    
        // -----------------  find inclination   -------------------
    
        double hk = hz/H;
        
	*i = acos(hk);

        // --------  determine type of orbit for later use  --------
        // ------ elliptical, parabolic, hyperbolic inclined -------
        int typeorbit = 0; // 0 - elliptical, inclined; 1 - circular equatorial, 2 - circular inclined; 3 - elliptical, equatorial
        if ( *e < 1e-12 ){
            // ----------------  circular equatorial ---------------
            if  ( *i < 1e-12 || fabs((*i)-PI) < 1e-12 )
                typeorbit = 1;
            else
            // --------------  circular inclined ---------------
                typeorbit = 2;
        }
        else{
            // - elliptical, parabolic, hyperbolic equatorial --
            if  ( *i < 1e-12 || fabs((*i)-PI) < 1e-12 )
                typeorbit = 3;
        }

        // ----------  find longitude of ascending node ------------
        if ( N > 1e-12 ){
            double temp = nx/N;
    	    // NaN failsafe
            if ( temp >  1.0 ) temp =  1.0;
    	    if ( temp < -1.0 ) temp = -1.0;
            *o = acos(temp);
            if ( ny < 0.0  ) *o = PI2 - (*o);
	}
	else{
            *o = 0.0; // random? Don't matter?
	}

        // ---------------- find argument of perigee ---------------
        if ( typeorbit == 0 ){
    	    if ( (*e)*N > 1e-24 ){
		double temp = (ex*nx + ey*ny)/(*e)/N;
	        // NaN failsafe
		if ( temp >  1.0 ) temp =  1.0;
	        if ( temp < -1.0 ) temp = -1.0;
		*w = acos(temp);
    	    }
    	    else *w = 0.0;	// will zero out - too close to zero

    	    if ( ez < 0.0  )       *w = PI2 - (*w);
	}
        else{ // Only inclined and elliptical orbit has it defined
            *w = 0.0;		// random? Don't matter?
	}

	// ------------  find true anomaly at epoch    -------------
	if ( typeorbit == 0 || typeorbit == 3 ){
    	    if ( (*e)*R > 1e-24 ){
		double temp = (ex*x + ey*y + ez*z)/(*e)/R;
		// NaN failsafe
		if ( temp >  1.0 ) temp =  1.0;
		if ( temp < -1.0 ) temp = -1.0;
		*nu = acos(temp);
    	    }
    	    else *nu = 0.0;

    	    if ( rv < 0.0  ) *nu = PI2 - (*nu);
	}
	else{ // Circular orbit - true anomaly not defined
            *nu = 0.0;		// random? Don't matter?
	}
    }
}

// Generate a random particle orbit within the disk
void genorb_disk(int j,double m,long *rseed,double *x,double *y,double *z,\
		 double *vx,double *vy,double *vz,int id){
    configdata *cfg = &cconfig;
    double a,e,i,w,o,nu,beta;

    long       seed = *rseed;
    
    do{
	e = sqrt(-2.0*log(ran3(&seed)))*cos(PI2*ran3(&seed))*cfg->disks[j].sigmae;
	if ( e<0 ) e = -e;
    } while( e>=1.0 );

    double sem_a_min = cfg->disks[j].Rin  / (1.0 - e );
    double sem_a_max = cfg->disks[j].Rout / (1.0 + e );

    a = pow(ran3(&seed)*(pow(sem_a_max,2.0+cfg->disks[j].gammas)-pow(sem_a_min,2.0+cfg->disks[j].gammas))+\
	pow(sem_a_min,2.0+cfg->disks[j].gammas),1.0/(2.0+cfg->disks[j].gammas));

    i = sqrt(-2.0*log(ran3(&seed)))*cos(PI2*ran3(&seed))*cfg->disks[j].dhdr;
    if ( i < 0) i = -i;

    nu = PI2 * ran3(&seed);
    w  = PI2 * ran3(&seed);
    o  = PI2 * ran3(&seed);


    if ( id < 0 ) beta = 0.0;
    else          beta = cfg->beta[id];

    if ( beta >= 0.5 ) beta = 0.0;	// It will get removed anyway, let's get it out quick!

    coe2rv(cfg->planets[0].m,m,a,e,i,w,o,nu,x,y,z,vx,vy,vz,beta);

    *rseed=seed;
}

// Generate a random particle orbit for a blob
void genorb_blob(int j,double m,long *rseed,double *x,double *y,double *z,double *vx,double *vy,double *vz,int id){
    configdata *cfg = &cconfig;
    double a,e,i,w,o,nu,beta;

    long       seed = *rseed;

    do{
        a = cfg->blobs[j].a+sqrt(-2.0*log(ran3(&seed)))*cos(PI2*ran3(&seed))*cfg->blobs[j].sigma;
    } while( a<=0.0 );
    e  = cfg->blobs[j].ecc;
    i  = cfg->blobs[j].inc;
    w  = cfg->blobs[j].w;
    o  = cfg->blobs[j].o;
    nu = cfg->blobs[j].nu + sqrt(-2.0*log(ran3(&seed)))*cos(PI2*ran3(&seed))*cfg->blobs[j].sigma/a;

    if ( id < 0 ) beta = 0.0;
    else          beta = cfg->beta[id];

    if ( beta >= 0.5 ) beta = 0.0;	// It will get removed anyway, let's get it out quick!

    coe2rv(cfg->planets[0].m,m,a,e,i,w,o,nu,x,y,z,vx,vy,vz,beta);

    // Collisional products - I'll be repeating a lot of the calculations here. Not too quick.
    if ( cfg->blobs[j].collcloud==0 ){

	double voff = cfg->blobs[j].vkepp*0.01*sqrt((*vx)*(*vx)+(*vy)*(*vy)+(*vz)*(*vz));
	double costheta = -1.0 + 2.0*ran3(&seed);
	double sintheta = sqrt(1.0-costheta*costheta);
	double gamma    = PI2*ran3(&seed);
	double cosgamma = cos(gamma);
	double singamma = sin(gamma);
	
	*vx += voff * sintheta*cosgamma;
	*vy += voff * sintheta*singamma;
	*vz += voff * costheta;
    }

    *rseed=seed;
}

// Integrate up the size distribution function until we get a complete N=1
void getbod(long double *mmax,double *m,double C,double ome){
    long double mmin = pow(pow(*mmax,ome)-ome/C,(1.0/ome));

    *m = C*(pow(*mmax,1.0+ome)-pow(mmin,1.0+ome))/(1.0+ome)/MSun*Grav;
    *mmax = mmin;
}

// Interpolates the dust optical properties to the grid defined in the calculations
void interpolate_dust_to_grid(){
    configdata *cfg = &cconfig;
    double *a,*Qabs,*Qsca,*Qpfunc;
    double f = log(cfg->amax/cfg->amin)/((double)cfg->n_grid-1.0);

    a      = (double *)malloc(sizeof(double)*cfg->n_grid);
    if ( !a )		exit_with_usage(50);
    cfg->mem += sizeof(double)*cfg->n_grid;
    Qabs   = (double *)malloc(sizeof(double)*cfg->n_grid*cfg->n_wav);
    if ( !Qabs )	exit_with_usage(50);
    cfg->mem += sizeof(double)*cfg->n_grid*cfg->n_wav;
    Qsca   = (double *)malloc(sizeof(double)*cfg->n_grid*cfg->n_wav);
    if ( !Qsca )	exit_with_usage(50);
    cfg->mem += sizeof(double)*cfg->n_grid*cfg->n_wav;
    Qpfunc = (double *)malloc(sizeof(double)*cfg->n_grid*cfg->n_wav*cfg->n_theta);
    if ( !Qpfunc )	exit_with_usage(50);
    cfg->mem += sizeof(double)*cfg->n_grid*cfg->n_wav*cfg->n_theta;
    cfg->q = (float *)malloc(sizeof(float)*cfg->n_grid);
    if ( !cfg->q )	exit_with_usage(50);
    cfg->mem += sizeof(float)*cfg->n_grid;

    int im=0;
    int ip=0;
    for ( int i=0 ; i<cfg->n_grid ; i++){
	a[i] = cfg->amin * exp((double)i*f);
	cfg->q[i] = cfg->kqe * a[i]*1e9;		// Could store this! cbrt takes a long time
	// Numerical safety switches!
	if ( i==cfg->n_grid-1 && (a[i]-cfg->dmax) > 1e-12 ){ 
	    printf("i: %d ngrid: %d a[i]: %.8e dmax: %.8e diff: %.8e\nLarger? Shouldn't happen!\n",i,cfg->n_grid,a[i],cfg->dmax,a[i]-cfg->dmax); 
	    exit_with_usage(54);
	}
	// Fail switch
	if ( i==cfg->n_grid-1 && (a[i]-cfg->dmax) < 1e-12 && (a[i]-cfg->dmax) > 0 ) a[i] = cfg->dmax;
	while ( cfg->s[ip] < a[i] ) ip++;
	im = ip-1;
	for ( int j=0 ; j<cfg->n_wav ; j++ ){
	    Qabs[i*cfg->n_wav+j] = exp(log(cfg->Qabs[im*cfg->n_wav+j])+log(a[i]/cfg->s[im])*log(cfg->Qabs[ip*cfg->n_wav+j]/cfg->Qabs[im*cfg->n_wav+j])/log(cfg->s[ip]/cfg->s[im]));
	    Qsca[i*cfg->n_wav+j] = exp(log(cfg->Qsca[im*cfg->n_wav+j])+log(a[i]/cfg->s[im])*log(cfg->Qsca[ip*cfg->n_wav+j]/cfg->Qsca[im*cfg->n_wav+j])/log(cfg->s[ip]/cfg->s[im]));
//	    if ( Qabs[i*cfg->n_wav+j] != Qabs[i*cfg->n_wav+j] )
//		printf("a: %.6e wav: %.6e Qam %.6e ai %.6e s_im %.6e Qap %.6e s_ip %.6e \n",a[i],cfg->w[j],cfg->Qabs[im*cfg->n_wav+j],a[i],cfg->s[im],cfg->Qabs[ip*cfg->n_wav+j],cfg->s[ip]);
	    for ( int k=0 ; k<cfg->n_theta ; k++ ){
		Qpfunc[cfg->n_wav*(cfg->n_theta*i+k)+j] = exp(log(cfg->Qpfunc[cfg->n_wav*(cfg->n_theta*im+k)+j])+log(a[i]/cfg->s[im])*log(cfg->Qpfunc[cfg->n_wav*(cfg->n_theta*ip+k)+j]/cfg->Qpfunc[cfg->n_wav*(cfg->n_theta*im+k)+j])/log(cfg->s[ip]/cfg->s[im]));
	    }
	}
	ip--;
    }

    free(cfg->s);
    free(cfg->Qabs);
    free(cfg->Qsca);
    free(cfg->Qpfunc);

    cfg->s=a;
    cfg->Qabs=Qabs;
    cfg->Qsca=Qsca;
    cfg->Qpfunc=Qpfunc;
    cfg->n_dust=cfg->n_grid;
}

void interpolate_stellar_to_grid(){
    configdata *cfg = &cconfig;

    cfg->in    = (double *)malloc(sizeof(double)*cfg->n_wav);
    if ( !cfg->in )	exit_with_usage(51);
    cfg->mem += sizeof(double)*cfg->n_wav;
    cfg->inexc = (double *)malloc(sizeof(double)*cfg->n_wav);
    if ( !cfg->inexc )	exit_with_usage(51);
    cfg->mem += sizeof(double)*cfg->n_wav;

    int im=0;
    int ip=0;
    for ( int j=0 ; j<cfg->n_wav ; j++ ){
	if ( cfg->w[j] < cfg->spec[0].wvl )
	    cfg->in[j] = 0.0;
	else if ( cfg->w[j] > cfg->spec[cfg->nspec-1].wvl )
	    cfg->in[j] = cfg->spec[cfg->nspec-1].F*pow(cfg->spec[cfg->nspec-1].wvl/cfg->w[j],4.0);
	else{
	    while ( cfg->spec[ip].wvl <= cfg->w[j] ) ip++;
	    im = ip-1;
	    if ( cfg->spec[im].F!=0 && cfg->spec[ip].F!=0 )
		cfg->in[j] = exp(log(cfg->spec[im].F)+log(cfg->w[j]/cfg->spec[im].wvl)*log(cfg->spec[ip].F/cfg->spec[im].F)/log(cfg->spec[ip].wvl/cfg->spec[im].wvl));
	    else
		cfg->in[j] = 0.0;
	    ip--;
	}	    
    }

    for ( int j=0 ; j<cfg->nwavs ; j++ ){
	int i=0;
	if ( cfg->wavs[j] < cfg->w[0] || cfg->wavs[j] > cfg->w[cfg->n_wav-1] ) exit_with_usage(43);
	while ( cfg->w[i] < cfg->wavs[j] && i<cfg->n_wav-2 ) i++;
	cfg->wj[j]=i-1;	
    }
}

void calculate_beta(){
    configdata *cfg = &cconfig;

    cfg->beta   = (float *)malloc(sizeof(float)*cfg->n_dust);
    if ( !cfg->beta )	exit_with_usage(52);
    cfg->betasw = (float *)malloc(sizeof(float)*cfg->n_dust);
    if ( !cfg->betasw )	exit_with_usage(52);
    cfg->mem += 2*sizeof(float)*cfg->n_dust;
    double *Qpr = (double *)malloc(sizeof(double)*cfg->n_wav);
    if ( !Qpr )		exit_with_usage(52);

    cfg->betaswconst = 3.0 * cfg->mloss * Cd * MSun * year * cfg->vsw / 32.0 / PI / cfg->planets[0].m / cfg->rho / AU / AU / AU;

    // Integrate star once
    double Star=0.0;
    for ( int j=0 ; j<cfg->n_wav-1 ; j++ ) Star += 0.5*(cfg->w[j+1]-cfg->w[j])*(cfg->in[j+1]+cfg->in[j]);

    cfg->Lum = Star*4.0*PI*pow(cfg->dist,2.0)/3.827e26;

    // Integrate Qpfunc for anisotropy parameter for all particles, all wavelengths
    for ( int i=0 ; i<cfg->n_dust ; i++ ){
	if ( cfg->radQ==1 ){
	    double Qsum = 0.0;
	    for ( int j=0 ; j<cfg->n_wav ; j++ ){
		double g=0.0;
		double S=0.0;
		double costm = -1.0;
		for ( int k=1 ; k<cfg->n_theta ; k++ ){
		    double costp = (double)k*2.0/(cfg->n_theta-1.0)-1.0;	
		    S += 0.5*(cfg->Qpfunc[cfg->n_wav*(cfg->n_theta*i+k-1)+j]+cfg->Qpfunc[cfg->n_wav*(cfg->n_theta*i+k)+j])*(costp-costm);
		    g += 0.5*(costm*cfg->Qpfunc[cfg->n_wav*(cfg->n_theta*i+k-1)+j]+costp*cfg->Qpfunc[cfg->n_wav*(cfg->n_theta*i+k)+j])*(costp-costm);
		    costm=costp;
		}
		g /= S;
		Qpr[j] = (cfg->Qabs[i*cfg->n_wav+j]+cfg->Qsca[i*cfg->n_wav+j]*(1.0-g))*cfg->in[j];
	    }
	    for ( int j=0 ; j<cfg->n_wav-1 ; j++ ) Qsum += 0.5*(cfg->w[j+1]-cfg->w[j])*(Qpr[j+1]+Qpr[j]);
	    cfg->beta[i] = 0.57 * Qsum/Star * (cfg->Lum/cfg->planets[0].m*Grav) / cfg->s[i]/1e6 / cfg->rho/1e-3;
	}
	else{
	    cfg->beta[i]  = 0.0;
	}
	if ( cfg->SWQ==1 ){
	    cfg->betasw[i] = cfg->betaswconst / cfg->s[i];
	}
	else{
	    cfg->betasw[i] = 0.0;
	}
    }
    free(Qpr);
}

double intens(double temp,double wvl) /* B_lambda */
{
 double 	b;

    b=(2*hplank*cspeed*cspeed)/(pow(wvl,5)*(exp(hplank*cspeed/(wvl*kboltz*temp))-1));

 return(b);
}

void checkmemusage(int n_all,int n_dust){
    configdata *cfg = &cconfig;

    cfg->memgpu     = sizeof(double)*41*n_all+sizeof(int)*n_all+sizeof(double)*(cfg->n_dust*(cfg->n_wav*(cfg->n_theta+2)+nR)+2*cfg->n_wav+cfg->xsize*cfg->ysize+n_dust);
    cfg->memgpudust = sizeof(double)*(cfg->n_dust*(cfg->n_wav*(cfg->n_theta+2)+nR)+2*cfg->n_wav);
}
