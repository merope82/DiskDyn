/******************************************************************************/
/* Configuration file processing program                                      */
/* Much thanks to Andras Pal for source codes on character processing         */
/******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ddyn.cuh"

int char_is_space(int c)
{
 if ( c==64 || c==32 || c==13 || c==10 || c==9 )        return(1);
 else                                   return(0);
}

void remove_quotes(char *buff)
{
 int k;
 while ( *buff )
  {     for ( k=0 ; buff[k]=='"' ; )    k++;
        if ( k )        memmove(buff,buff+k,strlen(buff)+1-k);
        else            buff++;
  }
}

int tokenize_spaces(char *buff,char **tokens,int max)
{
 int    intoken,inquota,n;
 char   **tsave;

 tsave=tokens;

 intoken=0,inquota=0;n=0;
 while ( *buff && n<max )
  {     if ( ( ! char_is_space(*buff) ) && ! intoken )
         {      *tokens=buff;
                intoken=!0,inquota=0;n++;
                if ( *buff=='"' )       inquota=!0;
                tokens++,buff++;
         }
        else if ( intoken && ( (char_is_space(*buff) && inquota) || (!char_is_space(*buff)) ) )
         {      if ( *buff=='"' )       inquota=!inquota;
                buff++;
         }
        else if ( intoken && ! inquota && char_is_space(*buff) )
         {      *buff=0,buff++;
                intoken=0;
         }
        else    buff++;
  };
 *tokens=NULL;

 while ( *tsave != NULL )
  {     remove_quotes(*tsave);
        tsave++;
  };

 return(n);
}

int add_star(configdata *cfg,double mass){ // Star is planet #0.
 
 cfg->n_planets++;
 cfg->planets    = (planet_cart *)malloc(sizeof(planet_cart));
 if ( cfg->planets == NULL ) exit_with_usage(44);
 cfg->mem += sizeof(planet_cart);

 cfg->planets[0].x      = 0.0;
 cfg->planets[0].y      = 0.0;
 cfg->planets[0].z      = 0.0;
 cfg->planets[0].vx     = 0.0;
 cfg->planets[0].vy     = 0.0;
 cfg->planets[0].vz     = 0.0;
 cfg->planets[0].m      = mass * Grav;	// Stellar mass in Solar Units * Grav constant

 return(0);
}

int add_planets(configdata *cfg,char **dat){
 int	j;
 double	i,a,e,nu,o,w,m;
 double x,y,z,vx,vy,vz;

 sscanf(dat[0],"%lg",&i);
 sscanf(dat[1],"%lg",&a);
 sscanf(dat[2],"%lg",&e);
 sscanf(dat[3],"%lg",&nu);
 sscanf(dat[4],"%lg",&o);
 sscanf(dat[5],"%lg",&w);
 sscanf(dat[6],"%lg",&m);

 j=cfg->n_planets; cfg->n_planets++;
 cfg->planets    = (planet_cart *)realloc(cfg->planets,    sizeof(planet_cart)*(j+1));
 cfg->mem += sizeof(planet_cart);

 i  *= M_PI/180.0;
 nu *= M_PI/180.0;
 o  *= M_PI/180.0;
 w  *= M_PI/180.0;
 m  *= Grav;

 coe2rv(cfg->planets[0].m,m,a,e,i,w,o,nu,&x,&y,&z,&vx,&vy,&vz,0.0);

 cfg->planets[j].x  = x;
 cfg->planets[j].y  = y;
 cfg->planets[j].z  = z;
 cfg->planets[j].vx = vx;
 cfg->planets[j].vy = vy;
 cfg->planets[j].vz = vz;
 cfg->planets[j].m  = m;

 return(0);
}

int add_disks(configdata *cfg,char **dat){		// review
 double	  Rin,Rout,dhdr,sigmae,amax,amin,etaa,mtot,gammas;
 int	  j;
 long int Ngrav,Nsmall;

 sscanf(dat[0], "%lg",&Rin);
 sscanf(dat[1], "%lg",&Rout);
 sscanf(dat[2], "%lg",&dhdr);
 sscanf(dat[3], "%lg",&sigmae);
 sscanf(dat[4], "%lg",&amax);
 sscanf(dat[5], "%lg",&amin);
 sscanf(dat[6], "%lg",&etaa);
 sscanf(dat[7], "%lg",&mtot);
 sscanf(dat[8], "%lg",&gammas);
 sscanf(dat[9], "%ld" ,&Ngrav);
 sscanf(dat[10],"%ld" ,&Nsmall);

 j=cfg->n_disks;
 if ( j==0 ){ 
    cfg->disks=(diskconf *)malloc(sizeof(diskconf));
    if ( !cfg->disks ) exit_with_usage(45);
    cfg->mem += sizeof(diskconf);
 }
 else{ 
    diskconf *tmp;
    tmp = (diskconf *)realloc(cfg->disks,sizeof(diskconf)*(j+1)); 
    if (tmp)	cfg->disks = tmp;
    else	exit_with_usage(45);
    cfg->mem += sizeof(diskconf);
 }

 cfg->disks[j].Rin    = Rin;			// Rin and Rout stays in AU
 cfg->disks[j].Rout   = Rout;			// Rin and Rout stays in AU
 cfg->disks[j].dhdr   = atan(dhdr);
 cfg->disks[j].sigmae = sigmae;
 cfg->disks[j].amax   = amax*1000.0;
 cfg->disks[j].amin   = amin/1.0e6;
 cfg->disks[j].eta    = (etaa+2.0)/3.0;
 cfg->disks[j].Mtot   = mtot*MEarth;		// in kg
 cfg->disks[j].mmax   = 4.0 * pow(amax*1000.0,3.0) * PI * cfg->rho / 3.0;	// in kg
 cfg->disks[j].mmin   = 4.0 * pow(amin/1.0e6,3.0) * PI * cfg->rho / 3.0;	// in kg
 cfg->disks[j].gammas = gammas;
 cfg->disks[j].Ngrav  = Ngrav;
 cfg->disks[j].Nsmall = Nsmall;

 cfg->n_disks++;
 return(0);
}

int add_blobs(configdata *cfg,char **dat,int t){	// review
 double	a,ecc,inc,w,o,nu,sigma,amax,amin,etaa,mtot,vkepp;
 int	j,collcloud;
 long	Ngrav,Nsmall;
 
 if ( t<14 ) exit_with_usage(58);

 sscanf(dat[0], "%lg", &a);
 sscanf(dat[1], "%lg", &ecc);
 sscanf(dat[2], "%lg", &inc);
 sscanf(dat[3], "%lg", &w);
 sscanf(dat[4], "%lg", &o);
 sscanf(dat[5], "%lg", &nu);
 sscanf(dat[6], "%lg", &amax);
 sscanf(dat[7], "%lg", &amin);
 sscanf(dat[8], "%lg", &etaa);
 sscanf(dat[9], "%lg", &mtot);
 sscanf(dat[10],"%ld" ,&Ngrav);
 sscanf(dat[11],"%ld" ,&Nsmall);
 if      ( strcmp(dat[12],"cloud"     )==0  || 
           strcmp(dat[12],"Cloud"     )==0 ){
	   collcloud = 1;
	   sscanf(dat[13],"%lg",&sigma);
	   vkepp = 0.0;
 }
 else if ( strcmp(dat[12],"collision" )==0  || 
	   strcmp(dat[12],"Collision" )==0 ){
	   collcloud = 0;
	   sscanf(dat[13],"%lg",&vkepp);
	   sigma = 0.0;
 }
 else						exit_with_usage(57);

 j=cfg->n_blobs; 
 if ( j==0 ){ 
    cfg->blobs=(blobconf *)malloc(sizeof(blobconf));
    if ( !cfg->blobs ) exit_with_usage(46);
    cfg->mem += sizeof(blobconf);
 }
 else{ 
    blobconf *tmp;
    tmp = (blobconf *)realloc(cfg->blobs,sizeof(blobconf)*(j+1));
    if (tmp)	cfg->blobs = tmp;
    else	exit_with_usage(46);
    cfg->mem += sizeof(blobconf);
 }

 cfg->blobs[j].a         = a;
 cfg->blobs[j].ecc       = ecc;
 cfg->blobs[j].inc       = inc*M_PI/180.0;
 cfg->blobs[j].w         = w *M_PI/180.0;
 cfg->blobs[j].o         = o*M_PI/180.0;
 cfg->blobs[j].nu        = nu*M_PI/180.0;
 cfg->blobs[j].sigma     = sigma;
 cfg->blobs[j].amax      = amax*1000.0;
 cfg->blobs[j].amin      = amin/1.0e6;
 cfg->blobs[j].eta       = (etaa+2.0)/3.0;
 cfg->blobs[j].Mtot      = mtot*MEarth;
 cfg->blobs[j].mmax      = 4.0 * pow(amax*1000.0,3.0) * PI * cfg->rho / 3.0;	// in kg
 cfg->blobs[j].mmin      = 4.0 * pow(amin/1.0e6,3.0) * PI * cfg->rho / 3.0;	// in kg
 cfg->blobs[j].Ngrav     = Ngrav;
 cfg->blobs[j].Nsmall    = Nsmall;
 cfg->blobs[j].collcloud = collcloud;
 cfg->blobs[j].vkepp     = vkepp;

 cfg->n_blobs++;

 return(0);
}

void confCdisks(){
 configdata  *cfg=&cconfig;

 // reconfiguring disks to same mass grid
 for ( int d=0 ; d<cfg->n_disks ; d++){

    cfg->disks[d].C = cfg->disks[d].Mtot;

    // N = C * r^(gammas+1) * exp(-theta^2/2.0/atan(dhdr)^2) * m^(-eta)

    if      ( cfg->disks[d].eta==1 ) cfg->disks[d].C /= (cfg->disks[d].mmax - cfg->disks[d].mmin);
    else if ( cfg->disks[d].eta==2 ) cfg->disks[d].C /= log(cfg->disks[d].mmax/cfg->disks[d].mmin);
    else                             cfg->disks[d].C /= (pow(cfg->disks[d].mmax,2.0-cfg->disks[d].eta)-pow(cfg->disks[d].mmin,2.0-cfg->disks[d].eta))/(2.0-cfg->disks[d].eta); // Division does not cancel when integrating for number (integral defined by mass)
 }

}

void confCblobs(){
 configdata  *cfg=&cconfig;

// reconfiguring blobs to same mass grid
 for ( int d=0 ; d<cfg->n_blobs ; d++){

    cfg->blobs[d].C = cfg->blobs[d].Mtot;

    if      ( cfg->blobs[d].eta==1 )     cfg->blobs[d].C /= (cfg->blobs[d].mmax - cfg->blobs[d].mmin);
    else if ( cfg->blobs[d].eta==2 )     cfg->blobs[d].C /= log(cfg->blobs[d].mmax/cfg->blobs[d].mmin);
    else                                 cfg->blobs[d].C /= (pow(cfg->blobs[d].mmax,2.0-cfg->blobs[d].eta)-pow(cfg->blobs[d].mmin,2.0-cfg->blobs[d].eta))/(2.0-cfg->blobs[d].eta);
 }

}

void detminmax(){		// minimum and maximum dust size! Don't care about gravitating bodies here!
 configdata  *cfg=&cconfig;
 double amin,amax,kappa;

 kappa = 4.0 * PI * cfg->rho / 3.0;

 // Set minimum to optical minimum if it is larger than the disk minimum
 
 for ( int i=0 ; i < cfg->n_disks ; i++ ){
    if ( cfg->disks[i].amin < cfg->dmin ){
	print_reset_min(i,0);
	cfg->disks[i].amin = cfg->dmin;
	cfg->disks[i].mmin = kappa * pow(cfg->disks[i].amin,3.0);
    }
    if ( cfg->disks[i].amax > cfg->dmax ){
	print_reset_max(i,0);
	cfg->disks[i].dmax = cfg->dmax;
    }
    else
	cfg->disks[i].dmax = cfg->disks[i].amax;

    cfg->disks[i].dmmax = kappa * pow(cfg->disks[i].dmax,3.0);
 }

 for ( int i=0 ; i < cfg->n_blobs ; i++ ){
    if ( cfg->blobs[i].amin < cfg->dmin ){
	print_reset_min(i,1);
	cfg->blobs[i].amin = cfg->dmin;
	cfg->blobs[i].mmin = kappa * pow(cfg->blobs[i].amin,3.0);
    }
    if ( cfg->blobs[i].amax > cfg->dmax ){
	print_reset_max(i,1);
	cfg->blobs[i].dmax = cfg->dmax;
    }
    else
	cfg->blobs[i].dmax  = cfg->blobs[i].amax;

    cfg->blobs[i].dmmax = kappa * pow(cfg->blobs[i].dmax,3.0);
 }

 // Set global minimum and maximum value

 if ( cfg->n_disks > 0 ){
    amin = cfg->disks[0].amin;
    amax = cfg->disks[0].dmax;
 }
 else if ( cfg->n_disks <= 0 && cfg->n_blobs > 0 ){
    amin = cfg->blobs[0].amin;
    amax = cfg->blobs[0].dmax;
 }
 else{
    amin = cfg->dmin;
    amax = cfg->dmax;
 }

 for ( int i=1 ; i < cfg->n_disks ; i++ ){
    if ( cfg->disks[i].amin < amin ) amin = cfg->disks[i].amin;
    if ( cfg->disks[i].dmax > amax ) amax = cfg->disks[i].dmax;
 }

 for ( int i=0 ; i < cfg->n_blobs ; i++ ){
    if ( cfg->blobs[i].amin < amin ) amin = cfg->blobs[i].amin;
    if ( cfg->blobs[i].dmax > amax ) amax = cfg->blobs[i].dmax;
 }

 // Set global values

 cfg->amin   = amin;
 cfg->amax   = amax;
 cfg->mmin   = kappa * pow(amin,3.0);
 cfg->mmax   = kappa * pow(amax,3.0);					// cfg mmax is for dust
 cfg->n_grid = (int)(log(cfg->mmax/cfg->mmin)/log(mgrid)+1.0);
 cfg->f      = log(cfg->mmax/cfg->mmin)/((double)cfg->n_grid-1.0);

}

int recalc_COM(configdata *cfg){
 int j;
 double Xc,Yc,Zc,Vx,Vy,Vz,Mtot;

 if ( cfg->n_planets == 1 ) return(0);

 Xc=Yc=Zc=Mtot=0;

 for ( j=0 ; j<cfg->n_planets ; j++ ){
    Xc   += cfg->planets[j].m * cfg->planets[j].x;
    Yc   += cfg->planets[j].m * cfg->planets[j].y;
    Zc   += cfg->planets[j].m * cfg->planets[j].z;
    Mtot += cfg->planets[j].m;
 }
 Xc/=Mtot;Yc/=Mtot;Zc/=Mtot;
 for ( j=0 ; j<cfg->n_planets ; j++ ){
    cfg->planets[j].x -= Xc;
    cfg->planets[j].y -= Yc;
    cfg->planets[j].z -= Zc;
 }
 Vx=Vy=Vz=0;
 for ( j=1 ; j<cfg->n_planets ; j++ ){
    Vx -= cfg->planets[j].m * cfg->planets[j].vx;
    Vy -= cfg->planets[j].m * cfg->planets[j].vy;
    Vz -= cfg->planets[j].m * cfg->planets[j].vz;
 }
 cfg->planets[0].vx = Vx/cfg->planets[0].m ;
 cfg->planets[0].vy = Vy/cfg->planets[0].m ;
 cfg->planets[0].vz = Vz/cfg->planets[0].m ;

 return(0);
}

void add_specfile(){
 FILE        *fr;
 int         i,t,n;
 double	     w,F;
 char        buff[512],*dat[3];
 configdata  *cfg=&cconfig;

 fr=fopen(cfg->specfile,"rb");
 if ( fr==NULL )        exit_with_usage(39);

 cfg->spec = (specdata *)malloc(sizeof(specdata));
 if ( !cfg->spec ) exit_with_usage(47);
 cfg->mem += sizeof(specdata);
 
 i=0;
 while ( ! feof(fr) ){
    if (fgets(buff,255,fr)==NULL)
        break;
    if ( buff[0]=='#' )
        continue;
    n=tokenize_spaces(buff,dat,2);
    if ( n<2 )
        continue;
    t=0;
    t+=sscanf(dat[0],"%lf",&w);
    t+=sscanf(dat[1],"%lf",&F);
    if ( t<2 )
	continue;
    specdata *tmp;
    tmp = (specdata *)realloc(cfg->spec,sizeof(specdata)*(i+1));
    if (tmp)	cfg->spec = tmp;
    else	exit_with_usage(47);
    cfg->mem += sizeof(specdata);
    cfg->spec[i].wvl = w;
    cfg->spec[i].F   = F*pow(cfg->rstar*AU,2)/pow(cfg->dist,2);
    i++;
 }
 if ( fileno(fr) != fileno(stdin) )	fclose(fr);

 cfg->nspec=i;
}

void add_optical_data(){
 FILE        *fr;
 int         i,t,itemp,j,nw,l;
 double	     dtemp;
 char        buff[256],*dat[13];
 configdata  *cfg=&cconfig;

 fr=fopen(cfg->comp,"rb");
 if ( fr==NULL )        exit_with_usage(39);

 i=j=t=nw=l=0;
 while ( ! feof(fr) ){
    if (fgets(buff,255,fr)==NULL)
        break;
    if (i==0){
	i++;
	continue;
    }
    t=tokenize_spaces(buff,dat,12);
    if ( i==1 ){	// Definition line for malloc
	sscanf(dat[1],"%d",&itemp);
	cfg->n_wav = itemp;
	sscanf(dat[2],"%d",&itemp);
	cfg->n_theta = itemp;
	sscanf(dat[3],"%d",&itemp);
	cfg->n_dust = itemp;
	sscanf(dat[4],"%lf",&dtemp);
	cfg->dmin = dtemp*1e-6;
	sscanf(dat[5],"%lf",&dtemp);
	cfg->dmax = dtemp*1e-6;
//	printf("nwav: %d ntheta: %d n_dust: %d\n",cfg->n_wav,cfg->n_theta,cfg->n_dust);
	cfg->w      = (double *)malloc(sizeof(double)*cfg->n_wav);
	if ( !cfg->w )		exit_with_usage(48);
	cfg->s      = (double *)malloc(sizeof(double)*cfg->n_dust);
	if ( !cfg->s )		exit_with_usage(48);
	cfg->Qabs   = (double *)malloc(sizeof(double)*cfg->n_wav*cfg->n_dust);
	if ( !cfg->Qabs )	exit_with_usage(48);
	cfg->Qsca   = (double *)malloc(sizeof(double)*cfg->n_wav*cfg->n_dust);
	if ( !cfg->Qsca )	exit_with_usage(48);
	cfg->Qpfunc = (double *)malloc(sizeof(double)*cfg->n_wav*cfg->n_dust*cfg->n_theta);
	if ( !cfg->Qpfunc )	exit_with_usage(48);
    }
    if ( i>1 && buff[0]=='#' ){
	sscanf(dat[1],"%lf",&dtemp);
	cfg->s[j] = dtemp*1e-6;
	print_status(j,cfg->n_dust-1);
	j++; nw=0; l=0;
    }
    if ( i>1 && buff[0]!='#' && t==3 && j==1 ){	// First set, let's gather wavelengths
	sscanf(dat[0],"%lf",&dtemp);
	cfg->w[nw] = dtemp * 1e-6;
    }
    if ( i>1 && buff[0]!='#' && t==3 ){
	sscanf(dat[1],"%lf",&dtemp);
	cfg->Qabs[cfg->n_wav*(j-1)+nw] = dtemp;
	sscanf(dat[2],"%lf",&dtemp);
	cfg->Qsca[cfg->n_wav*(j-1)+nw] = dtemp;
	nw++;
    }
    if ( i>1 && buff[0]!='#' && t==10 ){
	sscanf(dat[0],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+0+l*10] = dtemp;
	sscanf(dat[1],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+1+l*10] = dtemp;
	sscanf(dat[2],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+2+l*10] = dtemp;
	sscanf(dat[3],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+3+l*10] = dtemp;
	sscanf(dat[4],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+4+l*10] = dtemp;
	sscanf(dat[5],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+5+l*10] = dtemp;
	sscanf(dat[6],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+6+l*10] = dtemp;
	sscanf(dat[7],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+7+l*10] = dtemp;
	sscanf(dat[8],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+8+l*10] = dtemp;
	sscanf(dat[9],"%lf",&dtemp);
	cfg->Qpfunc[cfg->n_theta*cfg->n_wav*(j-1)+9+l*10] = dtemp;
	l++;
    }
    i++;
 }
 if ( fileno(fr) != fileno(stdin) )	fclose(fr);
 cfg->prev=0;

}

int read_in_param(char *param){
 FILE        *fr;
 int         i,t,itemp;
 double	     dtemp;
 char        buff[512],*dat[64];
 configdata  *cfg=&cconfig;

 cfg->n_planets=cfg->n_disks=cfg->n_blobs=0;

 fr=fopen(param,"rb");
 if ( fr==NULL )        exit_with_usage(3);

 cfg->ram  = getMemorySize();
 cfg->mem += sizeof(configdata);
 cfg->prev = 0;

 i=0;
 while ( ! feof(fr) ){
    if ( fgets(buff,512,fr)==NULL )     break;
    if ( buff[0]=='#' )                 continue;
    t=tokenize_spaces(buff,dat,16);
    if ( t==0 )				continue;
    if ( i==0 ){											// verbose
	sscanf(dat[0],"%d",&itemp);
	if ( itemp <0 ) exit_with_usage(4);
	else cfg->verb = itemp;	    
    }
    if ( i==1 ){											// GPU ID
	sscanf(dat[0],"%d",&itemp);
	if ( itemp <0 ) exit_with_usage(6);
	else cfg->gpuid = itemp;	    
    }
    if ( i==2 ){											// Graphical Output
	    if ( strcmp(dat[0],"yes")==0 || strcmp(dat[0],"Y")==0 || strcmp(dat[0],"YES")==0
	    || strcmp(dat[0],"y")==0 || strcmp(dat[0],"true")==0 || strcmp(dat[0],"TRUE")==0 ){
		cfg->gl=1;
	    sscanf(dat[1],"%lg",&dtemp); 
	    if ( dtemp <=0 ) exit_with_usage(7);
	    else cfg->t_gl=dtemp;
       }
       else if ( strcmp(dat[0],"no")==0 || strcmp(dat[0],"N")==0 || strcmp(dat[0],"NO")==0
	    || strcmp(dat[0],"n")==0 || strcmp(dat[0],"false")==0 || strcmp(dat[0],"FALSE")==0 ){
		cfg->gl=0; cfg->t_gl=0; }
          else exit_with_usage(8);
    }
    if ( i==3 ){											// Data output time
	sscanf(dat[0],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(9);
	else cfg->t_write = dtemp;
	sscanf(dat[1],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(9);
	else cfg->t_fits  = dtemp;
	sscanf(dat[2],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(9);
	else cfg->t_sed   = dtemp;
    }    
    if ( i==4 ){											// Evolution endtime
	sscanf(dat[0],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(10);
	else cfg->t_end=dtemp;
    }    
    if ( i==5 ){											// Dynamical timesteps per orbit
	sscanf(dat[0],"%d",&itemp); 
	if ( itemp <=0 ) exit_with_usage(11);
	else cfg->n_step=itemp;
    }    
    if ( i==6 ){ 											// Stellar params
	if ( t<8 ) exit_with_usage(13);
	sscanf(dat[0],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(12);
	else add_star(cfg,dtemp);

	cfg->sptype=(char *)malloc(sizeof(char)*(strlen(dat[1])+1));
	if ( !cfg->sptype )		exit_with_usage(49);
        cfg->mem += sizeof(char)*(strlen(dat[1])+1);

	strcpy(cfg->sptype,dat[1]);

	sscanf(dat[2],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(12);
	cfg->teff=dtemp;
	sscanf(dat[3],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(12);
	cfg->rstar=dtemp*RSun/AU;
	sscanf(dat[4],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(12);
	cfg->Bstar=dtemp*1.0e-4*year*year/MSun;
	sscanf(dat[5],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(12);
	cfg->Wz=365.2422/dtemp;
	sscanf(dat[6],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(12);
	cfg->mloss=dtemp;
	sscanf(dat[7],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(12);
	cfg->vsw=dtemp*1000.0;
    }
    if ( i==7 ){
	if ( t<7 ) exit_with_usage(36);

	sscanf(dat[0],"%lg",&dtemp);
	if ( dtemp <=0 ) exit_with_usage(37);
	cfg->dist = dtemp * pc;
	sscanf(dat[1],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(38);
	cfg->fov  = dtemp;
	sscanf(dat[2],"%lg",&dtemp);
	cfg->PA   = dtemp * PI/180.0;
	sscanf(dat[3],"%lg",&dtemp);
	cfg->inc  = dtemp * PI/180.0;

	sscanf(dat[4],"%lg",&dtemp);
	cfg->longitude  = dtemp * PI/180.0;
	sscanf(dat[5],"%d",&itemp);
	cfg->xsize  = itemp;
	sscanf(dat[6],"%d",&itemp);
	cfg->ysize  = itemp;

    }
    if ( i==8 ){
	if ( t<1 ) exit_with_usage(41);
	cfg->nwavs = t;
	cfg->wavs = (double *)malloc(sizeof(double)*t);
	if ( !cfg->wavs )	exit_with_usage(49);
        cfg->mem += sizeof(double)*t;
	cfg->wj   = (int *)malloc(sizeof(int)*t);
	if ( !cfg->wj )		exit_with_usage(49);
        cfg->mem += sizeof(int)*t;
	for ( int k=0 ; k<cfg->nwavs ; k++ ) sscanf(dat[k],"%lg",&cfg->wavs[k]);
	for ( int k=0 ; k<cfg->nwavs ; k++ ) cfg->wavs[k]*=1e-6;
    }
    if ( i==9 ){
	int length = 5;
        if ( cfg->teff < 9875 ) length = 4;
        if ( strcmp(dat[0],"Library")==0 || strcmp(dat[0],"library")==0 || strcmp(dat[0],"LIBRARY")==0 ){
	    if      ( cfg->teff >= 3375  && cfg->teff < 10000 )						// Library or file	
		cfg->teff = round( cfg->teff / 250.0  ) * 250;
	    else if ( cfg->teff >= 10000 && cfg->teff < 13000 )
		cfg->teff = round( cfg->teff / 500.0  ) * 500;
	    else if ( cfg->teff >= 13000 && cfg->teff < 35000 )
		cfg->teff = round( cfg->teff / 1000.0 ) * 1000;
	    else if ( cfg->teff >= 35000 && cfg->teff < 50000 )
		cfg->teff = round( cfg->teff / 2500.0 ) * 2500;
	    else
		exit_with_usage(40);
	    cfg->specfile = (char *)malloc(length+18+1);
	    if ( !cfg->specfile )	exit_with_usage(49);
    	    cfg->mem += length+18+1;
	    sprintf(cfg->specfile,"Kurucz/Kurucz_%.0f.dat",cfg->teff);
	}
	else{
	    cfg->specfile=(char *)malloc(strlen(dat[0])+1);
	    if ( !cfg->specfile )	exit_with_usage(49);
	    strcpy(cfg->specfile,dat[0]);
    	    cfg->mem += strlen(dat[0])+1;
	}
    }
    if ( i==10 ){											// Physics to include
	if ( t<3 ) exit_with_usage(30);
	    if ( strcmp(dat[0],"yes")==0 || strcmp(dat[0],"Y")==0 || strcmp(dat[0],"YES")==0		// Radiation?
	    || strcmp(dat[0],"y")==0 || strcmp(dat[0],"true")==0 || strcmp(dat[0],"TRUE")==0 )	
		cfg->radQ=1;
       else if ( strcmp(dat[0],"no")==0 || strcmp(dat[0],"N")==0 || strcmp(dat[0],"NO")==0
	    || strcmp(dat[0],"n")==0 || strcmp(dat[0],"false")==0 || strcmp(dat[0],"FALSE")==0 )
		cfg->radQ=0;
       else exit_with_usage(31);

	    if ( strcmp(dat[1],"yes")==0 || strcmp(dat[1],"Y")==0 || strcmp(dat[1],"YES")==0		// Stellar-Wind drag
	    || strcmp(dat[1],"y")==0 || strcmp(dat[1],"true")==0 || strcmp(dat[1],"TRUE")==0 )	
		cfg->SWQ=1;
       else if ( strcmp(dat[1],"no")==0 || strcmp(dat[1],"N")==0 || strcmp(dat[1],"NO")==0
	    || strcmp(dat[1],"n")==0 || strcmp(dat[1],"false")==0 || strcmp(dat[1],"FALSE")==0 )
		cfg->SWQ=0;
       else exit_with_usage(31);

	    if ( strcmp(dat[2],"yes")==0 || strcmp(dat[2],"Y")==0 || strcmp(dat[2],"YES")==0		// Magnetic forces
	    || strcmp(dat[2],"y")==0 || strcmp(dat[2],"true")==0 || strcmp(dat[2],"TRUE")==0 )	
		cfg->FLQ=1;
       else if ( strcmp(dat[2],"no")==0 || strcmp(dat[2],"N")==0 || strcmp(dat[2],"NO")==0
	    || strcmp(dat[2],"n")==0 || strcmp(dat[2],"false")==0 || strcmp(dat[2],"FALSE")==0 )
		cfg->FLQ=0;
       else exit_with_usage(31);
    }
    if ( i==11 ){ // Here come the planets!								// Planets
	if ( strcmp(dat[0],"begin_planets")==0 ) i--;
	else if ( strcmp(dat[0],"end_planets")==0 ){ recalc_COM(cfg); }					// Recalculate system to Center of Mass
	else{ if ( t<7 ) exit_with_usage(18);
	    add_planets(cfg,dat); i--; }	
    }
    if ( i==12 ){											// density and composition
	if ( t<3 ) exit_with_usage(15);
	sscanf(dat[0],"%lg",&dtemp); 
	if ( dtemp <=0 ) exit_with_usage(14);
	else cfg->rho=dtemp*1000;
	cfg->comp=(char *)malloc(strlen(dat[1])+1);
	if ( !cfg->comp )	exit_with_usage(49);
        cfg->mem += strlen(dat[1])+1;
	strcpy(cfg->comp,dat[1]);
	sscanf(dat[2],"%lg",&dtemp);
	if ( dtemp <=0 && cfg->FLQ==1 ) exit_with_usage(32);
	else cfg->kqe=dtemp/year;
    }
    if ( i==13 ){											// Disks
	if ( strcmp(dat[0],"begin_disks")==0 ) i--;
	else if ( strcmp(dat[0],"end_disks")==0 ) ; 
	else{ add_disks(cfg,dat); i--; }	
    }
    if ( i==14 ){
	if ( strcmp(dat[0],"begin_blobs")==0 ) i--; 							// Blobs
	else if ( strcmp(dat[0],"end_blobs")==0 ) ; 
	else{ add_blobs(cfg,dat,t); i--; }
    }
    if ( i==15 ){											// stubs
	cfg->name_stub=(char *)malloc(strlen(dat[0])+1);
	if ( !cfg->name_stub )	exit_with_usage(49);
	strcpy(cfg->name_stub,dat[0]);
        cfg->mem += strlen(dat[0])+1;
    }    
 i++;
 }
 if ( i<15 )            exit_with_usage(17);	// new number

 fclose(fr);
 
 return 0;
}

// Read in model file
int read_in_model(char * model,const char *ver){
 FILE        *fr;
 int         i,t,itemp;
 double	     dtemp;
 char        buff[2048],*dat[64];
 configdata  *cfg=&cconfig;

 fr=fopen(model,"rb");
 if ( fr==NULL )        exit_with_usage(3);

 cfg->ram  = getMemorySize();
 cfg->mem += sizeof(configdata);
 cfg->prev = 0;

 i=0;
 while ( ! feof(fr) ){
    if ( fgets(buff,2048,fr)==NULL )     break;
    t=tokenize_spaces(buff,dat,32);
    if ( t==0 )				continue;
    if ( i==1 ){
	if ( strcmp(ver,dat[1])!=0 )	exit_with_usage(55);
    }
    if ( i==4 ){
	sscanf(dat[1],"%d",&itemp);
	if ( itemp <0 ) exit_with_usage(4);
	else cfg->verb = itemp;
    }
    if ( i==5 ){
	sscanf(dat[1],"%d",&itemp);
	if ( itemp <0 ) exit_with_usage(6);
	else cfg->gpuid = itemp;
    }
    if ( i==6 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->gl = itemp;
    }
    if ( i==7 ){
        sscanf(dat[1],"%lg",&dtemp);
	cfg->t_gl = dtemp;
    }
    if ( i==8 ){
        sscanf(dat[1],"%lg",&dtemp);
	cfg->t_write = dtemp;
    }
    if ( i==9 ){
        sscanf(dat[1],"%lg",&dtemp);
	cfg->t_sed = dtemp;
    }
    if ( i==10 ){
        sscanf(dat[1],"%lg",&dtemp);
	cfg->t_fits = dtemp;
    }
    if ( i==11 ){
        sscanf(dat[1],"%lg",&dtemp);
	cfg->t_end = dtemp;
    }
    if ( i==12 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->n_step = itemp;
    }
    if ( i==13 ){
	cfg->sptype=(char *)malloc(sizeof(char)*(strlen(dat[1])+1));
	if ( !cfg->sptype )		exit_with_usage(49);
        cfg->mem += sizeof(char)*(strlen(dat[1])+1);
	strcpy(cfg->sptype,dat[1]);
    }
    if ( i==14 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->teff = dtemp;
    }
    if ( i==15 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->rstar = dtemp*RSun/AU;
    }
    if ( i==16 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->Bstar = dtemp*year*year*1.0e-4/MSun;
    }
    if ( i==17 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->Wz = 365.2422/dtemp;
    }
    if ( i==18 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->mloss = dtemp;
    }
    if ( i==19 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->vsw = dtemp;
    }
    if ( i==20 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->radQ = itemp;
    }
    if ( i==21 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->SWQ = itemp;
    }
    if ( i==22 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->FLQ = itemp;
    }
    if ( i==23 ){
	sscanf(dat[1],"%lg",&dtemp); 
	cfg->rho = dtemp;
    }
    if ( i==24 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->kqe = dtemp/year;
    }
    if ( i==25 ){
	cfg->name_stub=(char *)malloc(strlen(dat[1])+1);
	if ( !cfg->name_stub )	exit_with_usage(49);
	strcpy(cfg->name_stub,dat[1]);
        cfg->mem += strlen(dat[1])+1;
    }
    if ( i==26 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->dist = dtemp * pc;
    }
    if ( i==27 ){
	sscanf(dat[1],"%lg",&dtemp); 
	cfg->fov  = dtemp;
    }
    if ( i==28 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->PA   = dtemp * PI/180.0;
    }
    if ( i==29 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->inc  = dtemp * PI/180.0;
    }
    if ( i==30 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->longitude  = dtemp * PI/180.0;
    }
    if ( i==31 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->xsize  = itemp;
    }
    if ( i==32 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->ysize  = itemp;
    }
    if ( i==33 ){
        cfg->specfile=(char *)malloc(strlen(dat[1])+1);
        if ( !cfg->specfile )	exit_with_usage(49);
        strcpy(cfg->specfile,dat[1]);
	cfg->mem += strlen(dat[1])+1;
    }
    if ( i==34 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->amin = dtemp;
    }
    if ( i==35 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->amax = dtemp;
    }
    if ( i==36 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->mmin = dtemp;
    }
    if ( i==37 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->mmax = dtemp;
    }
    if ( i==38 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->n_grid = itemp;
    }
    if ( i==39 ){
	sscanf(dat[1],"%lg",&dtemp);
	cfg->f = dtemp;
    }
    if ( i==40 ){
	cfg->comp=(char *)malloc(strlen(dat[1])+1);
	if ( !cfg->comp )	exit_with_usage(49);
        cfg->mem += strlen(dat[1])+1;
	strcpy(cfg->comp,dat[1]);
    }
    if ( i==41 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->n_wav = itemp;
    }
    if ( i==42 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->n_theta = itemp;
    }
    if ( i==43 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->n_dust = itemp;
    }
    if ( i==44 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->n_planets = itemp;
	cfg->planets   = (planet_cart *)malloc(sizeof(planet_cart)*cfg->n_planets);
	if ( cfg->planets == NULL ) exit_with_usage(44);
	cfg->mem += sizeof(planet_cart)*cfg->n_planets;
    }
    if ( i==45 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->n_disks = itemp;
	cfg->disks   = (diskconf *)malloc(sizeof(diskconf)*cfg->n_disks);
	if ( cfg->disks == NULL ) exit_with_usage(45);
	cfg->mem += sizeof(diskconf)*cfg->n_disks;
    }
    if ( i==46 ){
	sscanf(dat[1],"%d",&itemp);
	cfg->n_blobs = itemp;
	cfg->blobs   = (blobconf *)malloc(sizeof(blobconf)*cfg->n_blobs);
	if ( cfg->blobs == NULL ) exit_with_usage(46);
	cfg->mem += sizeof(blobconf)*cfg->n_blobs;
    }
    if ( i==47 ){
	if ( t<2 ) exit_with_usage(41);
	cfg->nwavs = t-1;
	cfg->wavs = (double *)malloc(sizeof(double)*(t-1));
	if ( !cfg->wavs )	exit_with_usage(49);
        cfg->mem += sizeof(double)*(t-1);
	cfg->wj   = (int *)malloc(sizeof(int)*(t-1));
	if ( !cfg->wj )		exit_with_usage(49);
        cfg->mem += sizeof(int)*(t-1);
	for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->wavs[k-1]);
	for ( int k=1 ; k<t ; k++ ) cfg->wavs[k-1]*=1e-6;
    }
    if ( i==49 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].Rin);
    if ( i==50 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].Rout);
    if ( i==51 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].dhdr);
    if ( i==52 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].sigmae);
    if ( i==53 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].amax);
    if ( i==54 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].mmax);
    if ( i==55 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].amin);
    if ( i==56 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].mmin);
    if ( i==57 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].dmax);
    if ( i==58 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].dmmax);
    if ( i==59 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].eta);
    if ( i==60 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].gammas);
    if ( i==61 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].C);
    if ( i==62 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->disks[k-1].Mtot);
    if ( i==63 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%ld",&cfg->disks[k-1].Ngrav);
    if ( i==64 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%ld",&cfg->disks[k-1].Nsmall);
    if ( i==65 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].a);
    if ( i==66 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].ecc);
    if ( i==67 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].inc);
    if ( i==68 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].w);
    if ( i==69 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].o);
    if ( i==70 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].nu);
    if ( i==71 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].sigma);
    if ( i==72 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].amax);
    if ( i==73 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].mmax);
    if ( i==74 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].amin);
    if ( i==75 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].mmin);
    if ( i==76 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].dmax);
    if ( i==77 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].dmmax);
    if ( i==78 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].eta);
    if ( i==79 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].C);
    if ( i==80 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].Mtot);
    if ( i==81 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%ld",&cfg->blobs[k-1].Ngrav);
    if ( i==82 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%ld",&cfg->blobs[k-1].Nsmall);
    if ( i==83 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%d",&cfg->blobs[k-1].collcloud);
    if ( i==84 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->blobs[k-1].vkepp);
    if ( i==87 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->planets[k-1].x);
    if ( i==88 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->planets[k-1].y);
    if ( i==89 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->planets[k-1].z);
    if ( i==90 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->planets[k-1].vx);
    if ( i==91 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->planets[k-1].vy);
    if ( i==92 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->planets[k-1].vz);
    if ( i==93 ) for ( int k=1 ; k<t ; k++ ) sscanf(dat[k],"%lg",&cfg->planets[k-1].m);
    if ( i==94 ){
	sscanf(dat[1],"%d",&itemp);
	if ( itemp != nR ) exit_with_usage(56);
    }
    if ( i>94) break;
 i++;
 }
 fclose(fr);
 
 return 0;
}