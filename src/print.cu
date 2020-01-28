#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddyn.cuh"

void print_header(const char * Version,const char *VDate){
 configdata	*cfg = &cconfig;

 if ( cfg->verb>0 ){
 int i,length;

 length=strlen(Version)+strlen(VDate);

 printf("\n");
 printf("#########################################################################\n");
 printf("#                          Running DDyn!                                #\n");
 printf("#                        Ver %s, %s",Version,VDate);
 for ( i=0 ; i<41-length ; i++ ) printf(" ");printf("#\n");
 printf("#                                                                       #\n");
 printf("# Please contact Andras Gaspar if there seems to be a problem at:       #\n");
 printf("# agaspar@as.arizona.edu                                                #\n");
 printf("#                                                                       #\n");
 printf("#########################################################################\n");
 printf("\n");
 }
}


void print_cuda(){
 cudadata	*cuda = &ccuda;
 configdata	*cfg  = &cconfig;

 if ( cfg->verb>0 ){
    printf("Running on the following GPU:\n\n");
    printf("\tGPU id: ................................. %d\n",cuda->prop[cfg->gpuid].gpu);
    printf("\tGPU type: ............................... %s\n",cuda->prop[cfg->gpuid].name);
    printf("\tCuda version: ........................... %d.%d\n",cuda->prop[cfg->gpuid].major,cuda->prop[cfg->gpuid].minor);
    printf("\tNumber of CUDA streaming multiprocessors: %d\n",cuda->prop[cfg->gpuid].sm);
    printf("\tNumber of CUDA cores: ................... %d\n",cuda->prop[cfg->gpuid].cores);
    printf("\tTotal Global Memory: .................... %.3f Gb\n",cuda->prop[cfg->gpuid].totalGlobalMem/1.0e9);
    printf("\tTotal Constant Memory: .................. %ld bytes\n",cuda->prop[cfg->gpuid].totalConstMem);
    printf("\tMemory Pitch: ........................... %ld bytes\n",cuda->prop[cfg->gpuid].memPitch);
    printf("\tShared Memory Per Block: ................ %ld bytes\n",cuda->prop[cfg->gpuid].sharedMemPerBlock);
    printf("\tMax Threads Per Block: .................. %d bytes\n",cuda->prop[cfg->gpuid].maxThreadsPerBlock);
    printf("\n");
 }
}

void print_setup(){
 configdata	*cfg = &cconfig;

 if ( cfg->verb<2 ) return;
 printf("Setting up initial particle grid onto host and GPUs\n");

}

void print_optics(){
 configdata	*cfg = &cconfig;

 if ( cfg->verb<2 ) return;
 printf("Reading in optical data parameters\n");
 printf("This may take some time!\n");

}

void printmin(double m,int j,int s){
 configdata	*cfg = &cconfig;

    if ( cfg->verb<1 ) return;
    double r=pow(m*3.0/4.0/PI/cfg->rho,1.0/3.0);
    if ( s==0 ) printf("Smallest massive object followed in disk %d has a radius of %.4f km and a mass of %.4e kg\n",j,r/1000,m);
    if ( s==1 ) printf("Smallest massive object followed in blob %d has a radius of %.4f km and a mass of %.4e kg\n",j,r/1000,m);
}

void printstatus(double time){
    configdata	*cfg = &cconfig;
    if ( cfg->verb<3 ) return;

    printf("Model at time: %.3f yr\r",time);
    fflush(stdout);
}

void printout(double time){
    configdata	*cfg = &cconfig;
    if ( cfg->verb<2 ) return;

    printf("Printing at time: %.6f\n",time);
    fflush(stdout);
}

void printtimestep(double time){
    configdata	*cfg = &cconfig;
    if ( cfg->verb<1 ) return;

    printf("Step size: %.6f year\n",time);
    fflush(stdout);
}

void print_reset_min(int i,int s){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<1 ) return;

    if ( s==0 ) printf("Resetting smallest particles in disk %d to %.4f micron to match optical library\n",i,cfg->dmin*1e6);
    if ( s==1 ) printf("Resetting smallest particles in blob %d to %.4f micron to match optical library\n",i,cfg->dmin*1e6);
}

void print_reset_max(int i,int s){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<1 ) return;

    if ( s==0 ){
	if ( fabs(cfg->disks[i].amax-cfg->dmax)<1e-12 ) return;		// Numerical precision issue only;
	printf("Setting largest emitting dust particles in disk %d to %.2f micron according to optical library\n",i,cfg->dmax*1e6);
    }
    if ( s==1 ){
	if ( fabs(cfg->blobs[i].dmax-cfg->dmax)<1e-12 ) return;		// Numerical precision issue only;
	printf("Setting largest emitting dust particles in blob %d to %.2f micron according to optical library\n",i,cfg->dmax*1e6);
    }
}

void printmemall(){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<2 ) return;
    printf("System Memory usage Summary\n");
    printf("\tTotal System RAM: .......................... %.2f Gb\n",cfg->ram/1e9);
    printf("\tTotal system memory used: .................. %.2f Gb\n",cfg->mem/1e9);
    printf("\tTotal GPU memory used: ..................... %.2f Gb\n",cfg->memgpu/1e9);
    printf("\tTotal GPU memory used for optical constants: %.2f Gb\n",cfg->memgpudust/1e9);
    printf("\tTotal GPU memory used for particles: ....... %.2f Gb\n",(cfg->memgpu-cfg->memgpudust)/1e9);
    printf("\tTotal constant memory used: ................ %lu bytes\n",nsize*3*sizeof(double));
    printf("\n");
}

void printsystemparams(){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<2 ) return;
    printf("\n#####################################\n");
    printf("The following system is evaluated\n");
    printf("#####################################\n");

    printf("System distance:     %.2f pc\n\n",cfg->dist/pc);

    printf("Stellar parameters\n");
    printf("\tStellar Luminosity: ...................... %.2f LSun\n",cfg->Lum);
    printf("\tStellar Temperature: ..................... %.2f K\n",cfg->teff);
    printf("\tStellar Radius: .......................... %.2f RSun\n",cfg->rstar*AU/RSun);
    printf("\tStellar Mass: ............................ %.2f MSun\n\n",cfg->planets[0].m/Grav);

    if ( cfg->n_planets > 1)
	printf("Planetary parameters\n");
    for ( int i=1 ; i<cfg->n_planets ; i++){
	printf("\tPlanet %d at R: %.2f AU m: %.4f M_Earth\n",i,sqrt(cfg->planets[i].x*cfg->planets[i].x+cfg->planets[i].y*cfg->planets[i].y+cfg->planets[i].z*cfg->planets[i].z),cfg->planets[i].m*MSun/MEarth/Grav);
    }

    if ( cfg->n_disks > 0 )
	printf("Disk parameters\n");
    for ( int i=0 ; i<cfg->n_disks ; i++){
	printf("\tDisk %d Rin: %.2f AU Rout: %.2f AU\n",i+1,cfg->disks[i].Rin,cfg->disks[i].Rout);
    }
    if ( cfg->n_blobs > 0 )
	printf("Blob parameters\n");
    for ( int i=0 ; i<cfg->n_blobs ; i++){
	printf("\tBlob %d a: %.2f AU size: %.2f AU\n",i+1,cfg->blobs[i].a,cfg->blobs[i].sigma);
    }
    printf("Dust parameters\n");
    printf("\tSmallest dust size modeled: .............. %.3f micron\n",cfg->amin*1e6);
    printf("\tLargest dust size modeled: ............... %.3f micron\n",cfg->amax*1e6);
    printf("\tNumber of dust grids modeled: ............ %d\n\n",cfg->n_grid);

}

void print_status(int t,int l){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<3 ) return;
    int i;

    if ( (int)(100.0*(double)t/((double)l)) > cfg->prev ){
	printf("Status: [");
	for ( i=0 ; i < 50*t/l    ; i++ ) printf("#");
	for ( i=50*t/l ; i < 50 ; i++   ) printf(" ");
	printf("] %2.0f %%\r",100.0*(double)t/((double) l));

	if ( t==l ) printf("\n"); 
	fflush(stdout);
	cfg->prev = (int)(100.0*(double)t/((double)l));
    }

    return;
}

void printdata(double time){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<2 ) return;
    printf("Printing data file at %.2f year\n",time);
}

void printsed(double time){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<2 ) return;
    printf("Printing SED file at %.2f year\n",time);
}

void printcalcflux(int wj){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<3 ) return;
    printf("Calculting image at %.2f micron\n",cfg->w[wj]*1e6);
}

void printsedcalc(){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<3 ) return;
    printf("Calculting the complete SED of the system!\n");
}

void printnewline(){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<3 ) return;
    printf("\n");
}

void printdone(){
    configdata	*cfg = &cconfig;

    if ( cfg->verb<1 ) return;
    printf("Model is finished!\n");
}

