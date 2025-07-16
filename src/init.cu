#include <stdio.h>
#include "ddyn.cuh"

Grid::Grid(){						// Initialization
    host.x     = NULL;
    host.y     = NULL;
    host.z     = NULL;
    host.vx    = NULL;
    host.vy    = NULL;
    host.vz    = NULL;
    host.m     = NULL;
    host.N     = NULL;
    host.id    = NULL;

    host.k1x   = NULL;
    host.k2x   = NULL;
    host.k3x   = NULL;
    host.k4x   = NULL;
    host.k1y   = NULL;
    host.k2y   = NULL;
    host.k3y   = NULL;
    host.k4y   = NULL;
    host.k1z   = NULL;
    host.k2z   = NULL;
    host.k3z   = NULL;
    host.k4z   = NULL;

    host.k2rx  = NULL;
    host.k3rx  = NULL;
    host.k4rx  = NULL;
    host.k2ry  = NULL;
    host.k3ry  = NULL;
    host.k4ry  = NULL;
    host.k2rz  = NULL;
    host.k3rz  = NULL;
    host.k4rz  = NULL;

    host.k1vx  = NULL;
    host.k2vx  = NULL;
    host.k3vx  = NULL;
    host.k4vx  = NULL;
    host.k1vy  = NULL;
    host.k2vy  = NULL;
    host.k3vy  = NULL;
    host.k4vy  = NULL;
    host.k1vz  = NULL;
    host.k2vz  = NULL;
    host.k3vz  = NULL;
    host.k4vz  = NULL;

    gpu.x      = NULL;					// GPU memory pointer
    gpu.y      = NULL;
    gpu.z      = NULL;
    gpu.vx     = NULL;
    gpu.vy     = NULL;
    gpu.vz     = NULL;
    gpu.m      = NULL;
    gpu.N      = NULL;
    gpu.id     = NULL;

    gpu.k1x    = NULL;
    gpu.k2x    = NULL;
    gpu.k3x    = NULL;
    gpu.k4x    = NULL;
    gpu.k1y    = NULL;
    gpu.k2y    = NULL;
    gpu.k3y    = NULL;
    gpu.k4y    = NULL;
    gpu.k1z    = NULL;
    gpu.k2z    = NULL;
    gpu.k3z    = NULL;
    gpu.k4z    = NULL;

    gpu.k2rx   = NULL;
    gpu.k3rx   = NULL;
    gpu.k4rx   = NULL;
    gpu.k2ry   = NULL;
    gpu.k3ry   = NULL;
    gpu.k4ry   = NULL;
    gpu.k2rz   = NULL;
    gpu.k3rz   = NULL;
    gpu.k4rz   = NULL;

    gpu.k1vx   = NULL;
    gpu.k2vx   = NULL;
    gpu.k3vx   = NULL;
    gpu.k4vx   = NULL;
    gpu.k1vy   = NULL;
    gpu.k2vy   = NULL;
    gpu.k3vy   = NULL;
    gpu.k4vy   = NULL;
    gpu.k1vz   = NULL;
    gpu.k2vz   = NULL;
    gpu.k3vz   = NULL;
    gpu.k4vz   = NULL;

    Th         = NULL;
    Td         = NULL;
    w          = NULL;
    in	       = NULL;
    Qabs       = NULL;
    Qsca       = NULL;
    Qpfc       = NULL;
    Fh	       = NULL;
    Fd         = NULL;
    Fdust      = NULL;

    synced     = false;					// Not synced at this point
}

Grid::~Grid(){

    free(host.x);
    free(host.y);
    free(host.z);
    free(host.vx);
    free(host.vy);
    free(host.vz);
    free(host.m);
    free(host.N);
    free(host.id);

    free(host.k1x);
    free(host.k2x);
    free(host.k3x);
    free(host.k4x);
    free(host.k1y);
    free(host.k2y);
    free(host.k3y);
    free(host.k4y);
    free(host.k1z);
    free(host.k2z);
    free(host.k3z);
    free(host.k4z);

    free(host.k2rx);
    free(host.k3rx);
    free(host.k4rx);
    free(host.k2ry);
    free(host.k3ry);
    free(host.k4ry);
    free(host.k2rz);
    free(host.k3rz);
    free(host.k4rz);

    free(host.k1vx);
    free(host.k2vx);
    free(host.k3vx);
    free(host.k4vx);
    free(host.k1vy);
    free(host.k2vy);
    free(host.k3vy);
    free(host.k4vy);
    free(host.k1vz);
    free(host.k2vz);
    free(host.k3vz);
    free(host.k4vz);
    free(Th);
    free(Fh);

    gpuErrchk(cudaFree(gpu.x));
    gpuErrchk(cudaFree(gpu.y));
    gpuErrchk(cudaFree(gpu.z));
    gpuErrchk(cudaFree(gpu.vx));
    gpuErrchk(cudaFree(gpu.vy));
    gpuErrchk(cudaFree(gpu.vz));
    gpuErrchk(cudaFree(gpu.m));
    gpuErrchk(cudaFree(gpu.N));
    gpuErrchk(cudaFree(gpu.id));

    gpuErrchk(cudaFree(gpu.k1x));
    gpuErrchk(cudaFree(gpu.k2x));
    gpuErrchk(cudaFree(gpu.k3x));
    gpuErrchk(cudaFree(gpu.k4x));
    gpuErrchk(cudaFree(gpu.k1y));
    gpuErrchk(cudaFree(gpu.k2y));
    gpuErrchk(cudaFree(gpu.k3y));
    gpuErrchk(cudaFree(gpu.k4y));
    gpuErrchk(cudaFree(gpu.k1z));
    gpuErrchk(cudaFree(gpu.k2z));
    gpuErrchk(cudaFree(gpu.k3z));
    gpuErrchk(cudaFree(gpu.k4z));

    gpuErrchk(cudaFree(gpu.k2rx));
    gpuErrchk(cudaFree(gpu.k3rx));
    gpuErrchk(cudaFree(gpu.k4rx));
    gpuErrchk(cudaFree(gpu.k2ry));
    gpuErrchk(cudaFree(gpu.k3ry));
    gpuErrchk(cudaFree(gpu.k4ry));
    gpuErrchk(cudaFree(gpu.k2rz));
    gpuErrchk(cudaFree(gpu.k3rz));
    gpuErrchk(cudaFree(gpu.k4rz));

    gpuErrchk(cudaFree(gpu.k1vx));
    gpuErrchk(cudaFree(gpu.k2vx));
    gpuErrchk(cudaFree(gpu.k3vx));
    gpuErrchk(cudaFree(gpu.k4vx));
    gpuErrchk(cudaFree(gpu.k1vy));
    gpuErrchk(cudaFree(gpu.k2vy));
    gpuErrchk(cudaFree(gpu.k3vy));
    gpuErrchk(cudaFree(gpu.k4vy));
    gpuErrchk(cudaFree(gpu.k1vz));
    gpuErrchk(cudaFree(gpu.k2vz));
    gpuErrchk(cudaFree(gpu.k3vz));
    gpuErrchk(cudaFree(gpu.k4vz));

    gpuErrchk(cudaFree(w));
    gpuErrchk(cudaFree(in));
    gpuErrchk(cudaFree(Td));
    gpuErrchk(cudaFree(Qabs));
    gpuErrchk(cudaFree(Qsca));
    gpuErrchk(cudaFree(Qpfc));
    gpuErrchk(cudaFree(Fd));
    gpuErrchk(cudaFree(Fdust));
}

void Grid::allochost(){
    configdata *cfg  = &cconfig;

    if (!host.x   )  host.x      = (double *)malloc(sizeof(double)*n_all);
	if ( !host.x )		 exit_with_usage(53);
    if (!host.y   )  host.y      = (double *)malloc(sizeof(double)*n_all);
	if ( !host.y )		 exit_with_usage(53);
    if (!host.z   )  host.z      = (double *)malloc(sizeof(double)*n_all);
	if ( !host.z )		 exit_with_usage(53);
    if (!host.vx  )  host.vx     = (double *)malloc(sizeof(double)*n_all);
	if ( !host.vx )		 exit_with_usage(53);
    if (!host.vy  )  host.vy     = (double *)malloc(sizeof(double)*n_all);
	if ( !host.vy )		 exit_with_usage(53);
    if (!host.vz  )  host.vz     = (double *)malloc(sizeof(double)*n_all);
	if ( !host.vz )	         exit_with_usage(53);
    if (!host.m   )  host.m      = (double *)malloc(sizeof(double)*n_all);
	if ( !host.m )		 exit_with_usage(53);
    if (!host.N   )  host.N      = (double *)malloc(sizeof(double)*n_all);
	if ( !host.N )		 exit_with_usage(53);
    if (!host.id  )  host.id     = (int    *)malloc(sizeof(int)   *n_all);
	if ( !host.id )		 exit_with_usage(53);
    if (!Fh       )  Fh          = (double *)malloc(sizeof(double)*cfg->xsize*cfg->ysize);
	if ( !Fh )		 exit_with_usage(53);

    if (!host.k1x )  host.k1x    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k1x )	 exit_with_usage(53);
    if (!host.k2x )  host.k2x    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2x )	 exit_with_usage(53);
    if (!host.k3x )  host.k3x    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3x )	 exit_with_usage(53);
    if (!host.k4x )  host.k4x    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4x )	 exit_with_usage(53);
    if (!host.k1y )  host.k1y    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k1y )	 exit_with_usage(53);
    if (!host.k2y )  host.k2y    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2y )	 exit_with_usage(53);
    if (!host.k3y )  host.k3y    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3y )	 exit_with_usage(53);
    if (!host.k4y )  host.k4y    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4y )	 exit_with_usage(53);
    if (!host.k1z )  host.k1z    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k1z )	 exit_with_usage(53);
    if (!host.k2z )  host.k2z    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2z )	 exit_with_usage(53);
    if (!host.k3z )  host.k3z    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3z )	 exit_with_usage(53);
    if (!host.k4z )  host.k4z    = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4z )	 exit_with_usage(53);

    if (!host.k2rx)  host.k2rx   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2rx )	 exit_with_usage(53);
    if (!host.k3rx)  host.k3rx   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3rx )	 exit_with_usage(53);
    if (!host.k4rx)  host.k4rx   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4rx )	 exit_with_usage(53);
    if (!host.k2ry)  host.k2ry   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2ry )	 exit_with_usage(53);
    if (!host.k3ry)  host.k3ry   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3ry )	 exit_with_usage(53);
    if (!host.k4ry)  host.k4ry   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4ry )	 exit_with_usage(53);
    if (!host.k2rz)  host.k2rz   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2rz )	 exit_with_usage(53);
    if (!host.k3rz)  host.k3rz   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3rz )	 exit_with_usage(53);
    if (!host.k4rz)  host.k4rz   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4rz )	 exit_with_usage(53);

    if (!host.k1vx)  host.k1vx   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k1vx )	 exit_with_usage(53);
    if (!host.k2vx)  host.k2vx   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2vx )	 exit_with_usage(53);
    if (!host.k3vx)  host.k3vx   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3vx )	 exit_with_usage(53);
    if (!host.k4vx)  host.k4vx   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4vx )	 exit_with_usage(53);
    if (!host.k1vy)  host.k1vy   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k1vy )	 exit_with_usage(53);
    if (!host.k2vy)  host.k2vy   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2vy )	 exit_with_usage(53);
    if (!host.k3vy)  host.k3vy   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3vy )	 exit_with_usage(53);
    if (!host.k4vy)  host.k4vy   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4vy )	 exit_with_usage(53);
    if (!host.k1vz)  host.k1vz   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k1vz )	 exit_with_usage(53);
    if (!host.k2vz)  host.k2vz   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k2vz )	 exit_with_usage(53);
    if (!host.k3vz)  host.k3vz   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k3vz )	 exit_with_usage(53);
    if (!host.k4vz)  host.k4vz   = (double *)malloc(sizeof(double)*n_all);
	if ( !host.k4vz )	 exit_with_usage(53);

    cfg->mem +=(sizeof(double)*41+sizeof(int))*n_all+sizeof(double)*cfg->xsize*cfg->ysize;

    null_host();
}

void Grid::reallochost(int plus){
    configdata *cfg  = &cconfig;
    double *dtmp;
    int	   *itmp;

    dtmp        = (double *)realloc(host.x   ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.x = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.y   ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.y = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.z   ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.z = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.vx   ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.vx = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.vy   ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.vy = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.vz   ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.vz = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.m   ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.m = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.N   ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.N = dtmp;
    else	exit_with_usage(53);

    itmp        = (int *)realloc(host.id     ,sizeof(int)*(n_all+plus));
    if (itmp)	host.id = itmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k1x ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k1x = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k2x ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2x = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3x ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3x = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4x ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4x = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k1y ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k1y = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k2y ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2y = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3y ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3y = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4y ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4y = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k1z ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k1z = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k2z ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2z = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3z ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3z = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4z ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4z = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k2rx ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2rx = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3rx ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3rx = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4rx ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4rx = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k2ry ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2ry = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3ry ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3ry = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4ry ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4ry = dtmp;
    else	exit_with_usage(53);
                                    
    dtmp        = (double *)realloc(host.k2rz ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2rz = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3rz ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3rz = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4rz ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4rz = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k1vx ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k1vx = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k2vx ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2vx = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3vx ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3vx = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4vx ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4vx = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k1vy ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k1vy = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k2vy ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2vy = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3vy ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3vy = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4vy ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4vy = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k1vz ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k1vz = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k2vz ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k2vz = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k3vz ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k3vz = dtmp;
    else	exit_with_usage(53);

    dtmp        = (double *)realloc(host.k4vz ,sizeof(double)*(n_all+plus));
    if (dtmp)	host.k4vz = dtmp;
    else	exit_with_usage(53);

    cfg->mem += (sizeof(double)*41+sizeof(int))*plus;
}

void Grid::null_host(){

    for ( int j=0 ; j<n_all ; j++ ){
        host.x[j]    = 0;
	host.y[j]    = 0;
	host.z[j]    = 0;
	host.vx[j]   = 0;
	host.vy[j]   = 0;
	host.vz[j]   = 0;
	host.m[j]    = 0;
        host.N[j]    = 0;
        host.id[j]   = 0;
        
        host.k1x[j]  = 0;
        host.k2x[j]  = 0;
        host.k3x[j]  = 0;
        host.k4x[j]  = 0;
        host.k1y[j]  = 0;
        host.k2y[j]  = 0;
        host.k3y[j]  = 0;
        host.k4y[j]  = 0;
        host.k1z[j]  = 0;
        host.k2z[j]  = 0;
        host.k3z[j]  = 0;
        host.k4z[j]  = 0;
        
        host.k2rx[j] = 0;
        host.k3rx[j] = 0;
        host.k4rx[j] = 0;
        host.k2ry[j] = 0;
        host.k3ry[j] = 0;
        host.k4ry[j] = 0;
        host.k2rz[j] = 0;
        host.k3rz[j] = 0;
        host.k4rz[j] = 0;
        
        host.k1vx[j] = 0;
        host.k2vx[j] = 0;
        host.k3vx[j] = 0;
        host.k4vx[j] = 0;
        host.k1vy[j] = 0;
        host.k2vy[j] = 0;
        host.k3vy[j] = 0;
        host.k4vy[j] = 0;
        host.k1vz[j] = 0;
        host.k2vz[j] = 0;
        host.k3vz[j] = 0;
        host.k4vz[j] = 0;
    }
    synced   = false;
}

void Grid::allocgpu(){
    configdata *cfg  = &cconfig;

    gpuErrchk(cudaMalloc(&gpu.x   ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.y   ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.z   ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.vx  ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.vy  ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.vz  ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.m   ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.N   ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.id  ,sizeof(int)   *n_all));

    gpuErrchk(cudaMalloc(&gpu.k1x ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k2x ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3x ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4x ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k1y ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k2y ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3y ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4y ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k1z ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k2z ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3z ,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4z ,sizeof(double)*n_all)); //12

    gpuErrchk(cudaMalloc(&gpu.k2rx,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3rx,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4rx,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k2ry,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3ry,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4ry,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k2rz,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3rz,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4rz,sizeof(double)*n_all)); //9

    gpuErrchk(cudaMalloc(&gpu.k1vx,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k2vx,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3vx,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4vx,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k1vy,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k2vy,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3vy,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4vy,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k1vz,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k2vz,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k3vz,sizeof(double)*n_all));
    gpuErrchk(cudaMalloc(&gpu.k4vz,sizeof(double)*n_all)); //12 - 41 double 1 int,

    gpuErrchk(cudaMalloc(&w    ,sizeof(double)*cfg->n_wav));
    gpuErrchk(cudaMalloc(&in   ,sizeof(double)*cfg->n_wav));
    gpuErrchk(cudaMalloc(&Td   ,sizeof(double)*cfg->n_dust*nR));
    gpuErrchk(cudaMalloc(&Fd   ,sizeof(double)*cfg->xsize*cfg->ysize));
    gpuErrchk(cudaMalloc(&Qabs ,sizeof(double)*cfg->n_dust*cfg->n_wav));
    gpuErrchk(cudaMalloc(&Qsca ,sizeof(double)*cfg->n_dust*cfg->n_wav));
    gpuErrchk(cudaMalloc(&Qpfc ,sizeof(double)*cfg->n_dust*cfg->n_wav*cfg->n_theta));
    gpuErrchk(cudaMalloc(&Fdust,sizeof(double)*n_dust));

    synced   = false;
}


void Grid::add_to_host(double x,double y,double z,double vx,double vy,double vz,\
                         double N,double m,int id,int tot){

    print_status(n_all,tot-1);
    host.x[n_all]    = x;
    host.y[n_all]    = y;
    host.z[n_all]    = z;
    host.vx[n_all]   = vx;
    host.vy[n_all]   = vy;
    host.vz[n_all]   = vz;
    host.id[n_all]   = id;
    host.m[n_all]    = m;
    host.N[n_all]    = N;
    n_all++;

}

void Grid::sync_to_gpu(){
    if (!synced ){

    gpuErrchk(cudaMemcpy(gpu.x   ,host.x   ,sizeof(double)*n_all,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu.y   ,host.y   ,sizeof(double)*n_all,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu.z   ,host.z   ,sizeof(double)*n_all,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu.vx  ,host.vx  ,sizeof(double)*n_all,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu.vy  ,host.vy  ,sizeof(double)*n_all,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu.vz  ,host.vz  ,sizeof(double)*n_all,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu.m   ,host.m   ,sizeof(double)*n_all,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu.N   ,host.N   ,sizeof(double)*n_all,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu.id  ,host.id  ,sizeof(int)   *n_all,cudaMemcpyHostToDevice));

    synced = true;
    }
}

void Grid::sync_dust_to_gpu(){
    configdata *cfg  = &cconfig;

    gpuErrchk(cudaMemcpy(Qabs    ,cfg->Qabs       ,sizeof(double)*cfg->n_dust*cfg->n_wav,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(Qsca    ,cfg->Qsca       ,sizeof(double)*cfg->n_dust*cfg->n_wav,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(Qpfc    ,cfg->Qpfunc     ,sizeof(double)*cfg->n_dust*cfg->n_wav*cfg->n_theta,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(Td      ,Th              ,sizeof(double)*cfg->n_dust*nR,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(w       ,cfg->w          ,sizeof(double)*cfg->n_wav,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(in      ,cfg->in         ,sizeof(double)*cfg->n_wav,cudaMemcpyHostToDevice));

    double *psize;
    float *bvalue;
    float *bswvalue;
    float *qvalue;
    cudaGetSymbolAddress((void **)&psize,size); // get a pointer to size
    cudaMemcpy(psize,cfg->s,sizeof(double)*cfg->n_dust,cudaMemcpyHostToDevice);
    cudaGetSymbolAddress((void **)&bvalue,beta); // get a pointer to beta 
    cudaMemcpy(bvalue,cfg->beta,sizeof(float)*cfg->n_dust,cudaMemcpyHostToDevice);
    cudaGetSymbolAddress((void **)&bswvalue,betasw); // get a pointer to betasw 
    cudaMemcpy(bswvalue,cfg->betasw,sizeof(float)*cfg->n_dust,cudaMemcpyHostToDevice);
    cudaGetSymbolAddress((void **)&qvalue,q); // get a pointer to size
    cudaMemcpy(qvalue,cfg->q,sizeof(float)*cfg->n_dust,cudaMemcpyHostToDevice);

//    cudaMemcpyToSymbol(size,cfg->s,sizeof(double)*cfg->n_dust,0,cudaMemcpyHostToDevice);
//    cudaMemcpyToSymbol(beta,cfg->beta,sizeof(float)*cfg->n_dust,0,cudaMemcpyHostToDevice);
//    cudaMemcpyToSymbol(betasw,cfg->betasw,sizeof(float)*cfg->n_dust,0,cudaMemcpyHostToDevice);
//    cudaMemcpyToSymbol(q,cfg->q,sizeof(float)*cfg->n_dust,0,cudaMemcpyHostToDevice);

}


void Grid::sync_from_gpu(){
    if (!synced ){

    gpuErrchk(cudaMemcpy(host.x   ,gpu.x   ,sizeof(double)*n_all,cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(host.y   ,gpu.y   ,sizeof(double)*n_all,cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(host.z   ,gpu.z   ,sizeof(double)*n_all,cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(host.vx  ,gpu.vx  ,sizeof(double)*n_all,cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(host.vy  ,gpu.vy  ,sizeof(double)*n_all,cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(host.vz  ,gpu.vz  ,sizeof(double)*n_all,cudaMemcpyDeviceToHost));

    synced = true;
    }
}
