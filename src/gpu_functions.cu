#include <stdio.h>
#include "ddyn.cuh"
#include "ddyn.cuh"

__device__ double atomicAdd_double(double* address, double val){
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__global__ void null_coord_driver(int n_all,coord gpu){
    int t = threadIdx.x + blockIdx.x * blockDim.x;
    
    if ( t < n_all ){
	gpu.x[t]    = 0.0;
        gpu.y[t]    = 0.0;
        gpu.z[t]    = 0.0;
        gpu.vx[t]   = 0.0;
        gpu.vy[t]   = 0.0;
        gpu.vz[t]   = 0.0;
        gpu.m [t]   = 0.0;
        gpu.N[t]    = 0.0;

	gpu.k1x[t]  = 0.0;
	gpu.k2x[t]  = 0.0;
	gpu.k3x[t]  = 0.0;
	gpu.k4x[t]  = 0.0;
	gpu.k1y[t]  = 0.0;
	gpu.k2y[t]  = 0.0;
	gpu.k3y[t]  = 0.0;
	gpu.k4y[t]  = 0.0;
	gpu.k1z[t]  = 0.0;
	gpu.k2z[t]  = 0.0;
	gpu.k3z[t]  = 0.0;
	gpu.k4z[t]  = 0.0;

	gpu.k2rx[t]  = 0.0;
	gpu.k3rx[t]  = 0.0;
	gpu.k4rx[t]  = 0.0;
	gpu.k2ry[t]  = 0.0;
	gpu.k3ry[t]  = 0.0;
	gpu.k4ry[t]  = 0.0;
	gpu.k2rz[t]  = 0.0;
	gpu.k3rz[t]  = 0.0;
	gpu.k4rz[t]  = 0.0;

	gpu.k1vx[t]  = 0.0;
	gpu.k2vx[t]  = 0.0;
	gpu.k3vx[t]  = 0.0;
	gpu.k4vx[t]  = 0.0;
	gpu.k1vy[t]  = 0.0;
	gpu.k2vy[t]  = 0.0;
	gpu.k3vy[t]  = 0.0;
	gpu.k4vy[t]  = 0.0;
	gpu.k1vz[t]  = 0.0;
	gpu.k2vz[t]  = 0.0;
	gpu.k3vz[t]  = 0.0;
	gpu.k4vz[t]  = 0.0;
    }
}

// DYNAMICAL CALCULATIONS

// Will need q in (A yr) instead of (A s), Wz in (1/year) instead of (1/s), B in Msun/year^2/A instead of kg/s^2/A, Rstar in AU 
// Divide kq by "year", multiply Wz by "year", divide B by Msun/year^2, divide Rstar by AU (if in SI) - done

__global__ void RK4calc_massive0(int n_grav,coord gpu){
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if ( i < n_grav ){
        gpu.k1vx[i]=0.0;
	gpu.k1vy[i]=0.0;
	gpu.k1vz[i]=0.0; 
        gpu.k2vx[i]=0.0;
	gpu.k2vy[i]=0.0;
	gpu.k2vz[i]=0.0; 
        gpu.k3vx[i]=0.0;
	gpu.k3vy[i]=0.0;
	gpu.k3vz[i]=0.0;
        gpu.k4vx[i]=0.0;
	gpu.k4vy[i]=0.0;
	gpu.k4vz[i]=0.0; 
    }
}

__global__ void RK4calc_massive1(int n_grav,coord gpu){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    if ( i < n_grav ){
	// K1
        gpu.k1x[i]  = gpu.x[i];
        gpu.k1y[i]  = gpu.y[i];
        gpu.k1z[i]  = gpu.z[i];
    }
}

__global__ void RK4calc_massive2(int j0,int ngrav,int TILE,coord gpu,int n_planets){
    int sq = TILE*TILE;
    extern __shared__ double s[];
    double *kvx = s;
    double *kvy = &kvx[sq];
    double *kvz = &kvy[sq];
    struct array *tmparr = (array*) &kvz[sq];

    int t = threadIdx.x;

    int jstart = (j0+blockIdx.x)*TILE;
    int istart = blockIdx.x*TILE;

    int ig = istart+t;
    int jg = jstart+t-TILE;

    // Let's populate tmparr
    if ( t < TILE && ig < ngrav ){
	tmparr[t].x   = gpu.k1x[ig];
	tmparr[t].y   = gpu.k1y[ig];
	tmparr[t].z   = gpu.k1z[ig];
	tmparr[t].m   = gpu.m[ig];
    }
    else if ( t>=TILE && t<2*TILE && jg < ngrav ){
	tmparr[t].x   = gpu.k1x[jg];
	tmparr[t].y   = gpu.k1y[jg];
	tmparr[t].z   = gpu.k1z[jg];
	tmparr[t].m   = gpu.m[jg];
    }
    else if ( t<TILE && ig >= ngrav ){
	tmparr[t].x   = 0.0;
	tmparr[t].y   = 0.0;
	tmparr[t].z   = 0.0;
	tmparr[t].m   = 0.0;
    }
    else if ( t>=TILE && t<2*TILE && jg >= ngrav ){
	tmparr[t].x   = 0.0;
	tmparr[t].y   = 0.0;
	tmparr[t].z   = 0.0;
	tmparr[t].m   = 0.0;
    }
    else;
    
    __syncthreads();

    // Let's calculate the multiplication array

    int tj = int(t/TILE);	// local tile i,j index
    int ti = t-tj*TILE;

    ig = istart+ti;
    jg = jstart+tj;

    tj += TILE;

    if ( jg > ig && jg < ngrav && ig < ngrav ){
	double dX = tmparr[tj].x - tmparr[ti].x;
	double dY = tmparr[tj].y - tmparr[ti].y;
	double dZ = tmparr[tj].z - tmparr[ti].z;
	double RR = dX * dX + dY * dY + dZ * dZ;

        double f  = RR * sqrt(RR);
	if ( f!=0 ){
    	    kvx[t]    =  dX / f;
    	    kvy[t]    =  dY / f;
    	    kvz[t]    =  dZ / f;
	}
	else{
	    kvx[t] = 0.0;
	    kvy[t] = 0.0;
	    kvz[t] = 0.0;
	}
    }
    else{
	kvx[t] = 0.0;
	kvy[t] = 0.0;
	kvz[t] = 0.0;
    }
    
    __syncthreads();

    if ( t < TILE ){
	ig = istart+t;
	if ( ig < ngrav ){
	    double KX=0.0;
	    double KY=0.0;
	    double KZ=0.0;
	    for ( int j=t ; j<TILE*TILE ; j+=TILE ){
		tj = int(j/TILE)+TILE;
		KX += kvx[j] * tmparr[tj].m;
		KY += kvy[j] * tmparr[tj].m;
		KZ += kvz[j] * tmparr[tj].m;
	    }
    	    atomicAdd_double(&gpu.k1vx[ig],KX);
    	    atomicAdd_double(&gpu.k1vy[ig],KY);
    	    atomicAdd_double(&gpu.k1vz[ig],KZ);
	}
    }
    else if ( t>=TILE && t<2*TILE ){
	jg = jstart+t-TILE;
	if ( jg < ngrav ){
	    double KX=0.0;
	    double KY=0.0;
	    double KZ=0.0;
	    for ( int i=(t-TILE)*TILE ; i<(t-TILE+1)*TILE ; i++ ){
		ti = i - (t-TILE)*TILE;
		KX -= kvx[i] * tmparr[ti].m;
		KY -= kvy[i] * tmparr[ti].m;
		KZ -= kvz[i] * tmparr[ti].m;
	    }
    	    atomicAdd_double(&gpu.k1vx[jg],KX);
    	    atomicAdd_double(&gpu.k1vy[jg],KY);
    	    atomicAdd_double(&gpu.k1vz[jg],KZ);
	}
    }
    else ;
}

__global__ void RK4calc_massive3(int n_grav,coord gpu,double dt){ // dt=0.5 dt in call function
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    if ( i < n_grav ){ // Only evolving massive and empty tracers. Phantoms will get regenerated before collisions
        gpu.k2x[i]  = gpu.x[i]  + gpu.vx[i]   * dt;
        gpu.k2y[i]  = gpu.y[i]  + gpu.vy[i]   * dt;
        gpu.k2z[i]  = gpu.z[i]  + gpu.vz[i]   * dt;
        gpu.k2rx[i] = gpu.vx[i] + gpu.k1vx[i] * dt;
        gpu.k2ry[i] = gpu.vy[i] + gpu.k1vy[i] * dt;
        gpu.k2rz[i] = gpu.vz[i] + gpu.k1vz[i] * dt;
    }
}

__global__ void RK4calc_massive4(int j0,int ngrav,int TILE,coord gpu,int n_planets){
    int sq = TILE*TILE;
    extern __shared__ double s[];
    double *kvx = s;
    double *kvy = &kvx[sq];
    double *kvz = &kvy[sq];
    struct array *tmparr = (array*) &kvz[sq];

    int t = threadIdx.x;

    int jstart = (j0+blockIdx.x)*TILE;
    int istart = blockIdx.x*TILE;

    int ig = istart+t;
    int jg = jstart+t-TILE;

    // Let's populate tmparr
    if ( t < TILE && ig < ngrav ){
	tmparr[t].x   = gpu.k2x[ig];
	tmparr[t].y   = gpu.k2y[ig];
	tmparr[t].z   = gpu.k2z[ig];
	tmparr[t].m   = gpu.m[ig];
    }
    else if ( t>=TILE && t<2*TILE && jg < ngrav ){
	tmparr[t].x   = gpu.k2x[jg];
	tmparr[t].y   = gpu.k2y[jg];
	tmparr[t].z   = gpu.k2z[jg];
	tmparr[t].m   = gpu.m[jg];
    }
    else if ( t<TILE && ig >= ngrav ){
	tmparr[t].x   = 0.0;
	tmparr[t].y   = 0.0;
	tmparr[t].z   = 0.0;
	tmparr[t].m   = 0.0;
    }
    else if ( t>=TILE && t<2*TILE && jg >= ngrav ){
	tmparr[t].x   = 0.0;
	tmparr[t].y   = 0.0;
	tmparr[t].z   = 0.0;
	tmparr[t].m   = 0.0;
    }
    else;
    
    __syncthreads();

    // Let's calculate the multiplication array

    int tj = int(t/TILE);	// local tile i,j index
    int ti = t-tj*TILE;

    ig = istart+ti;
    jg = jstart+tj;

    tj += TILE;

    if ( jg > ig && jg < ngrav && ig < ngrav ){
	double dX = tmparr[tj].x - tmparr[ti].x;
	double dY = tmparr[tj].y - tmparr[ti].y;
	double dZ = tmparr[tj].z - tmparr[ti].z;
	double RR = dX * dX + dY * dY + dZ * dZ;

        double f  = RR * sqrt(RR);
	if ( f!=0 ){
    	    kvx[t]    =  dX / f;
    	    kvy[t]    =  dY / f;
    	    kvz[t]    =  dZ / f;
	}
	else{
	    kvx[t] = 0.0;
	    kvy[t] = 0.0;
	    kvz[t] = 0.0;
	}
    }
    else{
	kvx[t] = 0.0;
	kvy[t] = 0.0;
	kvz[t] = 0.0;
    }
    
    __syncthreads();

    if ( t < TILE ){
	ig = istart+t;
	if ( ig < ngrav ){
	    double KX=0.0;
	    double KY=0.0;
	    double KZ=0.0;
	    for ( int j=t ; j<TILE*TILE ; j+=TILE ){
		tj = int(j/TILE)+TILE;
		KX += kvx[j] * tmparr[tj].m;
		KY += kvy[j] * tmparr[tj].m;
		KZ += kvz[j] * tmparr[tj].m;
	    }
    	    atomicAdd_double(&gpu.k2vx[ig],KX);
    	    atomicAdd_double(&gpu.k2vy[ig],KY);
    	    atomicAdd_double(&gpu.k2vz[ig],KZ);
	}
    }
    else if ( t>=TILE && t<2*TILE ){
	jg = jstart+t-TILE;
	if ( jg < ngrav ){
	    double KX=0.0;
	    double KY=0.0;
	    double KZ=0.0;
	    for ( int i=(t-TILE)*TILE ; i<(t-TILE+1)*TILE ; i++ ){
		ti = i - (t-TILE)*TILE;
		KX -= kvx[i] * tmparr[ti].m;
		KY -= kvy[i] * tmparr[ti].m;
		KZ -= kvz[i] * tmparr[ti].m;
	    }
    	    atomicAdd_double(&gpu.k2vx[jg],KX);
    	    atomicAdd_double(&gpu.k2vy[jg],KY);
    	    atomicAdd_double(&gpu.k2vz[jg],KZ);
	}
    }
    else ;
}

__global__ void RK4calc_massive5(int n_grav,coord gpu,double dt){ // dt=0.5 dt in call function
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    if ( i < n_grav ){ // Only evolving massive and empty tracers. Phantoms will get regenerated before collisions
	// K3
        gpu.k3x[i]  = gpu.x[i]  + gpu.k2rx[i] * dt;
        gpu.k3y[i]  = gpu.y[i]  + gpu.k2ry[i] * dt;
        gpu.k3z[i]  = gpu.z[i]  + gpu.k2rz[i] * dt;
        gpu.k3rx[i] = gpu.vx[i] + gpu.k2vx[i] * dt;
        gpu.k3ry[i] = gpu.vy[i] + gpu.k2vy[i] * dt;
        gpu.k3rz[i] = gpu.vz[i] + gpu.k2vz[i] * dt;
    }
}

__global__ void RK4calc_massive6(int j0,int ngrav,int TILE,coord gpu,int n_planets){
    int sq = TILE*TILE;
    extern __shared__ double s[];
    double *kvx = s;
    double *kvy = &kvx[sq];
    double *kvz = &kvy[sq];
    struct array *tmparr = (array*) &kvz[sq];

    int t = threadIdx.x;

    int jstart = (j0+blockIdx.x)*TILE;
    int istart = blockIdx.x*TILE;

    int ig = istart+t;
    int jg = jstart+t-TILE;

    // Let's populate tmparr
    if ( t < TILE && ig < ngrav ){
	tmparr[t].x   = gpu.k3x[ig];
	tmparr[t].y   = gpu.k3y[ig];
	tmparr[t].z   = gpu.k3z[ig];
	tmparr[t].m   = gpu.m[ig];
    }
    else if ( t>=TILE && t<2*TILE && jg < ngrav ){
	tmparr[t].x   = gpu.k3x[jg];
	tmparr[t].y   = gpu.k3y[jg];
	tmparr[t].z   = gpu.k3z[jg];
	tmparr[t].m   = gpu.m[jg];
    }
    else if ( t<TILE && ig >= ngrav ){
	tmparr[t].x   = 0.0;
	tmparr[t].y   = 0.0;
	tmparr[t].z   = 0.0;
	tmparr[t].m   = 0.0;
    }
    else if ( t>=TILE && t<2*TILE && jg >= ngrav ){
	tmparr[t].x   = 0.0;
	tmparr[t].y   = 0.0;
	tmparr[t].z   = 0.0;
	tmparr[t].m   = 0.0;
    }
    else;
    
    __syncthreads();

    // Let's calculate the multiplication array

    int tj = int(t/TILE);	// local tile i,j index
    int ti = t-tj*TILE;

    ig = istart+ti;
    jg = jstart+tj;

    tj += TILE;

    if ( jg > ig && jg < ngrav && ig < ngrav ){
	double dX = tmparr[tj].x - tmparr[ti].x;
	double dY = tmparr[tj].y - tmparr[ti].y;
	double dZ = tmparr[tj].z - tmparr[ti].z;
	double RR = dX * dX + dY * dY + dZ * dZ;

        double f  = RR * sqrt(RR);
	if ( f!=0 ){
    	    kvx[t]    =  dX / f;
    	    kvy[t]    =  dY / f;
    	    kvz[t]    =  dZ / f;
	}
	else{
	    kvx[t] = 0.0;
	    kvy[t] = 0.0;
	    kvz[t] = 0.0;
	}
    }
    else{
	kvx[t] = 0.0;
	kvy[t] = 0.0;
	kvz[t] = 0.0;
    }
    
    __syncthreads();

    if ( t < TILE ){
	ig = istart+t;
	if ( ig < ngrav ){
	    double KX=0.0;
	    double KY=0.0;
	    double KZ=0.0;
	    for ( int j=t ; j<TILE*TILE ; j+=TILE ){
		tj = int(j/TILE)+TILE;
		KX += kvx[j] * tmparr[tj].m;
		KY += kvy[j] * tmparr[tj].m;
		KZ += kvz[j] * tmparr[tj].m;
	    }
    	    atomicAdd_double(&gpu.k3vx[ig],KX);
    	    atomicAdd_double(&gpu.k3vy[ig],KY);
    	    atomicAdd_double(&gpu.k3vz[ig],KZ);
	}
    }
    else if ( t>=TILE && t<2*TILE ){
	jg = jstart+t-TILE;
	if ( jg < ngrav ){
	    double KX=0.0;
	    double KY=0.0;
	    double KZ=0.0;
	    for ( int i=(t-TILE)*TILE ; i<(t-TILE+1)*TILE ; i++ ){
		ti = i - (t-TILE)*TILE;
		KX -= kvx[i] * tmparr[ti].m;
		KY -= kvy[i] * tmparr[ti].m;
		KZ -= kvz[i] * tmparr[ti].m;
	    }
    	    atomicAdd_double(&gpu.k3vx[jg],KX);
    	    atomicAdd_double(&gpu.k3vy[jg],KY);
    	    atomicAdd_double(&gpu.k3vz[jg],KZ);
	}
    }
    else ;
}

__global__ void RK4calc_massive7(int n_grav,coord gpu,double dt){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    if ( i < n_grav ){
	gpu.k4x[i]  = gpu.x[i]  + gpu.k3rx[i] * dt;
	gpu.k4y[i]  = gpu.y[i]  + gpu.k3ry[i] * dt;
	gpu.k4z[i]  = gpu.z[i]  + gpu.k3rz[i] * dt;
	gpu.k4rx[i] = gpu.vx[i] + gpu.k3vx[i] * dt;
	gpu.k4ry[i] = gpu.vy[i] + gpu.k3vy[i] * dt;
	gpu.k4rz[i] = gpu.vz[i] + gpu.k3vz[i] * dt;
    }
}

__global__ void RK4calc_massive8(int j0,int ngrav,int TILE,coord gpu,int n_planets){
    int sq = TILE*TILE;
    extern __shared__ double s[];
    double *kvx = s;
    double *kvy = &kvx[sq];
    double *kvz = &kvy[sq];
    struct array *tmparr = (array*) &kvz[sq];

    int t = threadIdx.x;

    int jstart = (j0+blockIdx.x)*TILE;
    int istart = blockIdx.x*TILE;

    int ig = istart+t;
    int jg = jstart+t-TILE;

    // Let's populate tmparr
    if ( t < TILE && ig < ngrav ){
	tmparr[t].x   = gpu.k4x[ig];
	tmparr[t].y   = gpu.k4y[ig];
	tmparr[t].z   = gpu.k4z[ig];
	tmparr[t].m   = gpu.m[ig];
    }
    else if ( t>=TILE && t<2*TILE && jg < ngrav ){
	tmparr[t].x   = gpu.k4x[jg];
	tmparr[t].y   = gpu.k4y[jg];
	tmparr[t].z   = gpu.k4z[jg];
	tmparr[t].m   = gpu.m[jg];
    }
    else if ( t<TILE && ig >= ngrav ){
	tmparr[t].x   = 0.0;
	tmparr[t].y   = 0.0;
	tmparr[t].z   = 0.0;
	tmparr[t].m   = 0.0;
    }
    else if ( t>=TILE && t<2*TILE && jg >= ngrav ){
	tmparr[t].x   = 0.0;
	tmparr[t].y   = 0.0;
	tmparr[t].z   = 0.0;
	tmparr[t].m   = 0.0;
    }
    else;
    
    __syncthreads();

    // Let's calculate the multiplication array

    int tj = int(t/TILE);	// local tile i,j index
    int ti = t-tj*TILE;

    ig = istart+ti;
    jg = jstart+tj;

    tj += TILE;

    if ( jg > ig && jg < ngrav && ig < ngrav ){
	double dX = tmparr[tj].x - tmparr[ti].x;
	double dY = tmparr[tj].y - tmparr[ti].y;
	double dZ = tmparr[tj].z - tmparr[ti].z;
	double RR = dX * dX + dY * dY + dZ * dZ;

        double f  = RR * sqrt(RR);
	if ( f!=0 ){
    	    kvx[t]    =  dX / f;
    	    kvy[t]    =  dY / f;
    	    kvz[t]    =  dZ / f;
	}
	else{
	    kvx[t] = 0.0;
	    kvy[t] = 0.0;
	    kvz[t] = 0.0;
	}
    }
    else{
	kvx[t] = 0.0;
	kvy[t] = 0.0;
	kvz[t] = 0.0;
    }
    
    __syncthreads();

    if ( t < TILE ){
	ig = istart+t;
	if ( ig < ngrav ){
	    double KX=0.0;
	    double KY=0.0;
	    double KZ=0.0;
	    for ( int j=t ; j<sq ; j+=TILE ){
		tj = int(j/TILE)+TILE;
		KX += kvx[j] * tmparr[tj].m;
		KY += kvy[j] * tmparr[tj].m;
		KZ += kvz[j] * tmparr[tj].m;
	    }
    	    atomicAdd_double(&gpu.k4vx[ig],KX);
    	    atomicAdd_double(&gpu.k4vy[ig],KY);
    	    atomicAdd_double(&gpu.k4vz[ig],KZ);
	}
    }
    else if ( t>=TILE && t<2*TILE ){
	jg = jstart+t-TILE;
	if ( jg < ngrav ){
	    double KX=0.0;
	    double KY=0.0;
	    double KZ=0.0;
	    for ( int i=(t-TILE)*TILE ; i<(t-TILE+1)*TILE ; i++ ){
		ti = i - (t-TILE)*TILE;
		KX -= kvx[i] * tmparr[ti].m;
		KY -= kvy[i] * tmparr[ti].m;
		KZ -= kvz[i] * tmparr[ti].m;
	    }
    	    atomicAdd_double(&gpu.k4vx[jg],KX);
    	    atomicAdd_double(&gpu.k4vy[jg],KY);
    	    atomicAdd_double(&gpu.k4vz[jg],KZ);
	}
    }
    else ;
}

__global__ void RK4calc_massive_final(int n_grav,coord gpu,double dt){ //dt/6 in call function
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if ( i < n_grav ){
	gpu.x[i]  += dt * ( gpu.vx[i]   + 2.0 * gpu.k2rx[i] + 2.0 * gpu.k3rx[i] + gpu.k4rx[i] );
	gpu.y[i]  += dt * ( gpu.vy[i]   + 2.0 * gpu.k2ry[i] + 2.0 * gpu.k3ry[i] + gpu.k4ry[i] );
	gpu.z[i]  += dt * ( gpu.vz[i]   + 2.0 * gpu.k2rz[i] + 2.0 * gpu.k3rz[i] + gpu.k4rz[i] );

	gpu.vx[i] += dt * ( gpu.k1vx[i] + 2.0 * gpu.k2vx[i] + 2.0 * gpu.k3vx[i] + gpu.k4vx[i] );
	gpu.vy[i] += dt * ( gpu.k1vy[i] + 2.0 * gpu.k2vy[i] + 2.0 * gpu.k3vy[i] + gpu.k4vy[i] );
	gpu.vz[i] += dt * ( gpu.k1vz[i] + 2.0 * gpu.k2vz[i] + 2.0 * gpu.k3vz[i] + gpu.k4vz[i] );
    }
}

__global__ void RK4calc_dust(int n_grav,int n_dust,coord gpu,double dt,int FLQ,double Wz,double Bstar,double Rstar,double rho,int SWQ,double vswAU,double Rmin){

    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if ( i < n_dust ){
	double k1vx,k1vy,k1vz;
	double k2vx,k2vy,k2vz,k2rx,k2ry,k2rz;
	double k3vx,k3vy,k3vz,k3rx,k3ry,k3rz;
	double k4vx,k4vy,k4vz,k4rx,k4ry,k4rz;
	double dt2,dt6,RR,RR3,m,Fg;
	double xt,yt,zt,x,y,z,vx,vy,vz,dx,dy,dz,f;
	double Fbase1,Fbase2,Fdrx,Fdry,Fdrz,Flx,Fly;
	float  bet,betsw,qe;

	int t = n_grav + i;

	if ( gpu.N[t] == 0 ) return;

	int id = gpu.id[t];

	dt2   = dt / 2.0;
	dt6   = dt / 6.0;
	x     = gpu.x[t];			// x in AU
	y     = gpu.y[t];
	z     = gpu.z[t];
	vx    = gpu.vx[t];			// vx in AU/yr
	vy    = gpu.vy[t];
	vz    = gpu.vz[t];
	if ( id <0 ){
	    bet   = 0;
	    betsw = 0;
	}
	else{
	    bet   = beta[id];
	    betsw = betasw[id];
	}

	m     = gpu.m[t];
	qe    = q[id];

	k1vx  = 0.0;
	k1vy  = 0.0;
	k1vz  = 0.0;
	k2vx  = 0.0;
	k2vy  = 0.0;
	k2vz  = 0.0;
	k3vx  = 0.0;
	k3vy  = 0.0;
	k3vz  = 0.0;
	k4vx  = 0.0;
	k4vy  = 0.0;
	k4vz  = 0.0;

	// Star
	dx    = x - gpu.k1x[0];
	dy    = y - gpu.k1y[0];
	dz    = z - gpu.k1z[0];
        RR    = dx * dx + dy * dy + dz * dz;
	RR3   = RR * sqrt(RR);
	if ( RR <= Rmin*Rmin ){
	    gpu.N[t]=0;		// Sublimated!
	    return;
	}

        Fg     = ( 1.0 - bet - betsw ) * gpu.m[0] / RR3;

//	Fbase1 = bet * gpu.m[0] / lightAU / RR;
	if ( SWQ == 1){
	    Fbase1 = gpu.m[0] / RR * ( bet / lightAU + betsw / vswAU );
	}
	else{
	    Fbase1 = bet * gpu.m[0] / lightAU / RR;
	}
	Fbase2 = ( dx * vx + dy * vy + dz * vz ) / RR;

        Fdrx  = Fbase1 * ( Fbase2 * dx + vx ); 
        Fdry  = Fbase1 * ( Fbase2 * dy + vy );
        Fdrz  = Fbase1 * ( Fbase2 * dz + vz );

	if ( FLQ == 1 ){
    	    double Bz = Bstar * Rstar * Rstar * Rstar / RR3;
	    Flx       =   qe * ( vy - Wz * dx ) * Bz / m;
	    Fly       = - qe * ( vx + Wz * dy ) * Bz / m;
	}
	else{
	    Flx     = 0;
	    Fly     = 0;
	}

	k1vx = -Fg * dx - Fdrx + Flx;
	k1vy = -Fg * dy - Fdry + Fly;
	k1vz = -Fg * dz - Fdrz;

	// Gravitating bodies
	// No magnetic fields and radiative forces calculated here!
	for ( int j = 1 ; j < n_grav ; j++ ){
	    dx    =  x - gpu.k1x[j];
	    dy    =  y - gpu.k1y[j];
	    dz    =  z - gpu.k1z[j];
    	    RR    =  dx * dx + dy * dy + dz * dz;
	    if ( RR == 0 ){
		gpu.N[t]=0;				// Accretion. Lots of zeros probably!
	        return;
	    }
	    f     =  gpu.m[j] / RR / sqrt(RR);
	    k1vx += -f * dx;
	    k1vy += -f * dy;
	    k1vz += -f * dz;
	}

	// K2
	k2rx  = vx + k1vx * dt2;
	k2ry  = vy + k1vy * dt2;
	k2rz  = vz + k1vz * dt2;
	xt    = x  +   vx * dt2;
	yt    = y  +   vy * dt2;
	zt    = z  +   vz * dt2; // good

	// Star
	dx    = xt - gpu.k2x[0];
	dy    = yt - gpu.k2y[0];
	dz    = zt - gpu.k2z[0];
        RR    = dx * dx + dy * dy + dz * dz;
	RR3   = RR * sqrt(RR);
	if ( RR <= Rmin*Rmin ){ // Sublimated
	    gpu.N[t]=0;
	    return;
	}

        Fg     = ( 1.0 - bet - betsw ) * gpu.m[0] / RR3;

	if ( SWQ == 1){
	    Fbase1 = gpu.m[0] / RR * ( bet / lightAU + betsw / vswAU );
	}
	else{
	    Fbase1 = bet * gpu.m[0] / lightAU / RR;
	}
	Fbase2 = ( dx * k2rx + dy * k2ry + dz * k2rz ) / RR;

        Fdrx  = Fbase1 * ( Fbase2 * dx + k2rx ); 
        Fdry  = Fbase1 * ( Fbase2 * dy + k2ry );
        Fdrz  = Fbase1 * ( Fbase2 * dz + k2rz );

	if ( FLQ == 1 ){
    	    double Bz = Bstar * Rstar * Rstar * Rstar / RR3;
	    Flx       =   qe * ( k2ry - Wz * dx ) * Bz / m;
	    Fly       = - qe * ( k2rx + Wz * dy ) * Bz / m;
	}
	else{
	    Flx     = 0;
	    Fly     = 0;
	}

	k2vx = -Fg * dx - Fdrx + Flx;
	k2vy = -Fg * dy - Fdry + Fly;
	k2vz = -Fg * dz - Fdrz;

	// Planets
	for ( int j = 1 ; j < n_grav ; j++ ){
	    dx    =  xt - gpu.k2x[j];
	    dy    =  yt - gpu.k2y[j];
	    dz    =  zt - gpu.k2z[j];
	    RR    =  dx * dx + dy * dy + dz * dz;
	    if ( RR == 0 ){
		gpu.N[t]=0;
	        return;
	    }
	    f     =  gpu.m[j] / RR / sqrt(RR);
	    k2vx += -f * dx;
	    k2vy += -f * dy;
	    k2vz += -f * dz;
	}

	// K3
	k3rx  = vx + k2vx * dt2;
	k3ry  = vy + k2vy * dt2;
	k3rz  = vz + k2vz * dt2;
	xt    = x  + k2rx * dt2;
	yt    = y  + k2ry * dt2;
	zt    = z  + k2rz * dt2;

	// Star
	dx    = xt - gpu.k3x[0];
	dy    = yt - gpu.k3y[0];
	dz    = zt - gpu.k3z[0];
        RR    = dx * dx + dy * dy + dz * dz;
	RR3   = RR * sqrt(RR);
	if ( RR <= Rmin*Rmin ){
	    gpu.N[t]=0;
	    return;
	}

        Fg     = ( 1.0 - bet - betsw ) * gpu.m[0] / RR3;

	if ( SWQ == 1){
	    Fbase1 = gpu.m[0] / RR * ( bet / lightAU + betsw / vswAU );
	}
	else{
	    Fbase1 = bet * gpu.m[0] / lightAU / RR;
	}
	Fbase2 = ( dx * k3rx + dy * k3ry + dz * k3rz ) / RR;

        Fdrx  = Fbase1 * ( Fbase2 * dx + k3rx ); 
        Fdry  = Fbase1 * ( Fbase2 * dy + k3ry );
        Fdrz  = Fbase1 * ( Fbase2 * dz + k3rz );

	if ( FLQ == 1 ){
    	    double Bz = Bstar * Rstar * Rstar * Rstar / RR3;
	    Flx       =   qe * ( k3ry - Wz * dx ) * Bz / m;
	    Fly       = - qe * ( k3rx + Wz * dy ) * Bz / m;
	}
	else{
	    Flx     = 0;
	    Fly     = 0;
	}

	k3vx = -Fg * dx - Fdrx + Flx;
	k3vy = -Fg * dy - Fdry + Fly;
	k3vz = -Fg * dz - Fdrz;

	// Planets
	for ( int j = 1 ; j < n_grav ; j++ ){
	    dx    =  xt - gpu.k3x[j];
	    dy    =  yt - gpu.k3y[j];
	    dz    =  zt - gpu.k3z[j];
	    RR    =  dx * dx + dy * dy + dz * dz;
	    if ( RR == 0 ){
		gpu.N[t]=0;
	        return;
	    }
	    f     =  gpu.m[j] / RR / sqrt(RR);
	    k3vx += -f * dx;
	    k3vy += -f * dy;
	    k3vz += -f * dz;
	}

	// K4
	k4rx  = vx + k3vx * dt;
	k4ry  = vy + k3vy * dt;
	k4rz  = vz + k3vz * dt;
	xt    = x  + k3rx * dt;
	yt    = y  + k3ry * dt;
	zt    = z  + k3rz * dt;

	// Star
	dx    = xt - gpu.k4x[0];
	dy    = yt - gpu.k4y[0];
	dz    = zt - gpu.k4z[0];
    	RR    = dx * dx + dy * dy + dz * dz;
	RR3   = RR * sqrt(RR);
	if ( RR <= Rmin*Rmin ){
	    gpu.N[i]=0;
	    return;
	}

        Fg     = ( 1.0 - bet - betsw ) * gpu.m[0] / RR3;

	if ( SWQ == 1){
	    Fbase1 = gpu.m[0] / RR * ( bet / lightAU + betsw / vswAU );
	}
	else{
	    Fbase1 = bet * gpu.m[0] / lightAU / RR;
	}
	Fbase2 = ( dx * k4rx + dy * k4ry + dz * k4rz ) / RR;

        Fdrx  = Fbase1 * ( Fbase2 * dx + k4rx );
        Fdry  = Fbase1 * ( Fbase2 * dy + k4ry );
        Fdrz  = Fbase1 * ( Fbase2 * dz + k4rz );

	if ( FLQ == 1 ){
    	    double Bz = Bstar * Rstar * Rstar * Rstar / RR3;
	    Flx       =   qe * ( k4ry - Wz * dx ) * Bz / m;
	    Fly       = - qe * ( k4rx + Wz * dy ) * Bz / m;
	}
	else{
	    Flx     = 0;
	    Fly     = 0;
	}

	k4vx = -Fg * dx - Fdrx + Flx;
	k4vy = -Fg * dy - Fdry + Fly;
	k4vz = -Fg * dz - Fdrz;

	// Planets
	for ( int j = 1 ; j < n_grav ; j++ ){	// First "planet" is star
	    dx    =  xt - gpu.k4x[j];
	    dy    =  yt - gpu.k4y[j];
	    dz    =  zt - gpu.k4z[j];
    	    RR    =  dx * dx + dy * dy + dz * dz;
	    if ( RR == 0 ){
		gpu.N[t]=0;
	        return;
	    }
	    f     =  gpu.m[j] / RR / sqrt(RR);
	    k4vx += -f * dx;
	    k4vy += -f * dy;
	    k4vz += -f * dz;
	}

	// Evolve
	x  += dt6 * (   vx + 2.0 * ( k2rx + k3rx ) + k4rx );
	y  += dt6 * (   vy + 2.0 * ( k2ry + k3ry ) + k4ry );
	z  += dt6 * (   vz + 2.0 * ( k2rz + k3rz ) + k4rz );

	vx += dt6 * ( k1vx + 2.0 * ( k2vx + k3vx ) + k4vx );
	vy += dt6 * ( k1vy + 2.0 * ( k2vy + k3vy ) + k4vy );
	vz += dt6 * ( k1vz + 2.0 * ( k2vz + k3vz ) + k4vz );

	gpu.x[t]  = x;
	gpu.y[t]  = y;
	gpu.z[t]  = z;

	gpu.vx[t] = vx;
	gpu.vy[t] = vy;
	gpu.vz[t] = vz;
    }
}

__global__ void null_flux_driver(int npix,double *Fd){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    if ( i < npix ){
	Fd[i] = 0.0;
    }
}

__global__ void null_flux_part_driver(int ndust,double *Fdust){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    if ( i < ndust ){
	Fdust[i] = 0.0;
    }
}

__global__ void sum_flux_driver(int ndust,double *Fdust){
    extern __shared__ double var[];
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int t = threadIdx.x;

    if ( i < ndust ){
	var[t] = Fdust[i];
    }
    else
	var[t] = 0;

    __syncthreads();
    for ( unsigned int j = blockDim.x / 2 ; j > 0 ; j >>= 1 ){
	if ( t < j ){
	    var[t] += var[t+j];
        }
	__syncthreads();
    }
    if ( t==0 ) atomicAdd_double(&doubletmp,var[0]);

}

__global__ void calc_flux_driver(int n_grav,int n_dust,coord gpu,double dist,double Rmin,double Rmax,\
                double *Qabs,double *Qsca,double *Qpfunc,double *Td,double *Fd,double *w,double *in,\
                int nwav,int ndust,int ntheta,int wj,double rot,double incl,double PA,int xsize,int ysize,double FOV){

    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if ( i < n_dust ){
	int t = n_grav + i;
	int id = gpu.id[t];

	// Check if outside of array of dust particles
	if ( id < 0 ) return;

	double x     = gpu.x[t];
	double y     = gpu.y[t];
	double z     = gpu.z[t];
	double n     = gpu.N[t];

	double phi   = atan2(y,x);
	double R     = sqrt(x*x+y*y+z*z);
	double Rd    = sqrt(x*x+y*y);
	double wvl   = w[wj];

	// Let's rotate system 
	x = Rd*cos(phi+rot);
	y = Rd*sin(phi+rot);

	// particle at x,y,z;
	// viewing location in star center coordinates: lxs=0; lys=-dist*sin(incl); lzs=dist*cos(incl)
	// viewing location in particle center coordinates: lx=-x; ly=-dist*sin(incl)-y; lz=dist*cos(incl)-z;
	// star at sx=-x; sy=-y; sz=-z;

	// cosscat = (sx*lx+sy*ly+sz*lz)/R/L;
	
	double ly = -dist*sin(incl)-y;
	double lz =  dist*cos(incl)-z;
	
	double cossc = (x*x+(-y)*ly+(-z)*lz)/R/sqrt(x*x+ly*ly+lz*lz);
	// Safety switch
	if ( cossc < -1  ) cossc = -1.0;
	if ( cossc >=  1.0-2.0/ntheta ) cossc = 1.0-2.0/ntheta;

	// Get image coordinates and check if in image
	// Let's incline the system now!
	y   = sqrt(y*y+z*z)*cos(atan2(z,y)-incl);

	// Let's rotate system by PA!
	phi = atan2(y,x);
	Rd  = sqrt(x*x+y*y);
	x   = Rd*cos(phi+PA);
	y   = Rd*sin(phi+PA);
	
	// Let's offset coordinates by half of the FOV! (FOV in AU)
	double psca;
    	if ( xsize >= ysize) psca = FOV/xsize;
	else		     psca = FOV/ysize;

	x   = x+0.5*psca*xsize;
	y   = y+0.5*psca*ysize;
    
	int xim = floor(x/psca);
	int yim = floor(y/psca);
        
	// If not in image, no need to calculate further!
	if ( xim < 0 || xim > xsize-1 || yim < 0 || yim > ysize-1 ) return;

	double f = log(Rmax/Rmin)/((double)nR-1.0);
	double T;

	if      ( R <= Rmin ){
	    T = Td[id*nR];
	}
	else if ( R >= Rmax ){
	    T = Td[id*nR+nR-1];
	}
	else{
	    int    Rid  = (int)(log(R/Rmin)/f)+1;
	    double Rin  = Rmin * exp((Rid-1)*f);
	    double Rout = Rmin * exp(Rid*f);
	    T = exp(log(Td[id*nR+Rid-1])+log(Td[id*nR+Rid]/Td[id*nR+Rid-1])*log(R/Rin)/log(Rout/Rin));
	}
	double I     = (2*hplank*cspeed*cspeed)/(pow(wvl,5)*(exp(hplank*cspeed/(wvl*kboltz*T))-1.0))*PI;	// Blambda*PI

	int cosinm   = floor(0.5*(cossc+1.0)*(ntheta-1));
	int cosinp   = cosinm+1;

	double cosm  = (double)(cosinm)*2.0/(ntheta-1)-1.0;
	double cosp  = (double)(cosinp)*2.0/(ntheta-1)-1.0;

	double P     = exp(log(Qpfunc[nwav*(ntheta*id+cosinm)+wj])+(log(Qpfunc[nwav*(ntheta*id+cosinp)+wj])-log(Qpfunc[nwav*(ntheta*id+cosinm)+wj]))*(cossc-cosm)/(cosp-cosm));

	double Fscat = in[wj]*PI*size[id]*size[id]*Qsca[nwav*id+wj]*P/R/R/AU/AU;	// Looks good
	double Fther = I*size[id]*size[id]*Qabs[nwav*id+wj]/dist/dist; 			// Looks good

	if ( R > Rmin )
    	    atomicAdd_double(&Fd[ysize*xim+yim],n*(Fther+Fscat)*1e26*wvl*wvl/cspeed);

    }
}

__global__ void calc_flux_sed_driver(int n_grav,int n_dust,coord gpu,double dist,double Rmin,double Rmax,\
                double *Qabs,double *Qsca,double *Qpfunc,double *Td,double *Fdust,double *w,double *in,\
                int nwav,int ndust,int ntheta,int wj,double rot,double incl,double PA){

    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if ( i < n_dust ){
	int t = n_grav + i;
	int id = gpu.id[t];

	// Check if outside of array of dust particles
	if ( id < 0 ) return;

	double x     = gpu.x[t];
	double y     = gpu.y[t];
	double z     = gpu.z[t];
	double n     = gpu.N[t];

	double phi   = atan2(y,x);
	double R     = sqrt(x*x+y*y+z*z);
	double wvl   = w[wj];

	// Let's rotate system 
	x = R*cos(phi+rot);
	y = R*sin(phi+rot);

	// particle at x,y,z;
	// viewing location in star center coordinates: lxs=0; lys=-dist*sin(incl); lzs=dist*cos(incl)
	// viewing location in particle center coordinates: lx=-x; ly=-dist*sin(incl)-y; lz=dist*cos(incl)-z;
	// star at sx=-x; sy=-y; sz=-z;

	// cosscat = (sx*lx+sy*ly+sz*lz)/R/L;
	double ly = -dist*sin(incl)-y;
	double lz =  dist*cos(incl)-z;

	double cossc = (x*x+(-y)*ly+(-z)*lz)/R/sqrt(x*x+ly*ly+lz*lz);
	// Safety switch
	if ( cossc < -1  ) cossc = -1.0;
	if ( cossc >=  1.0-2.0/ntheta ) cossc = 1.0-2.0/ntheta;

	double f = log(Rmax/Rmin)/((double)nR-1.0);
	double T;

	if      ( R <= Rmin ){
	    T = Td[id*nR];
	}
	else if ( R >= Rmax ){
	    T = Td[id*nR+nR-1];
	}
	else{
	    int    Rid  = (int)(log(R/Rmin)/f)+1;
	    double Rin  = Rmin * exp((Rid-1)*f);
	    double Rout = Rmin * exp(Rid*f);
	    T = exp(log(Td[id*nR+Rid-1])+log(Td[id*nR+Rid]/Td[id*nR+Rid-1])*log(R/Rin)/log(Rout/Rin));
	}
	double I     = (2*hplank*cspeed*cspeed)/(pow(wvl,5)*(exp(hplank*cspeed/(wvl*kboltz*T))-1.0))*PI;	// Blambda*PI

	int cosinm   = floor(0.5*(cossc+1.0)*(ntheta-1));
	int cosinp   = cosinm+1;

	double cosm  = (double)(cosinm)*2.0/(ntheta-1)-1.0;
	double cosp  = (double)(cosinp)*2.0/(ntheta-1)-1.0;

	double P     = exp(log(Qpfunc[nwav*(ntheta*id+cosinm)+wj])+(log(Qpfunc[nwav*(ntheta*id+cosinp)+wj])-log(Qpfunc[nwav*(ntheta*id+cosinm)+wj]))*(cossc-cosm)/(cosp-cosm));

	double Fscat = in[wj]*PI*size[id]*size[id]*Qsca[nwav*id+wj]*P/R/R/AU/AU;	// Looks good
	double Fther = I*size[id]*size[id]*Qabs[nwav*id+wj]/dist/dist; 			// Looks good

	if ( R > Rmin)
    	    Fdust[i] = n * ( Fther + Fscat ) * 1e26 *wvl * wvl / cspeed;
	else
	    Fdust[i] = 0.0;
    }
}
