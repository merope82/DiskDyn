#include <stdio.h>
#include "ddyn.cuh"

// Function taken from NVIDIA CUDA helper files
inline int ConvertSMVer2Cores(int major, int minor)
{
    configdata	*cfg = &cconfig;
    // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
    typedef struct
    {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] =
    {
        { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
        { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
        { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
        { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
        { 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
        { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
        { 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
        {   -1, -1 }
    };

    int index = 0;

    while (nGpuArchCoresPerSM[index].SM != -1)
    {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
        {
            return nGpuArchCoresPerSM[index].Cores;
        }

        index++;
    }

    // If we don't find the values, we default use the previous one to run properly
    if ( cfg->verb != 0 )
    printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[index-1].Cores);
    return nGpuArchCoresPerSM[index-1].Cores;
}

// Get properties of CUDA cards
int cudaprops(){
    configdata		*cfg = &cconfig;
    cudaDeviceProp	prop;
    cudadata		*cuda= &ccuda;
    int		      	mycount;
								// Number of GPU devices on node
    CudaError(cudaGetDeviceCount(&mycount)); 			// may not be what actually is used!

    cuda->ncuda=0;
    cuda->prop=(cudaconf*)malloc(sizeof(cudaconf)*mycount);
    cfg->mem += sizeof(cudaconf)*mycount;
    if ( !cuda->prop ) exit_with_usage(49);

    for (int i=0 ; i<mycount ; i++){
	cudaSetDevice(i);
	cudaDeviceReset();
	cudaDeviceSynchronize();
//	cudaThreadSynchronize();
        CudaError(cudaGetDeviceProperties(&prop,i));
	int length=strlen(prop.name);
	cuda->prop[i].name = (char *)malloc(length+1);
	cfg->mem += length+1;
	if ( !cuda->prop[i].name ) exit_with_usage(49);
	strcpy(cuda->prop[i].name,prop.name);
	cuda->prop[i].gpu=i;
	cuda->prop[i].sm = prop.multiProcessorCount;
	cuda->prop[i].cores=ConvertSMVer2Cores(prop.major, prop.minor) * prop.multiProcessorCount;
	cuda->prop[i].major=prop.major;
	cuda->prop[i].minor=prop.minor;
	cuda->prop[i].totalGlobalMem=prop.totalGlobalMem;
	cuda->prop[i].totalConstMem=prop.totalConstMem;
	cuda->prop[i].memPitch=prop.memPitch;
	cuda->prop[i].sharedMemPerBlock=prop.sharedMemPerBlock;
	cuda->prop[i].maxThreadsPerBlock=prop.maxThreadsPerBlock;
	cuda->ncuda++;
    }

    return 0;
}

// Check CUDA card
void checkcuda(){
    configdata	*cfg = &cconfig;
    cudadata	*cuda = &ccuda;

    if ( cuda->prop[cfg->gpuid].major*100 + cuda->prop[cfg->gpuid].minor*10 < compute ) exit_with_usage(21);
    cudaSetDevice(cfg->gpuid);
}
