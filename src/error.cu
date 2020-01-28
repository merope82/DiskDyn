#include <stdio.h>
#include <stdlib.h>
#include "ddyn.cuh"

const char *errmsg[] = { "",
/* 1 */ "please give input parameter filename,",
/* 2 */ "please give input model filename,",
/* 3 */ "missing/empty input/model file,",
/* 4 */ "please give a valid verbose parameter",
/* 5 */ "please give yes or no to auto gpu parameter",
/* 6 */ "number of GPUs requested is NaN or < 1",
/* 7 */ "graphics step size needs to be > 0",
/* 8 */ "please give yes or no to graphics parameter",
/* 9 */ "write out step size needs to be > 0",
/* 10 */ "simulation end time needs to be > 0",
/* 11 */ "number of steps per orbit needs to be > 0",
/* 12 */ "stellar parameters need to have values > 0",
/* 13 */ "program requires all stellar parameters, even if not used",
/* 14 */ "bulk density of particles needs to be > 0",
/* 15 */ "please provide all three variables for material composition",
/* 16 */ "please provide grid mass ratio in param file",\
/* 17 */ "problem with parameter file, not enough parameters",\
/* 18 */ "not enough variables for planetary orbit in param file",\
/* 19 */ "cell log spacing not defined in param file",\
/* 20 */ "Requested number of GPUs not available",\
/* 21 */ "GPU compute capability too small",\
/* 22 */ "Need cuda version 4.0 or higher for multi-GPU",\
/* 23 */ "Rank cannot be other than 0 for non-MPI version",\
/* 24 */ "cell Rin larger than disk parameter",\
/* 25 */ "cell Rout smaller than disk parameter",\
/* 26 */ "Number of rings has to be integer multiplicative of number of GPUs",\
/* 27 */ "cell dhdr smaller than disk parameter",\
/* 28 */ "Max variation needs to be between 0 and 1",\
/* 29 */ "The number of tracers per mass needs to be either 1,2,4,8, or 16",\
/* 30 */ "Need to answer all Physics to include parameters",\
/* 31 */ "Answer yes/no to physics included parameters",\
/* 32 */ "Electric charge coefficient needs to be >0",\
/* 33 */ "Subgridding parameter needs to be >=0",\
/* 34 */ "Overdensity parameter needs to be >=0",\
/* 35 */ "Cuda card requested not available",\
/* 36 */ "Missing a System geometry parameter",\
/* 37 */ "Negative distance",\
/* 38 */ "Negative FOV",\
/* 39 */ "Spectral file missing",\
/* 40 */ "Teff out of Library temperature range",\
/* 41 */ "Need to give wavelengths!",\
/* 42 */ "Optical data not well formed!",\
/* 43 */ "Requested image wavelength not possible!",\
/* 44 */ "RAM allocation for planets failed!\n",\
/* 45 */ "RAM allocation for disks failed!\n",\
/* 46 */ "RAM allocation for blobs failed!\n",\
/* 47 */ "RAM allocation for spectra failed!\n",\
/* 48 */ "RAM allocation for optical constants failed!\n",\
/* 49 */ "Out of memory!\n",\
/* 50 */ "RAM allocation for new optical constants failed!\n",\
/* 51 */ "RAM allocation for SED output failed!\n",\
/* 52 */ "RAM allocation for betas and Qprs failed!\n",\
/* 53 */ "RAM allocation for particles failed!\n",\
/* 54 */ "Interpolation error!\n",\
/* 55 */ "Read-in file version is different from DiskDyn version!\n",\
/* 56 */ "nR value used in original compilation is different!\n",\
/* 57 */ "Please give cloud or collision parameter for blob!\n",\
/* 58 */ "Not enough blob parameters!\n",\
	NULL }; 

void print_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\t./ddyn [-h] {help} [-i] <init.param> {initializing\n"
"\tparameter file} [-m] <model file> {model file to continue}\n");

 return;
}

void exit_with_usage(int t)
{
 if ( t>0  )         fprintf(stderr,"\nError:\t%s\n\n",errmsg[t]);
 if ( t<=2 )         print_usage(stderr);
 if ( t )            exit(1);
 else                exit(0);
}

void Cuda_Error( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( 1 );
    }
}
