/*****************************************************************************/
/**                             DiskDyn                                     **/
/**                                                                         **/
/** Code to model the 3D dynamical evolution of debris disks.               **/
/** Written for systems with CUDA cards with compute capabilities >= 2.1.   **/
/*****************************************************************************/
/** Release version & date:                                                 **/
/**/ const char *Version = "1.01";                                         /**/
/**/ const char *VDate   = "4.23.2025";                                    /**/
/** Authors: Andras Gaspar (agaspar@as.arizona.edu)                         **/
/*****************************************************************************/
/** More info: http://merope.as.arizona.edu/~agaspar/DiskDyn                **/
/*****************************************************************************/
/* Updates                                                                  **/
/* - Memory handling for newer GPU architectures                            **/
/* - Stellar Wind drag properly handled                                     **/
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ddyn.cuh"

configdata	cconfig;
cudadata	ccuda;

// trying some of these as floats instead of doubles for memory conservation in const
__constant__ double size[nsize];	// Size is kept double
__constant__ float  beta[nsize];	// beta can be float
__constant__ float  betasw[nsize];	// betasw can be float
__constant__ float  q[nsize];		// q can be float
__device__   double doubletmp;		// temporary double mem

/*****************************************************************************/

int main(int argc,char *argv[])
{
    configdata		*cfg = &cconfig;	// Array of configuration data
    Grid		*grid;			// Grid of cells of dust data

    // Read in config file on server node (or just single computer if non-mpi)
    int  		s=0;
    char 		*model,*param;

    model=param=NULL;

    if ( argc==1 )                            exit_with_usage(0);
    for ( int i=1 ; i<argc ; i++ ){
	if ( strcmp(argv[i],"-h")==0 )
            exit_with_usage(0);
	else if ( strcmp(argv[i],"-i")==0 ){
            i++;if ( i==argc )                exit_with_usage(1);
	    param=argv[i]; s=1; break;
	}
	else if ( strcmp(argv[i],"-m")==0 ){
    	    i++;if ( i==argc )                exit_with_usage(2);
    	    model=argv[i]; s=2; break;
	}
	else                                  exit_with_usage(0);
    }

    cfg->model = s;

    // Read in model either from parameter file or from existing run file
    switch( s ){
	case 1: read_in_param(param); break;
	case 2: read_in_model(model,Version); break;
    }

    // Print header to stdout 
    print_header(Version,VDate);

    // Get properties of the CUDA cards that are present
    cudaprops();

    // Assign CUDA cards to use. Exit if not available or not enough.
    // Identify compute capabilities and CUDA version. Exit, if less than necessary.
    print_cuda();
    checkcuda();

    // Add optical data from input file
    print_optics();
    add_optical_data();

    // Determine particle size scaling - only for param file
    if ( s==1 ){
	detminmax();
	confCdisks();
	confCblobs();
    }

    // Read in stellar Kurucz file
    add_specfile();

    // Set up particle distribution
    print_setup();
    grid = new Grid();

    switch( s ){
	case 1: grid->setup(); break;
	case 2: grid->read_in(model); break;
    }

    fflush(stdout);
    // Evolve system
    grid->evolve(Version);
    printdone();
 
    return(0);
}
