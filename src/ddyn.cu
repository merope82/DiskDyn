/*****************************************************************************/
/**                             DiskDyn                                     **/
/**                                                                         **/
/** Code to model the 3D dynamical evolution of debris disks.               **/
/** Written for systems with CUDA cards with compute capabilities >= 2.1.   **/
/*****************************************************************************/
/** Release version & date:                                                 **/
/**/ const char *Version = "0.01";                                         /**/
/**/ const char *VDate   = "10.12.2018";                                   /**/
/** Authors: Andras Gaspar (agaspar@as.arizona.edu)                         **/
/*****************************************************************************/
/** More info: http://merope.as.arizona.edu/~agaspar/DiskDyn                **/
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ddyn.cuh"

configdata	cconfig;
cudadata	ccuda;

__constant__ double size[nsize];	// If only these two are in constant 
__constant__ double beta[nsize];	// memory, this is the max size!
__constant__ double    q[nsize];	// memory, this is the max size!
__device__   double doubletmp;		// Temporary device double value

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
    
    // Evolve system
    grid->evolve(Version);    
    printdone();
 
    return(0);
}
