#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fitsio.h"
#include "ddyn.cuh"

void Grid::write_image(int wj,const char * Version){
    fitsfile *fptr;					/* pointer to the FITS file; defined in fitsio.h */
    configdata *cfg = &cconfig;
    char *fname;
    char *ver;
    int status,i,x,y,jlen;
    long naxis = 2;
    long naxes[2];
    long fpixel = 1;
    double *F[cfg->ysize];
    status = 0; 					/* initialize status before calling fitsio routines 	*/

    gpuErrchk(cudaMemcpy(Fh,Fd,sizeof(double)*cfg->xsize*cfg->ysize,cudaMemcpyDeviceToHost));

    int whole = floor(time);
    jlen  = ( whole == 0 ? 1 : (int)(log10(whole)+1));
        whole = floor(cfg->w[wj]*1e6);
    jlen += ( whole == 0 ? 1 : (int)(log10(whole)+1));
    fname = (char *)malloc(strlen(cfg->name_stub)+jlen+16);
    if ( !fname ) exit_with_usage(49);
    sprintf(fname,"%s_%.2f_%.3f.fits",cfg->name_stub,cfg->w[wj]*1e6,time);
    
    ver=(char *)malloc(strlen(Version)+1);
    if ( !ver ) exit_with_usage(49);
    strcpy(ver,Version);

    naxes[0]=cfg->xsize;
    naxes[1]=cfg->ysize;

    remove(fname);					/* Overwrite if exists					*/

    fits_create_file(&fptr,fname, &status); 		/* create new file 					*/
    fits_report_error(stderr, status);                  /* print out any error messages 			*/
    fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
    fits_report_error(stderr, status);                  /* print out any error messages 			*/
    
    F[0] = (double *)malloc(sizeof(double)*cfg->xsize*cfg->ysize);
    if ( !F[0] ) exit_with_usage(49);

    for ( y=1 ; y < cfg->ysize ; y++ ) F[y] = F[y-1] + cfg->xsize;

    double Ptot=0;
    for ( y=0 ; y<cfg->ysize ; y++ ){
	for ( x=0 ; x<cfg->xsize ; x++){
	    i = cfg->ysize * x + y;
	    F[y][x] = Fh[i];
	    Ptot += Fh[i];
	}
    }

    int size;
    if ( cfg->xsize >= cfg->ysize ) size = cfg->xsize;
    else			    size = cfg->ysize;
	
    float pxscale = (float)(cfg->fov/cfg->dist*pc/size);
    float lambda  = (float)(cfg->w[wj]*1e6);

    fits_update_key(fptr, TFLOAT, "WAVELENG", &lambda,"[micron] Wavelength of dataset", &status);
    fits_update_key(fptr, TDOUBLE, "FLUX", &Ptot,"[Jy] Total flux detected", &status);
    fits_update_key(fptr, TFLOAT, "PIXELSCA", &pxscale,"[''/px] pixel scale of instrument setup", &status);
    fits_update_key(fptr, TSTRING, "VERSION", ver,"DiskDyn version", &status);
    fits_write_history(fptr,"Created with DiskDyn",&status);
    fits_write_date(fptr,&status);

    fits_write_img(fptr, TDOUBLE, fpixel, cfg->xsize*cfg->ysize, F[0], &status);
    fits_report_error(stderr, status);                  /* print out any error messages */
    fits_close_file(fptr, &status); 			/* close the file */
    fits_report_error(stderr, status);                  /* print out any error messages */

    if ( cfg->verb>=2 ) printf("\tFile %s successfully written\n",fname);

    free(fname);
    free(ver);
    free(F[0]);

    return;
}
