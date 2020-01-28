#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef PI
#define PI 3.141592653589793238
#endif

void mie_single(float s, float *lambda, float *nindex, float *kindex, int nlambda, float *costheta, int numthetas, float *Qabs, float *Qsca,  float *Qpfunc){

  // mie_single.c
  // Uses the indeces of refraction to calculate Qsca, Qabs, and Qpfunc
  // 
  // Input:
  // s: particle size in microns
  // lambda: wavelength vector in microns
  // nindex,kindex: indices of refraction at each wavelength
  // nlambda: length of lambda, nindex, and kindex vectors
  // costheta: cosine of scattering angle
  // numthetas: length of costheta vector
  //
  // Output:
  // Qabs: the absorption efficiencies corresponding to each lambda
  // Qsca: the scattering efficiencies corresponding to each lambda
  // Qpfunc: the scattering phase function corresponding to each lambda as a function of costheta
  //
  // Adapted from mie_single.pro IDL code created by G. Thomas
  //
  // Translated from IDL and updated by Christopher Stark
  // NASA GSFC
  // Last updated 23 Apr 2014 by Christopher Stark

  
  int i;

  long maxx;
  long n, nmax;

  double cm_real, cm_imag;
  double Ir_real, Ir_imag;
  double Y_real, Y_imag;
  double nstop;

  double x;
  double *d_real, *d_imag;
  double Tnp1, Tnm1;
  double Psi0, Psi1, Chi0, Chi1, Psi, Chi, Rnx;
  double Dg, Dqsc, Dqxt;
  double temp, temp_real, temp_imag, temp_rn, temp_in, temp_rd, temp_id, temp_rd2, temp_id2, temp_denom;
  double A_real, A_imag, B_real, B_imag;
  double Anm1_real, Anm1_imag, Bnm1_real, Bnm1_imag;

  //SPF variables
  int ipfunc;
  double *Sm_real, *Sm_imag;
  double *Sp_real, *Sp_imag;
  double *Pi0_real, *Pi0_imag;
  double *Pi1_real, *Pi1_imag;
  double A2, DN, Turbo;
  double S_real, S_imag;
  double T_real, T_imag;
  double Taun_real, Taun_imag;
  double Xs1_real, Xs1_imag;
  double Xs2_real, Xs2_imag;
  double A2A_real, A2A_imag, A2B_real, A2B_imag;
  

  maxx = 7200000;

  //Variable initialization for SPF
  Sm_real = (double *) malloc (numthetas * sizeof(double));
  if (Sm_real==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Sm_real array.\n"); exit (1);}
  Sm_imag = (double *) malloc (numthetas * sizeof(double));
  if (Sm_imag==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Sm_imag array.\n"); exit (1);}
  Sp_real = (double *) malloc (numthetas * sizeof(double));
  if (Sp_real==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Sp_real array.\n"); exit (1);}
  Sp_imag = (double *) malloc (numthetas * sizeof(double));
  if (Sp_imag==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Sp_imag array.\n"); exit (1);}
  Pi0_real = (double *) malloc (numthetas * sizeof(double));
  if (Pi0_real==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Pi0_real array.\n"); exit (1);}
  Pi0_imag = (double *) malloc (numthetas * sizeof(double));
  if (Pi0_imag==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Pi0_imag array.\n"); exit (1);}
  Pi1_real = (double *) malloc (numthetas * sizeof(double));
  if (Pi1_real==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Pi1_real array.\n"); exit (1);}
  Pi1_imag = (double *) malloc (numthetas * sizeof(double));
  if (Pi1_imag==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Pi1_imag array.\n"); exit (1);}


  for(i=0;i<nlambda;i++) { //The big loop over all lambda values

    x = 2.0 * PI * ((double) s) / ((double) lambda[i]);
    if (x > maxx) {fprintf(stderr,"ERROR: Size Parameter Overflow in Mie\n"); exit (1);}

    cm_real = (double) nindex[i];
    cm_imag = (double) -kindex[i];
    Ir_real = cm_real / (cm_real*cm_real + cm_imag*cm_imag);
    Ir_imag = - cm_imag / (cm_real*cm_real + cm_imag*cm_imag);
    Y_real =  x * cm_real;
    Y_imag =  x * cm_imag;

    if (x < 0.02)
      nstop = 2.0;
    else if (x <= 8.0)
      nstop = x + 4.0 * pow(x, 0.3333333333333333) + 2.0;
    else if (x < 4200.0)
      nstop = x + 4.05 * pow(x, 0.3333333333333333) + 2.0;
    else
      nstop = x + 4.0 * pow(x, 0.3333333333333333) + 2.0;
    
    temp = sqrt(Y_real * Y_real + Y_imag * Y_imag);
    temp = (nstop > temp) ? nstop : temp; // temp = greater of nstop or magnitude of Y
    temp += 15.0;
    nmax = (long) floor(temp);

    d_real = (double *) malloc ((nmax+1) * sizeof(double));
    if (d_real==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating d_real array.\n"); exit (1);}
    d_imag = (double *) malloc ((nmax+1) * sizeof(double));
    if (d_imag==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating d_imag array.\n"); exit (1);}
    for(n = 0; n < (nmax+1); n++) {d_real[n]=0.0; d_imag[n]=0.0;} //initialize to zero

    for(n = nmax-1; n > 0; n--) {
      temp_real = ( ((double) (n+1)) * Y_real ) / (Y_real*Y_real + Y_imag*Y_imag);
      temp_imag = (- ((double) (n+1)) * Y_imag) / (Y_real*Y_real + Y_imag*Y_imag);
      d_real[n] = temp_real - ( (temp_real + d_real[n+1]) / ( (temp_real + d_real[n+1])*(temp_real + d_real[n+1]) + (temp_imag + d_imag[n+1])*(temp_imag + d_imag[n+1]) ) );
      d_imag[n] = temp_imag - ( - (temp_imag + d_imag[n+1]) / ( (temp_real + d_real[n+1])*(temp_real + d_real[n+1]) + (temp_imag + d_imag[n+1])*(temp_imag + d_imag[n+1]) ) );
    }


    for(n = 0; n < numthetas; n++) {Sm_real[n]=0.0; Sm_imag[n]=0.0;} //initialize to zero
    for(n = 0; n < numthetas; n++) {Sp_real[n]=0.0; Sp_imag[n]=0.0;} //initialize to zero
    for(n = 0; n < numthetas; n++) {Pi0_real[n]=0.0; Pi0_imag[n]=0.0;} //initialize to zero
    for(n = 0; n < numthetas; n++) {Pi1_real[n]=1.0; Pi1_imag[n]=0.0;} //initialize real to 1, imaginary to zero

    Psi0 = cos(x);
    Psi1 = sin(x);
    Chi0 = -Psi1;
    Chi1 = Psi0;

    Dg = 0.0;
    Dqsc = 0.0;
    Dqxt = 0.0;
    Tnp1 = 1.0;

    for(n=1;n<=nstop;n++){ //loop over all n

      Tnp1 += 2.0;
      Tnm1 = Tnp1 - 2.0;

      DN = (double) n;
      A2 = Tnp1 / (DN*(DN+1.0));
      Turbo = (DN+1.0) / DN;

      Rnx = n / x;
      Psi = Tnm1 * (Psi1 / x) - Psi0;
      Chi = Tnm1 * (Chi1 / x) - Chi0;

      temp_rn = ((d_real[n] * Ir_real - d_imag[n] * Ir_imag) + Rnx) * Psi - Psi1; // real part of numerator
      temp_in = (d_real[n] * Ir_imag + d_imag[n] * Ir_real) * Psi; // imag part of numerator
      temp_rd = ((d_real[n] * Ir_real - d_imag[n] * Ir_imag) + Rnx); // real part of denominator
      temp_id = (d_real[n] * Ir_imag + d_imag[n] * Ir_real); // imag part of denominator
      temp_rd2 = temp_rd * Psi - temp_id * Chi;
      temp_id2 = temp_rd * Chi + temp_id * Psi;
      temp_rd = temp_rd2 - Psi1;
      temp_id = temp_id2 - Chi1; //Now I have the numerator and denominator, so I can do complex division
      temp_denom = 1. / (temp_rd * temp_rd + temp_id * temp_id);
      A_real = (temp_rn * temp_rd + temp_in * temp_id) * temp_denom;
      A_imag = (temp_in * temp_rd - temp_rn * temp_id) * temp_denom;


      temp_rn = ((d_real[n] * cm_real - d_imag[n] * cm_imag) + Rnx) * Psi - Psi1; // real part of numerator
      temp_in = (d_real[n] * cm_imag + d_imag[n] * cm_real) * Psi; // imag part of numerator
      temp_rd = ((d_real[n] * cm_real - d_imag[n] * cm_imag) + Rnx); // real part of denominator
      temp_id = (d_real[n] * cm_imag + d_imag[n] * cm_real); // imag part of denominator
      temp_rd2 = temp_rd * Psi - temp_id * Chi;
      temp_id2 = temp_rd * Chi + temp_id * Psi;
      temp_rd = temp_rd2 - Psi1;
      temp_id = temp_id2 - Chi1; //Now I have the numerator and denominator, so I can do complex division
      temp_denom = 1. / (temp_rd * temp_rd + temp_id * temp_id);
      B_real = (temp_rn * temp_rd + temp_in * temp_id) * temp_denom;
      B_imag = (temp_in * temp_rd - temp_rn * temp_id) * temp_denom;

      Dqxt = Tnp1 * (A_real + B_real) + Dqxt;
      Dqsc = Tnp1 * (A_real * A_real + A_imag * A_imag + B_real * B_real + B_imag * B_imag) + Dqsc;

      if(n > 1) Dg = Dg + ((DN * DN - 1.0) / DN) * ((Anm1_real * A_real + Anm1_imag * A_imag) + (Bnm1_real * B_real + Bnm1_imag * B_imag)) + Tnm1 * (Anm1_real * Bnm1_real + Anm1_imag * Bnm1_imag) / ( DN * DN - DN );

      Anm1_real = A_real;
      Anm1_imag = A_imag;
      Bnm1_real = B_real;
      Bnm1_imag = B_imag;

      for(ipfunc=0;ipfunc<numthetas;ipfunc++) {
	S_real = costheta[ipfunc] * Pi1_real[ipfunc];
	S_imag = costheta[ipfunc] * Pi1_imag[ipfunc];
	T_real = S_real - Pi0_real[ipfunc];
	T_imag = S_imag - Pi0_imag[ipfunc];
	Taun_real = n*T_real - Pi0_real[ipfunc];
	Taun_imag = n*T_imag - Pi0_imag[ipfunc];
	
	A2A_real = A2 * A_real;
	A2A_imag = A2 * A_imag;
	A2B_real = A2 * B_real;
	A2B_imag = A2 * B_imag;
	Sp_real[ipfunc] = (A2A_real + A2B_real) * (Pi1_real[ipfunc] + Taun_real) - (A2A_imag + A2B_imag) * (Pi1_imag[ipfunc] + Taun_imag) + Sp_real[ipfunc];
	Sp_imag[ipfunc] = (A2A_imag + A2B_imag) * (Pi1_real[ipfunc] + Taun_real) + (A2A_real + A2B_real) * (Pi1_imag[ipfunc] + Taun_imag) + Sp_imag[ipfunc];

	Sm_real[ipfunc] = (A2A_real - A2B_real) * (Pi1_real[ipfunc] - Taun_real) - (A2A_imag - A2B_imag) * (Pi1_imag[ipfunc] - Taun_imag) + Sm_real[ipfunc];
	Sm_imag[ipfunc] = (A2A_imag - A2B_imag) * (Pi1_real[ipfunc] - Taun_real) + (A2A_real - A2B_real) * (Pi1_imag[ipfunc] - Taun_imag) + Sm_imag[ipfunc];

	Pi0_real[ipfunc] = Pi1_real[ipfunc];
	Pi0_imag[ipfunc] = Pi1_imag[ipfunc];
	Pi1_real[ipfunc] = S_real + Turbo * T_real;
	Pi1_imag[ipfunc] = S_imag + Turbo * T_imag;
      }
      
      Psi0 = Psi1;
      Psi1 = Psi;
      Chi0 = Chi1;
      Chi1 = Chi;

    } // end of n loop

    free(d_real);
    free(d_imag);

    if (Dg > 0) Dg = 2.0 * Dg / Dqsc;

    for(ipfunc=0;ipfunc<numthetas;ipfunc++) {
      Xs1_real = (Sp_real[ipfunc] + Sm_real[ipfunc]) * 0.5;
      Xs1_imag = (Sp_imag[ipfunc] + Sm_imag[ipfunc]) * 0.5;
      Xs2_real = (Sp_real[ipfunc] - Sm_real[ipfunc]) * 0.5;
      Xs2_imag = (Sp_imag[ipfunc] - Sm_imag[ipfunc]) * 0.5;
      //Note: for a 2D array entry c[ix,iy], the 1D index = ix + iy*xsize
      //Here, for Qpfunc: ix = i, iy = ipfunc, and xsize = nlambda
      Qpfunc[i+nlambda*ipfunc] = (float) (((Xs1_real*Xs1_real) + (Xs1_imag*Xs1_imag) + (Xs2_real*Xs2_real) + (Xs2_imag*Xs2_imag)) / (4*PI*Dqsc));
    }

    Dqsc =  2.0 * Dqsc / (x * x);
    Dqxt =  2.0 * Dqxt / (x * x);

    Qabs[i] = (float) Dqxt - Dqsc;
    Qsca[i] = (float) Dqsc;

  } // end of lambda loop

  free(Sm_real);
  free(Sm_imag);
  free(Sp_real);
  free(Sp_imag);
  free(Pi0_real);
  free(Pi0_imag);
  free(Pi1_real);
  free(Pi1_imag);

}
