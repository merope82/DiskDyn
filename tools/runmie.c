#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265

void mie_single(float,float *,float *,float *,int,float *,int,float *,float *,float *);

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

void add_lnkfile(float **la,float **na,float **ka,int *nlnk,char *fname){
 FILE        *fr;
 int         i,t,j;
 double	     lr,nr,kr;
 float	     *l,*n,*k;
 char        buff[512],*dat[16];

 fr=fopen(fname,"rb");
 if ( fr==NULL )        exit(1);

 l=*la;
 n=*na;
 k=*ka;

 i=0;
 while ( ! feof(fr) ){
    if (fgets(buff,255,fr)==NULL)
        break;
    if ( buff[0]=='#' )
        continue;
    j=tokenize_spaces(buff,dat,15);
    if ( j<3 )
        continue;
    t=0;
    t+=sscanf(dat[0],"%lf",&lr);
    t+=sscanf(dat[1],"%lf",&nr);
    t+=sscanf(dat[2],"%lf",&kr);
    if ( t<3 )
	continue;
    l = (float *)realloc(l,sizeof(float)*(i+1));
    n = (float *)realloc(n,sizeof(float)*(i+1));
    k = (float *)realloc(k,sizeof(float)*(i+1));
    l[i] = (float)(lr);
    n[i] = (float)(nr);
    k[i] = (float)(kr);
    i++;
 }
 fclose(fr);

 *la = l;
 *na = n;
 *ka = k;

 *nlnk=i;
}

const char *errmsg[] = { "Missing input parameter(s)",
/* 1 */ "Missing number of angles",\
/* 2 */ "Missing number of dust particles",\
/* 3 */ "Missing minimum dust size",\
/* 4 */ "Missing maximum dust size",\
/* 5 */ "Missing LNK file",\
/* 6 */ "Missing output filename",\
/* 7 */ "Missing input parameter(s)",\
	NULL};

void print_help(FILE *fw){

    fprintf(fw,"Usage:\t\n");
    fprintf(fw,"\t\t./runmie -Ntheta <number of angles> -Ndust <number of dust sizes>\n");
    fprintf(fw,"\t\t         -amin <min dust size in micron> -amax <max dust size in micron>\n");
    fprintf(fw,"\t\t         -lnkfile <L,N,K input file> -output <output file name>\n");
}

void error(int t){
    
 fprintf(stderr,"ERROR:\t%s\n\n",errmsg[t]);
 print_help(stderr);
 exit(1);

}

int main(int argc,char *argv[]){
 FILE	*fw;
 float	s,*l,*n,*k,*costheta,*Qabs,*Qsca,*Qpfunc;
 char   *fname,*output;
 int	nlnk; 
 int	numthetas	=	90;
 int	ndust		= 	250;
 float	amin		=	0.001;
 float	amax		=	1e5;

 fname=output=NULL;

 if ( argc==1 )                       error(0);
 for ( int i=1 ; i<argc ; i++ ){
 if ( strcmp(argv[i],"-h")==0 ){
    print_help(stdout);
    return(0);
 }
 else if ( strcmp(argv[i],"-Ntheta")==0 ){
    i++;if ( i==argc )                error(1);
    sscanf(argv[i],"%d",&numthetas);
 }
 else if ( strcmp(argv[i],"-Ndust")==0 ){
    i++;if ( i==argc )                error(2);
    sscanf(argv[i],"%d",&ndust);
 }
 else if ( strcmp(argv[i],"-amin")==0 ){
    i++;if ( i==argc )                error(3);
    sscanf(argv[i],"%f",&amin);
 }
 else if ( strcmp(argv[i],"-amax")==0 ){
    i++;if ( i==argc )                error(4);
    sscanf(argv[i],"%f",&amax);
 }
 else if ( strcmp(argv[i],"-lnkfile")==0 ){
    i++;if ( i==argc )                error(5);
    fname=argv[i];
 }
 else if ( strcmp(argv[i],"-output")==0 ){
    i++;if ( i==argc )                error(6);
    output=argv[i];
 }
 else                                 error(7);
 }

 if ( fname ==  NULL )		      error(5);
 if ( output == NULL )		      error(6);

 l = (float *)malloc(sizeof(float));
 n = (float *)malloc(sizeof(float));
 k = (float *)malloc(sizeof(float));

 add_lnkfile(&l,&n,&k,&nlnk,fname);

 Qabs   = (float *)malloc(sizeof(float)*nlnk);
 Qsca   = (float *)malloc(sizeof(float)*nlnk);

 costheta = (float *)malloc(sizeof(float)*numthetas);
 Qpfunc   = (float *)malloc(sizeof(float)*numthetas*nlnk);
 
 for ( int i=0 ; i<numthetas ; i++ ){ costheta[i] = (double)i*2.0/(numthetas-1)-1.0; }

 fw=fopen(output,"wb");

 fprintf(fw,"#\tn_lnk\tnthetas\tndust\tamin\t\tamax\n");
 fprintf(fw,"#\t%d\t%d\t%d\t%.6e\t%.6e\n",nlnk,numthetas,ndust,amin,amax);
 for ( int z=0 ; z<ndust ; z++ ){
    s=amin*exp(log(amax/amin)*z/((float)ndust-1.0));
    printf("Working on dust size %.6f micron\n",s);
    mie_single(s,l,n,k,nlnk,costheta,numthetas,Qabs,Qsca,Qpfunc);
    fprintf(fw,"#s %.6e\n",s);
    for ( int i=0 ; i<nlnk           ; i++ ){
	if ( Qabs[i] < 0 ) Qabs[i]=0;
	if ( Qsca[i] < 0 ) Qsca[i]=0;
	fprintf(fw,"%.6e %.6e %.6e\n",l[i],Qabs[i],Qsca[i]);
    }
    for ( int i=0 ; i<numthetas*nlnk/10 ; i++ ){
	for ( int j=0; j<10 ; j++ ){
	    if ( Qpfunc[i*10+j] < 0 ) Qpfunc[i*10+j] = 0;
	    fprintf(fw,"%.6e ",Qpfunc[i*10+j]);
	}
	fprintf(fw,"\n");
    }
    fflush(fw);
 }
 
 free(l);
 free(n);
 free(k);
 free(costheta);
 free(Qabs);
 free(Qsca);
 free(Qpfunc);

 return(0);
}