#ifndef __DDYN_CUH_INCLUDED
#define __DDYN_CUH_INCLUDED

/****************************************************************************************/
/*                                                                                      */
/* The constants used in the code. The dynamical variables in the code are in AU, Msun, */
/* and yr units, while the collisional units are in SI. This is because the timescales  */
/* of the problem are vastly different as are the mass ranges considered.               */
/*                                                                                      */
/****************************************************************************************/

/****************************************************************************************/
/* 				Physical constants					*/
/****************************************************************************************/

#define MEarth  5.97e24                 /* Mass of Earth in kg                          */
#define MMoon   7.36e22                 /* Mass of Moon in kg                           */
#define AU      1.49598e11              /* AU in m                                      */
#define MSun    1.98892e30              /* Mass of Sun in kg                            */
#define RSun	695700000.0		/* Radius of the Sun in m			*/
#define year    3.1556926e7             /* A year in s                                  */
#define lightAU 63239.6717367           /* Speed of light in AU/yr                      */
#define GravSI	6.67384e-11		/* Gravitational constant in SI			*/
#define cspeed  299792458.0   		/* Speed of Light				*/
#define kboltz  1.3806503E-23		/* Boltzmann constant 				*/
#define hplank  6.626068E-34 		/* Plank constant				*/
#define Vconv	4740.57580894		/* 1 AU/yr in m/s                               */
#define PI      3.14159265358979323846	/* value of pi					*/
#define PI2     6.28318530717958623200	/* value of 2pi					*/
#define Grav	(4.0*PI*PI)		/* Grav constant in AU,MSun,yr			*/
#define pc	3.08567758149137e16	/* Parsec in m					*/
#define Tsub	2000.0			/* Sublimation temperature			*/

#define nR	100			/* Number of Radial divides			*/
#define nsize	2700			/* Size bins in constant memory			*/
#define compute	210			/* Minimum GPU compute capability		*/
#define mgrid	1.1			/* Mass seperation for dust particles		*/
#define Cd	2			/* Drag coefficient factor			*/

/****************************************************************************************/
/*				 Data structures					*/
/****************************************************************************************/

struct coord{				/* Particle data structures			*/
    double	*m,*N;			/* mass*Grav (in Solar mass), Number per tracer	*/
    double	*x,*y,*z,*vx,*vy,*vz;	/* Coodinate in AU, velocity in AU/yr		*/
    double	*k1x,*k2x,*k3x,*k4x;	/* RK4 variables				*/
    double	*k1y,*k2y,*k3y,*k4y;
    double 	*k1z,*k2z,*k3z,*k4z;
    double	*k2rx,*k2ry,*k2rz;
    double	*k3rx,*k3ry,*k3rz;
    double	*k4rx,*k4ry,*k4rz;
    double	*k1vx,*k1vy,*k1vz;
    double	*k2vx,*k2vy,*k2vz;
    double	*k3vx,*k3vy,*k3vz;
    double	*k4vx,*k4vy,*k4vz;
    int		*id;			/* size ID for emission output			*/
};

typedef struct{
    double	vx,vy,vz;		/* Only for storage in config file		*/
    double	x,y,z,m;		/* Cartesian coordinates and			*/
} planet_cart;				/* velocities of planets			*/

struct array{				/* Array for shared memory on GPU		*/
    double	m,x,y,z;
};

typedef struct{
    double	wvl,F;			/* Spectral input data				*/
}specdata;				/* This gets interpolated			*/

typedef struct{				/* Disk parameters				*/
    double	Rin;			/* Inner radius in AU				*/
    double	Rout;			/* Outer radius in AU				*/
    double	dhdr;			/* Total disk height				*/
    double	sigmae;			/* Distribution sigma of eccentricities		*/
    double	amax;			/* Largest body in distribution			*/
    double	mmax;			/* Largest body mass in distribution		*/
    double	amin;			/* Smallest dust size in distribution		*/
    double	mmin;			/* Smallest dust mass in distribution		*/
    double	dmax;			/* Largest dust size in distribution		*/
    double	dmmax;			/* Largest dust mass in distribution		*/
    double	eta;			/* Mass distribution slope (read in in size)	*/
    double	gammas;			/* Surface density profile (0: constant)	*/
    double	C;			/* Mass distribution scaling			*/
    double	Mtot;			/* Total mass in distribution			*/
    long int	Ngrav;			/* Number of gravitating bodies			*/
    long int	Nsmall;			/* Number of small emitting particles		*/
} diskconf;

typedef struct{
    double	a;	 		/* Initial semi-major axis (AU) of blob		*/
    double	ecc;	 		/* eccentricity	of particles in blob		*/
    double	inc;			/* Inclination of particles in blob		*/
    double	w;	 		/* Argument of periapsis of particles in blob	*/
    double	o;	 		/* Initial longitude of ascending node of blob	*/
    double	nu;	 		/* True anomaly of particles in blob		*/
    double	sigma;	 		/* Initial sigma (AU) of blob			*/
    double	amax;	 		/* Initial amax of blob				*/
    double	mmax;	 		/* Initial mmax of blob				*/
    double	amin;	 		/* Initial amin of blob				*/
    double	mmin;	 		/* Initial mmin of blob				*/
    double	dmax;			/* Largest dust size in distribution		*/
    double	dmmax;	 		/* Largest dust mass in distribution		*/
    double	eta;	 		/* Initial eta of blob (read in in size)	*/
    double	C;       		/* Mass distribution scaling			*/
    double	Mtot;			/* Total mass in distribution			*/
    double	vkepp;			/* Velocity percent of Keplerian for products	*/
    long int	Ngrav;			/* Number of gravitating bodies in distribution	*/
    long int	Nsmall;			/* Number of small emitting particles		*/
    int		collcloud;		/* Collision (=0) outcome or Cloud (=1)		*/
} blobconf;

typedef struct{
    int		verb;			/* Verbose? 0=none, 1=some, 2=more, 3=all	*/
    int		gpuid;			/* GPUID					*/
    int		gl;			/* 0=no, 1=yes					*/
    double	t_gl;			/* time for graphics				*/
    double	t_write;		/* time for data output				*/
    double	t_sed;			/* time for SED output				*/
    double	t_fits;			/* time for fits output				*/
    double	t_end;			/* time for end					*/
    int		n_step;			/* n_steps per orbit				*/
    char*	sptype;			/* Spectral-type				*/
    double	teff;			/* Stellar Temperature				*/
    double	rstar;			/* Stellar Radius				*/
    double	Bstar;			/* Stellar magnetic field			*/
    double	Wz;			/* Stellar rotation frequency			*/
    double	mloss;			/* Stellar mass loss rate			*/
    double	vsw;			/* Stellar wind speed				*/
    int		radQ;			/* Inc rad forces? 0=no, 1=yes			*/
    int		SWQ;			/* Inc S-W drag? 0=no, 1=yes			*/
    int		FLQ;			/* Inc Lorentz? 0=no, 1=yes			*/
    planet_cart	*planets;		/* Cartesian coordinates			*/
    int		n_planets;		/* Number of planets				*/
    int		n_disks;		/* Number of disks				*/
    diskconf	*disks;			/* Disk parameters				*/
    int		n_blobs;		/* Number of blobs				*/
    blobconf	*blobs;			/* Blob parameters				*/
    double	rho;			/* Bulk density of particles and bodies		*/
    double	kqe;			/* electric charge coefficient  		*/
    char*	name_stub;		/* Name stuf for data outpu			*/
    int		model;			/* Input switch					*/
    double	dist;			/* Distance of system				*/
    double	fov;			/* Field of view of model images		*/
    double	PA;			/* Position angle of image			*/
    double	inc;			/* Inclination of system			*/
    double	longitude;		/* Argument of periapsis offset in images	*/
    int		xsize;			/* x-size of output fits images			*/
    int		ysize;			/* y-size of output fits images			*/
    char*	specfile;		/* Spectral file name				*/
    specdata	*spec;			/* Spectral data read in			*/
    int		nspec;			/* Number of data lines in spectral data file	*/
    double	*wavs;			/* Wavelengths of output data files		*/
    int		*wj;			/* Index of output wavelengths			*/
    int		nwavs;			/* Number of output images			*/
    double	amin;			/* Global smallest dust particle modeled	*/
    double	amax;			/* Global largest dust particle modeled		*/
    double	mmin;			/* Mass of smallest mass dust particles		*/
    double	mmax;			/* Mass of largest mass dust particles		*/
    int		n_grid;			/* Number of grids for dust particles		*/
    double	f;			/* Mass grid ratio scaling			*/
    char*	comp;			/* Optivcal composition file			*/
    double	*w;			/* Interpolation wavelengths			*/
    double	*s;			/* Sizes of dust particles in interpolation	*/
    double	*q;			/* Charges of dust particles in interpolation	*/
    double	*in;			/* Stellar flux in interpolation		*/
    double	*inexc;			/* Excess fluxes from system in interpolation	*/
    double	*beta;			/* beta values of dust particles - to const mem	*/
    double	*Qabs;			/* Absorption coefficients in interpolation	*/
    double	*Qsca;			/* Total scattering coefficients interpolated	*/
    double	*Qpfunc;		/* Scattering phase function array		*/
    double	betasw;			/* Beta of stellar wind factor			*/
    int		n_wav;			/* Number of wavelengths in interpolation	*/
    int		n_theta;		/* Number of angles for SPF in interpolation	*/
    int		n_dust;			/* Number of dust grids				*/
    int		prev;			/* Previous output percentage			*/
    double	dmin;			/* Dust min size in optical constant file	*/
    double	dmax;			/* Dust max size in optical constant file	*/
    double	Lum;			/* Stellar Luminosity				*/
    size_t	ram;			/* Total RAM of system				*/
    size_t	memgpu;			/* Total GPU global memory used			*/
    size_t	memgpudust;		/* Total GPU global memory used for dust	*/
    size_t	mem;			/* Total memory used in RAM			*/
    long	seed;			/* Original random seed				*/
} configdata;

typedef struct{
    int		gpu;			/* GPU id					*/
    char	*name;			/* Name of GPU					*/
    int		major;			/* GPU computation compabality major value	*/
    int		minor;			/* GPU computation compabality minor value	*/
    int		cores;			/* Total number of GPU cores			*/
    int		sm;			/* Total number of multiprocessors		*/
    long int	totalGlobalMem;		/* Total Global memory of GPU			*/
    long int	totalConstMem;		/* Total Constant memory of GPU			*/
    long int	memPitch;		/* memory pitch					*/
    long int	sharedMemPerBlock;	/* Total available shared memory per block	*/
    int		maxThreadsPerBlock;	/* Maximum number of threads per block		*/
} cudaconf;

typedef struct{
    cudaconf	*prop;			/* Property array of cuda cards			*/
    int		ncuda;			/* Total number of CUDA cards in system		*/
} cudadata;

/****************************************************************************************/
/* 				Data Grid class						*/
/****************************************************************************************/

class Grid{
   public:
    int    	n_all;			/* Total number of particles			*/
    int    	n_dust;			/* Number of dust particles			*/
    int    	n_mass;			/* Number of massive particles			*/
    coord  	host;			/* Pointer to array on host			*/
    coord  	gpu;			/* Pointer to array on gpu			*/
    bool	synced;			/* Boolean to whether array is synced 		*/
    long	seed;			/* Seed for random number generator		*/
    double	time;			/* Time of model				*/
    double	dt;			/* timestep for model				*/
    double	tlastd;			/* Last time data was written out		*/
    double	tlastf;			/* Last time fits images were written out	*/
    double	tlasts;			/* Last time an SED was written out		*/
    double	*Qabs;			/* Pointer for absorption coefficients on dev	*/
    double	*Qsca;			/* Pointer for scattering coefficients on dev 	*/
    double	*Qpfc;			/* Pointer for SPF coefficients on dev		*/
    double	*w;			/* Pointer for wavelength array on dev		*/
    double	*in;			/* Pointer for stellar intensity array on dev	*/
    double	*Th;			/* T array between Rmin and Rmax - precalc host */
    double	*Td;			/* T array between Rmin and Rmax - precalc dev	*/
    double	Rmin;			/* Rmin for precalculated T array		*/
    double	Rmax;			/* Rmax for precalculated T array		*/
    double	*Fh;			/* Pointer for image flux array on host		*/
    double	*Fd;			/* Pointer for image flux array on device	*/
    double	*Fdust;			/* Pointer for dust fluxes on device for SED	*/

    Grid();				/* Constructor for class			*/
    ~Grid();				/* Destructor for class				*/

    /* Particle allocation and initiation functions					*/

    void allochost();			/* Allocate memory for host coord array		*/
    void reallochost(int);		/* Re-allocate memory for host coord array	*/
    void allocgpu();			/* Allocate GPU memory for dev coord array	*/
    void sync_to_gpu();			/* Sync particle data to GPU			*/
    void sync_dust_to_gpu();		/* Sync optical data to GPU			*/
    void sync_from_gpu();		/* Sync particle data from GPU			*/
    void null_gpu();			/* NULL out coord array on GPU			*/
    void null_host();			/* NULL out coord array on host			*/
    void setup();			/* Setup particle array for param file		*/
    void read_in(char *);			/* Setup particle array for model file		*/
    void add_to_host(double,double,\	/* Add random particle to host coord array	*/
	 double,double,double,double,\
	 double,double,int,int);
    void null_flux();			/* NULL out fluxes for image array on GPU	*/
    void null_flux_part();		/* NULL out fluxes for dust particles for SED	*/
    void calc_flux(int);		/* Calculate fluxes for image array		*/
    void calc_flux_sed(int);		/* Calculate fluxes for SED at particular wvl	*/
    void write_image(int,const char *);	/* Write out generated image			*/
    void copy_sed(int);			/* Calculate SED at all wavelengths		*/
    void write_sed(const char *);	/* Write out final calculated SED		*/
    void writegrid(const char *);	/* Write out data file				*/

    /* Particles evolution and other calculation functions				*/

    void evolve(const char *);		/* Evolve model					*/
    void massive_RK4(double);		/* Evolve massive particles dynamically	w/ RK4	*/
    void massive_final(double);		/* Final evolution step in RK4 for massive	*/
    void small_RK4(double);		/* Evolve "massless" particles w/ RK4		*/
    void calculate_dust_temp();		/* Pre-calculate dust temperatures		*/
    double det_step();			/* Determine evolutionary steps			*/
};

/****************************************************************************************/
/*				External variables and structures			*/
/****************************************************************************************/

extern __constant__ 	double		size[];
extern __constant__ 	double		beta[];
extern __constant__ 	double		q[];
extern __device__   	double		doubletmp;
extern			configdata 	cconfig;
extern			cudadata 	ccuda;
extern  		cudaconf	*prop;

/****************************************************************************************/
/*				Global device functions					*/
/****************************************************************************************/

/* Null dust particles on GPU								*/
__global__ void null_coord_driver(int,coord);
__global__ void null_flux_driver(int,double *);
__global__ void null_flux_part_driver(int,double *);

/* Runge-Kutta 4th order step drivers							*/
__global__ void RK4calc_driver(int start,int n_gpu,coord gpu,double dt,\
		double Bstar,double Rstar,double Wz);
__global__ void RK4calc_massive0(int,coord);
__global__ void RK4calc_massive1(int,coord);
__global__ void RK4calc_massive2(int,int,int,coord,int);
__global__ void RK4calc_massive3(int,coord,double);
__global__ void RK4calc_massive4(int,int,int,coord,int);
__global__ void RK4calc_massive5(int,coord,double);
__global__ void RK4calc_massive6(int,int,int,coord,int);
__global__ void RK4calc_massive7(int,coord,double);
__global__ void RK4calc_massive8(int,int,int,coord,int);
__global__ void RK4calc_massive_final(int,coord,double);
__global__ void RK4calc_dust(int,int,coord,double,int,double,\
		double,double,double);

/* Calculate dust fluxes on GPU								*/
__global__ void calc_flux_driver(int,int,coord,double,double,double,\
                double *,double *,double *,double *,double *,double *w,double *,\
                int,int,int,int,double,double,double,int,int,double);
__global__ void calc_flux_sed_driver(int,int,coord,double,double,double,\
                double *,double *,double *,double *,double *,double *w,double *,\
                int,int,int,int,double,double,double);
__global__ void sum_flux_driver(int,double *);

/* Double atomic addition definition							*/
__device__ double atomicAdd_double(double* address, double val);

/****************************************************************************************/
/*	 				All function headers				*/
/****************************************************************************************/

/* Read in parameters file								*/
int	read_in_param(char *param);
/* Read in model file									*/
int	read_in_model(char *model,const char *);
/* Convert orbital elements to cartesian data						*/
void	coe2rv(double,double,double,double,double,double,double,double,double *,double *,\
	double *,double *,double *,double *,double);
/* Convert cartesian data to orbital elements						*/
void	rv2coe(double,double,double,double,double,double,double,double,double *,double *,\
	double *,double *,double *,double *,double);

/* System functions									*/
int	cudaprops();
void	checkcuda();
void	checkmemusage(int,int);
double	ran3(long *seed);
size_t	getMemorySize();

/* Calculation functions								*/
void 	add_specfile();
void 	add_lnkfile();
void 	add_optical_data();
void 	getpro(double,double *,double *);
void 	getbod(long double *,double *,double,double);
void 	genorb_disk(int,double,long *,double *,double *,double *,\
	double *,double *,double *,int);
void 	genorb_blob(int,double,long *,double *,double *,double *,\
	double *,double *,double *,int);
void 	interpolate_dust_to_grid();
void 	interpolate_stellar_to_grid();
void 	calculate_beta();
void	detminmax();
void 	confCdisks();
void 	confCblobs();
double	intens(double,double);

/* Print functions									*/
void	print_header(const char *Version,const char *VDate);
void	print_cuda(void);
void	print_setup();
void	printmin(double,int,int);
void	printstatus(double);
void	print_optics();
void	printout(double);
void	printtimestep(double);
void	print_reset_min(int,int);
void	print_reset_max(int,int);
void	printmemall();
void	print_status(int,int);
void	printsystemparams();
void	printdata(double);
void	printsed(double);
void	printsedcalc();
void	printcalcflux(int);
void	printnewline();
void	printdone();

int	tokenize_spaces(char *,char **,int);
void	remove_quotes(char *);
int	char_is_space(int);

/* Error functions									*/
void	exit_with_usage(int t);
void	Cuda_Error(cudaError_t err,const char *file, int line);

/****************************************************************************************/
/*				 Defined simple functions				*/
/****************************************************************************************/

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define CudaError( err ) (Cuda_Error(err, __FILE__, __LINE__ ))

#endif
