#include <blitz/array.h>
#include<iostream>
#include<stdio.h>
#include<fstream>
#include<cmath>
#include<time.h>
#include<unistd.h>
#include<stdlib.h>
#include<fftw3.h>
#include<vector>
#include<map>
#include<set>
#include<bitset>

#ifdef GRAPHICS
#include <g2.h>  
#include <g2_X11.h>
#include <g2_gd.h>
#include <g2_PS.h>
#endif

//#define EPS

#define BIGNUMBER 999999999.9
#define SMALLNUMBER 0.00000001
#define PRECISIONFC 0.005
#define FORWARD 1
#define BACKWARD -1

#ifdef MAIN
char *option = NULL;
int c;   
char creep_flag, creep_flag2,depinning_flag, gs_pbc_flag, 
gs_obc_flag, eq_pbc_flag, eq_obc_flag, x_flag, out_of_equilibrium_ground_state_flag,
chaos_pbc_flag, gs_pbc_parabola_flag, gs_obc_parabola_flag, depinning_avalanches_flag;
int IHMAX;
#else
extern int IHMAX;
#endif


#define BLACK	1
#define CYAN	5
#define MARRON	10
#define BLUE	15
#define RED		19
#define YELLOW	25

#define USING_BITSET_VECTOR
//#define USING_CHAR_ARRAY

using namespace blitz; 

typedef int u_type; 
typedef bitset<2> int_mcut; 
//typedef bit_vector<3> int_mcut; 
typedef double Real; 

class key_Archive
{	
	public:
	Real E;	
	int i0; 
#ifdef USING_CHAR_ARRAY	
	char *u;
#endif	
#ifdef USING_BITSET_VECTOR	
	vector<int_mcut> u;
#endif	
	int idist, istart; 
};

class cmp_Archive
{
	public:
	bool operator()(const key_Archive &x,const key_Archive &y) const 
	{
		if(x.E==y.E){
			//cout << "degeneration energy" << endl;
			if(x.i0 == y.i0)
			{
				int i=0, N = IHMAX-1;				
#ifdef USING_CHAR_ARRAY	
				while(i<N && x.u[i] == y.u[i]) i++;
#endif					
#ifdef USING_BITSET_VECTOR	
				while(i<N && x.u[i].to_ulong() == y.u[i].to_ulong()) i++;
#endif					

#ifdef USING_CHAR_ARRAY	
				return (i < N)?((x.u)[i] < (y.u)[i]):(0>1);	
#endif					
#ifdef USING_BITSET_VECTOR	
				return (i < N)?((x.u)[i].to_ulong() < (y.u)[i].to_ulong()):(0>1);	
#endif					

				}
			else return (x.i0 < y.i0);
		}
		else return (x.E < y.E);		
	}
};


class Tri
{
	public:
	int j,k;  	// coordinates
	Real E; 	// Energy 
	int n;   	// number of points below
	int best;   // number of points below
	Real w; 	// thermal weight 

	
#ifdef OPTIMIZE1
	int * list_below;
#else
	vector<int> list_below;
	vector<Real> bond_below;
	vector<Real> bond_boltz_below;
#endif
		
/*	Tri();
	~Tri();
*/
};

class Fourier_transform
{
	fftw_complex *Sin, *Sout;
	fftw_plan Sp;

	public:
	void init(int);
	void faire(Array<double,1>);
	void faire(Array<double,1> &,Array<double,1> &);

	Array<Real,1> S, S2, Smax, Smaxtot, Si, Imi, Rei;
	int N, N2;	
	void reset();
	void print(const char *);
};

class avalanche
{
	public:
	int r,ri,ai,a;
	Real ritozi;
	int N,Ntot;
	Array<int,1> histo_a, histo_r;
	Array<Real,1> histo_ritoz;
		
	void acum(int,int);
	void reset(int);
	void init(int);
//	void print(int);
};



class dprop
{
	public:
	Real v,vi;
	int N;
	
	void acum(Real);
	void reset();
	void print();
};

class Fils
{
	public:
	Real E;
	int *Is;
	Fils(int n) {Is=(int *)malloc(sizeof(int)*n);};	
	void copy(Real e, Array<int,1> is)
	{
		E=e;
		for(int i=0;i<is.size();i++) Is[i] = is(i);
	};
};


class sistemavmc
{
	public:

	// general functions
	void init();
	void init(int);
	~sistemavmc();
	
	void find_critical();
	void find_gs(int);
	void find_gs(int,int);
	void find_next_meta();
	Real find_next_meta(Real,Real);
	Real find_Besc(Real, Array<int,1> &);
	Real find_Besc2(Real, Array<int,1> &,Array<int,1> &,Array<int,1> &, int, int&, int&);
	Real Besc(Real);
	Real Besc_new(Real);
	Real Besc_new2(Real);
	int grow_Tree(int, Real);
	int grow_Tree_new(int, Real);
	int grow_Tree_new2(int, Real);
	void adjust_mcut();
	void check(vector<Fils>);
	u_type * create_fils(Array<int,1> &);
	key_Archive create_fils_set(Array<int,1> &, Real,int,int);
	key_Archive create_mere_set(Array<int,1> & Y, Real E);
	map <Real,u_type*> Archive;
	set <key_Archive,cmp_Archive> Archive_set;

	void create_fils_set_new(Array<int,1> &, Real,int,int,key_Archive &);
	
	// polymer variables
	Array<int,1> Is,Is_replic,Is_replic2,Is0,Is_old,
	Is_critical,Is_gs_pbc, Is_biggest_barrier, Is_biggest_barrier_excited,
 	Is_biggest_barrier_nextmeta,Is_limit, Is_metastable_now;
	Array<Real,1> Is_diffusion;
	Array<Real,1> Is_diffusion_inst;
	Array<Real,1> Is_thermal_av;
	Array<Real,1> Is_energy;
	Array<Real,1> Is_energy2;
	Array<Real,1> Is_entropy;
	Array<int,2> Is_density;
	Array<int,1> width_zero_T;
	int Ihmax;
	Array<int,1> Ih_table;
	Fourier_transform Fou;
	avalanche Ropt,Ropt2;
	dprop Ene;	
	int Period;
	vector<Fils> Tree; 	
	int tiempo;
	double shift;
	Array<double,1> tilted;


	// polymer functions
	Real VMC(int);
	Real VMC(int,float,float);
	Real VMC_Besc(int);
	int Minimize(int);	
	int Minimize_GS(int,int);	
	int Minimize_EQ_avec_T(int);	
	int Minimize_EQ_avec_T0(int);	
	int Minimize_local(int, int, Array<int,1> &);
	int find_Is_limit(int);
	int find_Is_limit_plus(int);
	int Minimize_local_constrained(int , int, Array<int,1> &, Array<int,1> &,Array<int,1> &, int);
	int place_for_fils(Real, int);

	Real Energy();
	Real Energy_conf(Array<int,1> );
	Real Energy_elastic(Array<int,1> );
	Real Energy_pinning(Array<int,1> );
	Real Width_conf(Array<int,1> );
	void Density();
	void reset_Density();
	void print_Density();
	void print_conf(ofstream &,Array<int,1>);
	void print_conf(ofstream &,Array<double,1>);
	void print_conf2(Array<int,1>);
	void fftw(int);
	void print_config(Array<int,1>);
	void check_mcut(Array<int,1>);

			
	// disordered medium variables
	Array<Real,2> V;
	Array<Real,2> V_boltz;
	int M,Mnet,Mcut,ilmax, Ntot;
	Array<int,1> IM_table;
	Real Adis, F, Fc, Cel, beta, F_fraction;
	Array<Tri,1> Net;
	unsigned long initial_seed;
	Real parabola_position;
	Real mass_squared;

	// disordered medium functions
	Real potV(int i, int j)
	{
		Real dw=parabola_position-IM(j);
		if(dw>M*0.5) dw-=M;
		if(dw<-M*0.5) dw+=M;
		#ifndef TIPDRIVEN
		return V(Ih(i),IM(j))+mass_squared*dw*dw*0.5;
		#else
		return V(Ih(i),IM(j))+(i==0)?(dw*dw*0.5):(0.0);
		#endif
	}
	Real potV_boltz(int i, int j)
	{
		Real dw=parabola_position-IM(j);
		if(dw>M*0.5) dw-=M;
		if(dw<-M*0.5) dw+=M;

		#ifndef TIPDRIVEN
		return V_boltz(Ih(i),IM(j))*exp(-beta*mass_squared*dw*dw*0.5);
		#else
		return V_boltz(Ih(i),IM(j))*(i==0)?(exp(-beta*dw*dw*0.5)):(1.0);
		#endif	
	}

	int IM(int i)
	{
		#ifdef OPTIMIZE5
		if(i>M*Mcut*6) {cout << "error: IM_table(big), big = " << i << endl; exit(1);}
		if(i<-M*Mcut) {cout << "error: IM_table(-big), big = " << i << endl; exit(1);}
		return IM_table(i);
		#else
		return ((i+100*M)%M);
		#endif
	};
	
	int Ih(int i)	
	{
		#ifdef OPTIMIZE5
		if(i>2*Ihmax) {cout << "error: Ih_table(big)" << endl; exit(1);}
		if(i<-Ihmax) {cout << "error: Ih_table(-big)" << endl; exit(1);}
		return Ih_table(i);
		#else
		return ((i+100*Ihmax)%Ihmax);
		#endif
	}; 

	int INet(int,int);
	void genera_Net();
	void change_disorder();
	void set_parabola(double,double);
	void copy_disorder(sistemavmc &,double);
	void print_avalanches(int);
	void print_diffusion(int);
	void print_diffusion_inst(void);
	void correlate_disorder();

	// movie
	int Ok;
	int movie, movie_eps, movie_png;
	void cinema(int,int, Array<int,1>);	
	void cinema(int,int, Array<double,1>);	
	void cinema_set_png_name();
	void cinema_close_png_name();

	// output files
	ofstream Difi;
	ofstream out_meta;	
	ofstream out_conf_biggest_barrier;	
	ofstream out_conf_biggest_barrier_excited;	
	ofstream out_conf_biggest_barrier_nextmeta;	
	ofstream out_conf_ground_state;	
	ofstream out_conf_critical;	
	ofstream logfile;	
	ofstream out_dis;
};


/*----- generadores de numeros random --------*/

/* initializes state[N] with a seed */
extern void init_genrand(unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
extern void init_by_array(unsigned long *init_key, unsigned long key_length);

/* generates a random number on [0,0xffffffff]-interval */
extern unsigned long genrand_int32(void);

/* generates a random number on [0,0x7fffffff]-interval */
extern long genrand_int31(void);

/* generates a random number on [0,1]-Real-interval */
extern double genrand_real1(void);

/* generates a random number on [0,1)-Real-interval */
extern double genrand_real2(void);

/* generates a random number on (0,1)-Real-interval */
extern double genrand_real3(void);

/* generates a random number on [0,1) with 53-bit resolution*/
extern double genrand_res53(void); 

/* generates gaussian deviation */
extern double gasdev_real3(void);

