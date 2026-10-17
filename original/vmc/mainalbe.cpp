/* 
Ultima modificacion: Jan21 2014.
Ultima modificacion: Sep27 2014.
*/

#include <unistd.h>
#include <string.h>
#define MAIN
#include "vmc.hpp"
#undef MAIN

int out_of_equilibrium_ground_state();
int critical_configuration();
int ground_state_pbc();
int ground_state_pbc_and_critical();
int ground_state_chaos_pbc();
int ground_state_pbc_with_parabola();
int ground_state_obc_with_parabola();
int depinning_avalanches();

sistemavmc P, P1, P2;

class main_vmc{};

int main (int argc, char **argv)
{	
	int command_line_options(int, char **);
	command_line_options(argc, argv);
	
	if(out_of_equilibrium_ground_state_flag==1)
	out_of_equilibrium_ground_state();

	if(depinning_flag==1 && gs_pbc_flag==0)
	critical_configuration();

	if(gs_pbc_flag==1 && depinning_flag==0)
	ground_state_pbc();
	
	if(gs_pbc_flag==1 && depinning_flag==1)
	{
		cout << "calculating ground-state and critcal" << endl;
		ground_state_pbc_and_critical();
	}

	if(chaos_pbc_flag==1 && depinning_flag==0)
	ground_state_chaos_pbc();

	if(gs_pbc_parabola_flag==1 && depinning_flag==0)
	ground_state_pbc_with_parabola();

	if(gs_obc_parabola_flag==1 && depinning_flag==0)
	ground_state_obc_with_parabola();

	if(depinning_avalanches_flag==1)
	depinning_avalanches();

	return 0;
}


int critical_configuration()
{
	cout << "calculating critical configuration\n";
	P.init(1);
	P.tiempo = 0;

	Fourier_transform fft_critical;
	fft_critical.init(P.Ihmax);
	fft_critical.reset();

	Array<double,1> A(P.Ihmax);

	ofstream fc("fc.dat");
	for(int i=0;i<P.Ntot;i++){
		cout << "changing disorder...";
		P.change_disorder();
		//P.correlate_disorder();

		cout << "done\n";

		cout << "searching critical..."; cout.flush();
		P.find_critical();
		cout << "done\n";

		for(int j=0;j<P.Ihmax;j++) A(j) = P.Is_critical(j)-P.tilted(j);
		//P.adjust_mcut();
		P.print_conf(P.out_conf_critical,A);
		fc << P.Fc << " " << mean(pow2(P.Is_critical-mean(P.Is_critical))) << endl;

		P.check_mcut(P.Is_critical);
		//cout << A << endl;
		fft_critical.faire(A);
		fft_critical.print("critical_partial.fft.dat");

		#ifdef GRAPHICS
		P.cinema(0,RED,A);
		#endif
	}
	fft_critical.print("critical.fft.dat");
	//sleep(5);
	return 0;
}

int ground_state_pbc()
{
	cout << "calculating ground state configuration\n";
	P.init();
	ofstream Ene("GS_Energy.dat");
	ofstream maxANDwidth("maxANDwidth.dat");

	Fourier_transform fft_ground;
	fft_ground.init(P.Ihmax);
	fft_ground.reset();

	Array<double,1> A(P.Ihmax);

	//Array<int,1> dx(8); 
	//dx= 2, 4, 8, 16, 32, 64, 128, 256;
	Array<double,1> w2(P.Ihmax); 
	Array<double,1> w2b(P.Ihmax); 

	w2b=0.0;
	w2=0.0;
	int half = (int)(P.Ihmax*0.5)-1;

	for(int i=0;i<P.Ntot;i++){
		P.change_disorder();
		P.find_gs(0);

		for(int j=0;j<P.Ihmax;j++) A(j) = P.Is_gs_pbc(j);

		fft_ground.faire(A);

		P.cinema(2,BLACK,P.Is_gs_pbc); sleep(1);

		P.print_conf(P.out_conf_ground_state,P.Is_gs_pbc);
		//Ene << P.Energy_conf(P.Is_gs_pbc) << endl;
		Ene << P.Energy_conf(P.Is_gs_pbc) << " " << P.Width_conf(P.Is_gs_pbc) << endl; 
		fft_ground.print("ground_partial.fft.dat");

/*
		maxANDwidth << max(P.Is_gs_pbc-mean(P.Is_gs_pbc))
					<< " " << mean(pow2(P.Is_gs_pbc-mean(P.Is_gs_pbc)))
		<< endl;

		for(int n=0;n<P.Ihmax;n++){
			w2(n)+= pow2(A(n)-A(0));
			w2b(n)+= pow2(A((half+n)%P.Ihmax)-A(half));
		}

		if(i%100){
			ofstream albe("alberto_partial.dat");
			for(int n=0;n<P.Ihmax;n++)
			{
				albe << n << " " << w2(n)/i << " " << w2b(n)/i << endl;
			}
		}
*/
	}

	fft_ground.print("ground.fft.dat");

/*
	ofstream albe("alberto.dat");
	for(int n=0;n<P.Ihmax;n++)
	{
		albe << n << " " << w2(n)/P.Ntot << " " << w2b(n)/P.Ntot << endl;
	}
*/
	return 0;
}

/////////////////////////////////////////////////////////////////
int ground_state_pbc_and_critical()
{
	P.init();
	P.tiempo = 0;	

	ofstream Ene("GS_Energy.dat");
	ofstream fc("fc.dat");

	cout << "ground state and critical chosen" << endl;

	Fourier_transform fft_critical;
	Fourier_transform fft_ground;
	
	fft_critical.init(P.Ihmax);
	fft_ground.init(P.Ihmax);

	fft_critical.reset();
	fft_ground.reset();

	for(int i=0;i<P.Ntot;i++){
		P.change_disorder();

		P.find_gs(0);
		Array<double,1> A(P.Is_gs_pbc-P.tilted);
		P.print_conf(P.out_conf_ground_state,A);
		Ene << P.Energy_conf(P.Is_gs_pbc) << " " << P.Width_conf(P.Is_gs_pbc) << endl; 

		Ene.flush();
		fft_ground.faire(A);
		P.check_mcut(P.Is_gs_pbc);

		P.find_critical();
		Array<double,1> B(P.Is_critical-P.tilted);
		P.print_conf(P.out_conf_critical,B);
		fc << P.Fc << " " << P.Width_conf(P.Is_critical) << endl;
		fc.flush();
		fft_critical.faire(B);
		P.check_mcut(P.Is_critical);

		#ifdef GRAPHICS
		P.cinema(2,RED,A);           
		P.cinema(2,BLUE,B);           
		#endif
		fft_critical.print("critical_partial.fft.dat");
		fft_ground.print("ground_partial.fft.dat");


	}
	fft_critical.print("critical.fft.dat");
	fft_ground.print("ground.fft.dat");

	return 0;
}
/////////////////////////////////////////////////////////////////

int out_of_equilibrium_ground_state()
{
	P.init();
	
	ofstream fc("fc.dat");
	ofstream Ene("GS_Energy.dat");

	ofstream ub("ub.dat");
	ofstream ubi("ubi.dat");
	ofstream ubi2("ubi2.dat");

	ofstream statesout("states.dat");
	ofstream eventsout("events.dat");
	eventsout << "Leverel Severel zlowrel zhighrel barrier lopt\n"; 

	vector<Real> ub_vector;
	vector<Real> lopt_vector;
	vector<Real> lrel_vector;
	vector<Real> cm_vector;

	Fourier_transform fft_critical;
	Fourier_transform fft_ground;
	Fourier_transform fft_creep;
	Fourier_transform fft_meta;
	
	fft_critical.init(P.Ihmax);
	fft_ground.init(P.Ihmax);
	fft_creep.init(P.Ihmax);
	fft_meta.init(P.Ihmax);

	fft_critical.reset();
	fft_ground.reset();
	fft_creep.reset();
	fft_meta.reset();
	
	Array<double,1> A(P.Ihmax);

	for(int i=0;i<P.Ntot;i++){
	try{
		P.change_disorder();
		//P.correlate_disorder();


		cout << "finding critical state" << endl; 
		P.find_critical();
		fc << P.Fc << endl;
		cout << "fc = " << P.Fc << endl;
		#ifdef ABSOLUTEDRIVINGFORCE
		P.Fc=1.0;
		#endif

		#ifdef GRAPHICS 
		//P.cinema(0,19,P.Is_critical);          
		#endif

		cout << "finding ground state (PBC)..." << endl;
		P.find_gs(0);
		Ene << P.Energy_conf(P.Is_gs_pbc) << endl;

		//P.Is = P.Is_gs_pbc; //sets initial condition
		P.Is = 0;

		P.F = P.F_fraction*P.Fc;
		ofstream fuerza("force.dat");
		fuerza.precision(6);fuerza.width(12);
		fuerza.flags(ios_base::fixed);
		fuerza << P.F;
		#ifdef GRAPHICS 
		//P.cinema(2,19,P.Is_gs_pbc);          
		#endif

		cout << "finding next metastable state...";
		P.find_next_meta();	
		cout << "OK" << endl;
		#ifdef GRAPHICS 
		P.cinema(2,BLACK,P.Is);          
		#endif

		Array<int,1> Is_meta, Is_excited, Is_meta_next, Is_first;
		Is_meta.resize(P.Ihmax);
		Is_excited.resize(P.Ihmax);
		Is_meta_next.resize(P.Ihmax);
		Is_first.resize(P.Ihmax);
		
		Is_first = P.Is;		

#ifdef EXACT
		for(int smax=4; smax<5; smax++)
#else
		for(int smax=0; smax<8; smax++)		
#endif			
		{
			Real barrier=0, best_barrier=-1; 
			int best_size, best_size_big, size_ava, size_ava_big;

			int flag_next = 0;
			P.Is = Is_first;
			#ifdef GRAPHICS 
			//P.cinema(2,RED,P.Is_replic2);          
			P.cinema(0,BLACK,P.Is);          
			P.cinema(1,BLACK,P.Is);          
			//sleep(1);
			#endif

			P.print_conf2(P.Is);
			do
			{
				Is_meta = P.Is;
				
				/*--- find Besc, and move to the next metastable state ---*/
				cout << "searching Besc...";
				
				barrier = P.find_Besc2(best_barrier,Is_meta, Is_excited, Is_meta_next, 
									smax, size_ava, size_ava_big);
				cout << "OK Besc" << endl;
				cout.flush();

				/* --  check periodicity ---*/
				int isperiodic=-1;   
				for(int i=0;i<ub_vector.size();i++)
				{				
					if(mean(Is_meta)==cm_vector[i]+P.M)
					{ 
					    isperiodic=i;
					}  
				}
				cout.flush();
	
				A = Is_meta;		
				fft_meta.faire(A);
				fft_meta.print("meta_partial.fft.dat.dat");	
			
				// if periodicity is not reached yet
				if(isperiodic<0){
				    // records current's metastable state data
				    ub_vector.push_back(barrier);
				    lopt_vector.push_back(size_ava);
				    lrel_vector.push_back(size_ava_big);
				    cm_vector.push_back(mean(Is_meta));

				    ubi2.precision(6);ubi2.width(12);ubi2.flags(ios_base::fixed);
				    ubi2 << barrier << " " << size_ava << " " << size_ava_big << " " 
				    << mean(Is_meta) << " " << sum(Is_meta_next-Is_meta) << "\n";

				    int zlow, zhigh;
				    for(int i=0; i < P.Ihmax ; i++) {
						//ubi2 << "\t" << P.IM(Is_meta(i));
				    	//ubi2 << endl;
						//statesout << P.IM(Is_meta(i)) << " " << P.IM(Is_meta_next(i)) << "\n";

						///////
						int ip1=(i+1)%P.Ihmax;
						if(
						P.IM(Is_meta(ip1))!=P.IM(Is_meta_next(ip1)) && 
						P.IM(Is_meta(i))==P.IM(Is_meta_next(i))
						)
						zlow=i;//zlow

						int im1=(i-1)%P.Ihmax;
						if(
						P.IM(Is_meta(im1))!=P.IM(Is_meta_next(im1)) && 
						P.IM(Is_meta(i))==P.IM(Is_meta_next(i))
						)
						zhigh=i;//zlow

						/*if(
						P.IM(Is_meta(ip1))==P.IM(Is_meta_next(ip1)) && 
						P.IM(Is_meta(i))!=P.IM(Is_meta_next(i))
						)
						eventsout << (i+1)%P.Ihmax << " ";//zhigh*/
						///////						
				    }
				    eventsout 
				    << size_ava_big << " " <<  sum(Is_meta_next-Is_meta) 
				    << " " << zlow-(zlow>zhigh)*P.Ihmax << " " << zhigh 
				    << " " << barrier << " " << size_ava << std::endl;

	  			    //statesout << "\n\n";	

				    #ifdef GRAPHICS 
				    #ifndef MOVIE_PNG
				    //P.cinema(2,RED,P.Is_replic2);          
				    P.cinema(2,CYAN,Is_meta);          
				    P.cinema(2,RED,P.Is);          
				    //sleep(1);
				    #else 
				    P.cinema_set_png_name();
				    P.cinema(0,BLACK,Is_meta);          
				    P.cinema(1,BLACK,Is_meta);          
				    P.cinema_close_png_name();
  
				    P.cinema_set_png_name();
				    P.cinema(0,BLACK,Is_meta);          
				    P.cinema(1,BLACK,Is_meta);          
				    P.cinema(2,RED,Is_excited);          
				    //P.cinema(2,BLACK,Is_meta_next);          
				    P.cinema_close_png_name();
  
				    P.cinema_set_png_name();
				    P.cinema(0,BLACK,Is_meta);          
				    P.cinema(1,BLACK,Is_meta);          
				    P.cinema(2,RED,Is_excited);          
				    P.cinema(2,BLACK,Is_meta_next);          
				    P.cinema_close_png_name();
				    #endif
				    #endif
  
				    if(barrier > best_barrier) 
				    {
					  best_barrier = barrier;
					  P.Is_biggest_barrier = Is_meta;
					  P.Is_biggest_barrier_excited = Is_excited;
					  P.Is_biggest_barrier_nextmeta = Is_meta_next;
					  best_size = size_ava;
					  best_size_big = size_ava_big;
				    }
				}
				else if(mean(Is_meta-Is_first)>P.M*0.5)
				{
					Real best_roughness;
					cout << "Biggest barrier dominant = " << best_barrier << endl;
					cout << "Energy of dominant = " << P.Energy_conf(P.Is_biggest_barrier) << endl;
					cout << "Roughness of dominant = " << (best_roughness=P.Width_conf(P.Is_biggest_barrier)) << endl;
					
					ub.precision(6);ub.width(12);ub.flags(ios_base::fixed);
					ub << best_barrier << " " << best_size << " " << best_size_big << " " << best_roughness << endl;

					for(int i=isperiodic;i<ub_vector.size();i++){
					    ubi << ub_vector[i]-best_barrier << " " << lopt_vector[i] << " " 
						<< lrel_vector[i] << " " << cm_vector[i] << endl;
					}

					A = P.Is_biggest_barrier-P.tilted;		
					P.print_conf(P.out_conf_biggest_barrier,A);
					fft_creep.faire(A);					

					P.print_conf(P.out_conf_biggest_barrier_excited,P.Is_biggest_barrier_excited);
					P.print_conf(P.out_conf_biggest_barrier_nextmeta,P.Is_biggest_barrier_nextmeta);

					A = P.Is_critical-P.tilted;		
					P.print_conf(P.out_conf_critical,A);
					fft_critical.faire(A);					

					A = P.Is_gs_pbc-P.tilted;		
					P.print_conf(P.out_conf_ground_state,A);
					fft_ground.faire(A);					
							
					#ifdef GRAPHICS 
					P.cinema(2,YELLOW,P.Is_biggest_barrier);          			
					P.cinema(2,YELLOW,P.Is_biggest_barrier_excited);          			
					P.cinema(2,YELLOW,P.Is_biggest_barrier_nextmeta);          			
					//sleep(3);

					cout << "configuration with biggest barrier (yellow)" << endl;
	
					P.cinema(2,YELLOW,P.Is_biggest_barrier); 		         		
					P.cinema(2,YELLOW,P.Is_biggest_barrier_excited);          			
					P.cinema(2,YELLOW,P.Is_biggest_barrier_nextmeta);          			
					//sleep(1);

					cout << "critical configuration (red)" << endl;
					P.cinema(2,RED,P.Is_critical);          
					//sleep(1);

					cout << "(f=0) ground state (black)" << endl;
					P.cinema(2,BLACK,P.Is_gs_pbc);          
					//system(" xmms ~/mnt/xp/Documents_Ale/jazz/Libertango/Libertango.mp3");
					//while(getc(stdin)!='c');
					//exit(1);
					//sleep(5);
					#endif
					
					flag_next = 1;
				} 
			} while(flag_next==0);

			ubi << "######" << endl;
			ubi2 << "######" << endl;
			ub_vector.clear();
			lopt_vector.clear();
			lrel_vector.clear();
			cm_vector.clear();
			eventsout << "\n\n";
		} 

	
		fft_critical.print("critical_partial.fft.dat");
		fft_ground.print("ground_partial.fft.dat");
		fft_creep.print("creep_partial.fft.dat");
		fft_meta.print("meta_partial.fft.dat");
	}
	catch (bad_alloc&)
	{
		cout << "Error allocating memory." << endl;
		cout << P.Archive_set.size() << " elements." << endl;
		cout.flush();
		exit(1);
	}
	}		

	fft_critical.print("critical.fft.dat");
	fft_ground.print("ground.fft.dat");
	fft_creep.print("creep.fft.dat");
	
	return 1;
}



/*------ Command line options ----------*/
int command_line_options(int argc, char **argv)
{
	creep_flag=x_flag=creep_flag2=depinning_flag=depinning_avalanches_flag=gs_pbc_flag
		    	 =gs_obc_flag=eq_pbc_flag=eq_obc_flag=out_of_equilibrium_ground_state_flag=0;
	/*---- COMMAND LINE OPTIONS ----*/
	while ((c = getopt (argc, argv, "hdefgicxCojkDl")) != -1)
	switch (c)
	{
    	   case 'c':
            	 option = optarg;
		 creep_flag = 1;
            	 break;
    	   case 'C':
            	 option = optarg;
		 creep_flag = 1;
		 creep_flag2 = 1;
            	 break;
    	   case 'd':
            	 depinning_flag = 1;
            	 break;//depinning_avalanches_flag
    	   case 'D':
            	 depinning_avalanches_flag = 1;
            	 break;//depinning_avalanches_flag
    	   case 'e':
            	 gs_pbc_flag = 1;
            	 break;
    	   case 'f':
            	 gs_obc_flag = 1;
            	 break;
    	   case 'g':
            	 eq_pbc_flag = 1;
            	 break;
    	   case 'o':
            	 out_of_equilibrium_ground_state_flag = 1;
            	 break;
    	   case 'i':
            	 eq_obc_flag = 1;
            	 break;		   
    	   case 'j':
            	 chaos_pbc_flag = 1;
            	 break;
		   case 'x':
            	 x_flag = 1;
            	 break;
    	   case 'k':
            	 gs_pbc_parabola_flag = 1;
            	 break;
    	   case 'l':
            	 gs_obc_parabola_flag = 1;
            	 break;
   	   case 'h':
		cout << "choose an option: -c creep, -e,f ground state pbc/obc, "
				"-g,i equilibrium pbc/obc, -d depinning, -o out of eq GS, -k gs pbc parabola" << endl;
		cout << "-D depinning avalanches" << endl;
		return 1;
    	   case '?':
            	 if (isprint (optopt)) fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            	 else fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
            	 return 1;
    	   default:
 		 abort();
	}
	
	if(creep_flag+x_flag+depinning_flag+gs_pbc_flag+gs_obc_flag+eq_pbc_flag
			+eq_obc_flag+out_of_equilibrium_ground_state_flag+chaos_pbc_flag+gs_pbc_parabola_flag
			+depinning_avalanches_flag+gs_obc_parabola_flag==0)
	{
		cout << "choose an option: \n -c creep, \n -e,f ground state pbc/obc, "
				"\n -g,i equilibrium pbc/obc," <<
		" -d depinning, \n -o out of eq GS, \n -j chaos pbc, -k/l ground state pbc/obc parabola\n" 
		" -D depinning avalanches" 
		<< endl;
		exit(1);
	}			
	return 1;	
}



/////////////////////////////////////////////////////////////////////
int ground_state_chaos_pbc()
{
	P1.init();
	P2.init();

	Array<double,1> A1(P1.Ihmax);
	Array<double,1> A2(P2.Ihmax);
	Array<double,1> B12(P1.Ihmax);
	Array<double,1> B11(P1.Ihmax);
	Array<double,1> B22(P1.Ihmax);

	ofstream Ene1("GS_Energy1.dat");
	ofstream Ene2("GS_Energy2.dat");
	ofstream Conf1("GS_CONF1.dat");
	ofstream Conf2("GS_CONF2.dat");

	Fourier_transform fft_ground;
	fft_ground.init(P1.Ihmax);
	fft_ground.reset();

	Fourier_transform fft_ground2;
	fft_ground2.init(P2.Ihmax);
	fft_ground2.reset();

	int half = (int)(P1.Ihmax*0.5)-1;
	B12=0.0;
	B11=0.0;
	B22=0.0;

	double delta;
	cin >> delta;
	cout << delta << endl;

	for(int i=0;i<P1.Ntot;i++){
		P1.change_disorder();
		P2.copy_disorder(P1, delta);
		//for(int i0;i0<P1.M;i0++)
		{
			P1.find_gs(0);
			P2.find_gs(0);
			for(int j=0;j<P1.Ihmax;j++) A1(j) = P1.Is_gs_pbc(j);
			for(int j=0;j<P2.Ihmax;j++) A2(j) = P2.Is_gs_pbc(j);

			/*
			for(int j=0;j<P1.Ihmax;j++) A1(j) = P1.width_zero_T(j);
			for(int j=0;j<P2.Ihmax;j++) A2(j) = P2.width_zero_T(j);
			A1=A1-A1(0);
			A2=A2-A2(0);
			*/

			// critical configuration
			/*P1.find_critical();
			P2.find_critical();
			for(int j=0;j<P1.Ihmax;j++) A1(j) = P1.Is_critical(j);
			for(int j=0;j<P2.Ihmax;j++) A2(j) = P2.Is_critical(j);
			*/

			// fft overlap
			fft_ground.faire(A1,A2);
			fft_ground2.faire(A1);
			fft_ground.print("ground_partial.fft.dat");
			fft_ground2.print("ground_partial2.fft.dat");
			system("paste ground_partial.fft.dat ground_partial2.fft.dat > zzz.fft.dat");


			/*B12 = B12 + A1*A2;
			B11 = B11 + A1*A1;
			B22 = B22 + A2*A2;
			ofstream Bout("B12.dat");
			for(int j=0;j<P1.Ihmax;j++) {
				Bout << B12(j)/(i+1) << " " << B11(j)/(i+1)<< " " << B22(j)/(i+1) << " " << endl;
			}
			Bout << endl << endl;
			*/

			//P1.print_conf(P1.out_conf_critical,A1);
			//P1.print_conf(P1.out_conf_critical,A2);
			//Ene1 << P1.Energy_conf(P1.Is_gs_pbc) << endl;
			//Ene2 << P2.Energy_conf(P2.Is_gs_pbc) << endl;
			//P.cinema(2,BLACK,P.Is_gs_pbc); sleep(1);
			//P1.print_conf(Conf1,P1.Is_gs_pbc);
			//P2.print_conf(Conf2,P2.Is_gs_pbc);
		}
}

	fft_ground.print("ground.fft.dat");
	fft_ground2.print("ground2.fft.dat");

	return 0;
}

/////////////////////////////////////////////////////////////////////
#define FACTORMASS 3.0

double get_mass2(int L){
	double mass2=(10.0*FACTORMASS/L)*(10.0*FACTORMASS/L);
	return mass2;
}

int depinning_avalanches()
{
	P.init();
	ofstream lopt("lava.dat");
	ofstream cluster("clusters.dat");
	Fourier_transform fft_meta;
	fft_meta.init(P.Ihmax);
	fft_meta.reset();

	// lateral size avalanche cut off equal to system size 
	Real mass2=get_mass2(P.Ihmax);
	// choose otherwise (ex 40 points):
	std::cout << "mass2=" << mass2 << std::endl;

	// parabola's step
	//double dw=double(P.M)/double(P.Ntot);
	Real dw=((0.5/P.Ihmax)/mass2)/5.0;
	Real w=10.0;
	
	Array<double,1> A(P.Ihmax);

	P.change_disorder();

	// moving the parabola
	P.Is=0;
	P.F=0; 
	Real tolerance=0.0000001;
	Real delta_min_positive=P.find_next_meta(w, mass2);

	int clusterstart, clusterend, tcluster,lcluster;
	clusterstart=-10000;
	clusterend=-10000;
	tcluster=lcluster=0;

	bool equilibrated=0;

    int inc=0;
	for(int i=0;i<P.Ntot;i+=inc){
	    P.Is0=P.Is;
		delta_min_positive=P.find_next_meta(w, mass2); 
		assert(delta_min_positive>0);
	
		#ifndef TIPDRIVENDEP
		dw = (delta_min_positive+tolerance)/mass2;
		#else
		dw=0.1; //P.Is(0)+=1;
		#endif

		w+=dw;

		cout << "parabola at " << w << ", dw=" << dw << ", delta_min_positive="<< delta_min_positive << endl;

		int nava=0;
		int avastart,avaend;
		for(int j=0;j<P.Ihmax;j++){
			int jp1=(j+1)%P.Ihmax; 
			// start of avalanche detection 000>0<1111 (3), 11100>0< (5), 1110>0<1 (4)
			if( ((P.Is(j)-P.Is0(j))==0) && ((P.Is(jp1)-P.Is0(jp1))!=0) ){
				nava++;
				avastart=j; // border is fixed
			}
			// end of avalanche detection >0<0001111 (0), 111>0<00 (4), 111>0<01 (3)
			if( ((P.Is(j)-P.Is0(j))!=0) && ((P.Is(jp1)-P.Is0(jp1))==0) )
			avaend=jp1; // border is fixed
		}
		//assert(nava>1);
		avastart=((avastart>avaend)?(avastart-P.Ihmax):(avastart));

		int area=sum(P.Is - P.Is0);
		int lenght=count(P.Is > P.Is0)+2 ;

	    if(area>0){
			lopt << (lenght) << " " << (area) << " " << (lenght) << " " << (area) << " "
			<< " " << ((avastart>avaend)?(avastart):(avastart-P.Ihmax)) << " " << avaend << " " 
			<< ((avastart>avaend)?(avastart):(avastart-P.Ihmax)) << " " << avaend << " " << nava ;
			lopt << endl;

			//EZE format
			//fin >> lopt >> sopt >> lrelax >> srelax >> activate_zlow >> activate_zhigh >> relax_zlow >> relax_zhigh;	

			if(!equilibrated) equilibrated=(min(P.Is)>1);
			else
			{	
			A = P.Is;		
			fft_meta.faire(A);		
			fft_meta.print("meta_partial.fft.dat");
			}
		}

		// if event overlaps cluster and cluster not larger than system size
		if(nava==1){
		if((avastart<clusterend && avastart>clusterstart) || (avaend<clusterend && avaend>clusterstart)
		&& lcluster<P.Ihmax)
		{
			// does the cluster grow?
			if(avastart<clusterstart) clusterstart=avastart;	
			if(avaend>clusterend) clusterend=avaend;	

			// lifetime of the cluster
			tcluster++;

			// size of the cluster
			lcluster=clusterend-clusterstart+1;			
		}
		else{
			clusterstart=avastart;
			clusterend=avaend;
			tcluster=0;
			lcluster=avaend-avastart;			
		}
		}
		
		cluster << lcluster << " " << tcluster << endl;




		//for(int j=0;j<P.Ihmax;j++) A(j) = P.Is(j); fft_ground.faire(A);

		// print configuration only in case of a change...
		if(area>0){
			P.cinema(2,BLACK,P.Is); //sleep(1);
			P.print_conf(P.out_conf_critical,P.Is);
			//fft_ground.print("avalanche_partial.fft.dat");
		}
		inc=(area>0);
	}
	fft_meta.print("meta.fft.dat");

}

/////////////////////////////////////////////////////////////////////
int ground_state_with_parabola(int bc)
{
	P.init();
	ofstream Ene("GS_Energy.dat");
	ofstream lopt("lava.dat");
	ofstream cluster("clusters.dat");

	Fourier_transform fft_ground;
	fft_ground.init(P.Ihmax);
	fft_ground.reset();

	Array<double,1> A(P.Ihmax);

	//Array<int,1> dx(8); 
	//dx= 2, 4, 8, 16, 32, 64, 128, 256;
	Array<double,1> w2(P.Ihmax); 
	Array<double,1> w2b(P.Ihmax); 

	w2b=0.0;
	w2=0.0;
	int half = (int)(P.Ihmax*0.5)-1;

	// lateral size avalanche cut off equal to system size 
	double mass2=get_mass2(P.Ihmax);
	std::cout << "mass2=" << mass2 << std::endl;

	// parabola's step
	//double dw=double(P.M)/double(P.Ntot);
	double dw=((0.5/P.Ihmax)/mass2);

	Range R1(0,P.M-2);
	Range R2(1,P.M-1);
	
	// variables for cluster growth
	int clusterstart, clusterend, tcluster,lcluster;
	clusterstart=-10000;
	clusterend=-10000;
	tcluster=lcluster=0;

	P.change_disorder();

	// Ntot is the number of avalanches here
	for(int i=0;i<P.Ntot;i++){

		// moves the parabola
		P.set_parabola(mass2,i*dw);

		P.find_gs(bc,i*dw);

/*
		int nava=0;
		for(int j=0;j<P.Ihmax;j++){
			int jp1=(j+1)%P.Ihmax; 
			if( ((P.Is_gs_pbc(j)-P.Is(j))==0) && ((P.Is_gs_pbc(jp1)-P.Is(jp1))!=0) ){
				nava++;
			}
		}

/////////////////////////
		vector<int> ava_start;
		vector<int> ava_end;
		int nava=0;
		for(int j=0;j<P.Ihmax;j++){
			int jp1=(j+1)%P.Ihmax; 
			int jm1=(j-1+P.Ihmax)%P.Ihmax; 
			if( ((P.Is_gs_pbc(j)-P.Is(j))==0) && ((P.Is_gs_pbc(jp1)-P.Is(jp1))!=0) ){
				ava_start.push_back(j);			
				nava++;
			}
			if( ((P.Is_gs_pbc(j)-P.Is(j))!=0) && ((P.Is_gs_pbc(j)-P.Is(j))==0) ){
				ava_end.push_back(j);			
			}
		} 
		if(ava_start.size()!=nava) {
			std::cout << "ERROR\n";
			exit(1);
		}

		vector<int> sopts;
		vector<int> lopts;
		int shift=(ava_start[0]>ava_end[0]);
		int avas;
		for(avas=0;avas<nava;avas++)
		{
			lopts[avas]=(ava_end[(avas+shift)%nava]-ava_start[avas]+P.Ihmax)%P.Ihmax;

			if(ava_end[(avas+shift)%nava]-ava_start[avas]>0)
			{
				Range Rava(ava_start[avas],ava_end[(avas+shift)%nava]);
				sopts[avas]=sum(P.Is_gs_pbc(Rava) - P.Is(Rava));
			}
			else
			{
				Range Rava1(ava_start[avas],P.Ihmax-1);
				sopts[avas]=sum(P.Is_gs_pbc(Rava1) - P.Is(Rava1));
	
				Range Rava2(0,ava_end[(avas+shift)%nava]);
				sopts[avas]+=sum(P.Is_gs_pbc(Rava2) - P.Is(Rava2));
			}		
		}

/////////////////////////*/

		int nava=0;
		int avastart,avaend;
		for(int j=0;j<P.Ihmax;j++){
			int jp1=(j+1)%P.Ihmax; 
			// start of avalanche detection 000>0<1111 (3), 11100>0< (5), 1110>0<1 (4)
			if( ((P.Is_gs_pbc(j)-P.Is(j))==0) && ((P.Is_gs_pbc(jp1)-P.Is(jp1))!=0) ){
				nava++;
				avastart=j;
			}
			// end of avalanche detection >0<0001111 (0), 111>0<00 (4), 111>0<01 (3)
			if( ((P.Is_gs_pbc(j)-P.Is(j))!=0) && ((P.Is_gs_pbc(jp1)-P.Is(jp1))==0) )
			avaend=jp1;
		}
		//assert(nava>1);
		avastart=((avastart>avaend)?(avastart-P.Ihmax):(avastart));

		int area, lenght;
		lopt << (lenght=count(P.Is_gs_pbc > P.Is)+2) << " " 
		<< (area=sum(P.Is_gs_pbc - P.Is)) << " "
		<< nava << " [ " << avastart << " , " << avaend << " ]";
		lopt << endl;

		// if event overlaps cluster and cluster not larger than system size
		if(nava==1){
		if((avastart<clusterend && avastart>clusterstart) || (avaend<clusterend && avaend>clusterstart)
		&& lcluster<P.Ihmax)
		{
			// does the cluster grow?
			if(avastart<clusterstart) clusterstart=avastart;	
			if(avaend>clusterend) clusterend=avaend;	

			// lifetime of the cluster
			tcluster++;

			// size of the cluster
			lcluster=clusterend-clusterstart+1;			
		}
		else{
			clusterstart=avastart;
			clusterend=avaend;
			tcluster=0;
			lcluster=avaend-avastart;			
		}
		}
		
		cluster << lcluster << " " << tcluster << endl;

		
		// unlike VMC, minimize does not modify the initial state, so we have to force it
		P.Is=P.Is_gs_pbc;

		//for(int j=0;j<P.Ihmax;j++) A(j) = P.Is_gs_pbc(j); fft_ground.faire(A);

		// print configuration only in case of a change...
		if(area>0){
			P.cinema(2,BLACK,P.Is_gs_pbc); //sleep(1);
			P.print_conf(P.out_conf_ground_state,P.Is_gs_pbc);
			Ene << P.Energy_conf(P.Is_gs_pbc) << endl;
			//fft_ground.print("ground_partial.fft.dat");
		}

/*		for(int n=0;n<P.Ihmax;n++){
			w2(n)+= pow2(A(n)-A(0));
			w2b(n)+= pow2(A((half+n)%P.Ihmax)-A(half));
		}

		if(i%100){
			ofstream albe("alberto_partial.dat");
			for(int n=0;n<P.Ihmax;n++)
			{
				albe << n << " " << w2(n)/i << " " << w2b(n)/i << endl;
			}
		}
*/
	}

	//fft_ground.print("ground.fft.dat");

	/*
	ofstream albe("alberto.dat");
	for(int n=0;n<P.Ihmax;n++)
	{
		albe << n << " " << w2(n)/P.Ntot << " " << w2b(n)/P.Ntot << endl;
	}*/

	return 0;
}
int ground_state_pbc_with_parabola()
{
	ground_state_with_parabola(0);
}
int ground_state_obc_with_parabola()
{
	ground_state_with_parabola(1);
}




