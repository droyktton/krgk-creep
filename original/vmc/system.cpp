//#define INP
#include "vmc.hpp"
#include <exception>

sistemavmc::~sistemavmc()
{
#ifdef GRAPHICS	
	if(Ok!=0)
	g2_close(movie);
#endif
}

void sistemavmc::init()
{
	init(-1);
}

void sistemavmc::init(int mode)
{
	/*---- INPUT PARAMETERS ----*/
	#ifdef INP
	ifstream input("inp");
	input >> M >> Mcut >> Ihmax >> Adis >> F >> Ntot >> Cel >> beta >> shift >> Ok;     
	#else   
	cin >> M >> Mcut >> Ihmax >> Adis >> F >> Ntot >> Cel >> beta >> shift >> Ok;     
	#endif
	
	IHMAX = Ihmax;
	
	try{
		Is.resize(Ihmax);
		Is_replic.resize(Ihmax);
		Is_replic2.resize(Ihmax);
		Is0.resize(Ihmax);
		Is_old.resize(Ihmax);
		Is_critical.resize(Ihmax);
		Is_gs_pbc.resize(Ihmax);
		Is_biggest_barrier.resize(Ihmax);
		Is_biggest_barrier_excited.resize(Ihmax);
		Is_biggest_barrier_nextmeta.resize(Ihmax);
		Is_limit.resize(Ihmax);
		Is_metastable_now.resize(Ihmax);
	
		width_zero_T.resize(Ihmax);
		Is_diffusion.resize(Ihmax);
		Is_diffusion_inst.resize(Ihmax);
		Is_thermal_av.resize(Ihmax);
		Is_energy.resize(Ihmax);
		Is_energy2.resize(Ihmax);
		Is_entropy.resize(Ihmax);

		Ropt.init(Ihmax);
		Ropt2.init(Ihmax);
		tilted.resize(Ihmax);

		// super costosos...
		V.resize(Ihmax,M);	
		if(mode!=1){
			V_boltz.resize(Ihmax,M);	
			Is_density.resize(Ihmax,M);
		}
	}
	catch(exception& e)
  	{
    		cout << "Error allocating arrays: " << e.what() << endl;
  	}

	firstIndex i; tilted = i*shift;
	//cout << tilted << endl;
	
	// initial polymer
	Is=(int)(M*0.5); 

	Ropt.reset(1);        
	Ropt2.reset(1);        
	Fou.reset();  
	Is_diffusion=0.0;
	Is_thermal_av=0.0;
	Is_energy=0.0;
	Is_energy2=0.0;
	Is_entropy=0.0;
	Is_density=0;
	Archive_set.clear();
	mass_squared=0.0;

	// initialize Fourier transforms...
	Fou.init(Ihmax);
	
	// initialize random number generator
	// if no seed is given in idum.dat
	// seed = time
	unsigned long seed=0;
	ifstream seedfile("idum.dat");
	if(seedfile.good()) seedfile >> seed;
	else{ 
		cout << "RN seed=time" << endl;	
		seed = (unsigned long)time (NULL);
	}
	init_genrand(seed);	  
	cout << "seed = " << seed << endl;
	initial_seed=seed;	
	
	#ifdef OPTIMIZE5 
	// initialize tables of reductions IM and Ih
	IM_table.resize(Range(-M*Mcut,8*M*Mcut));
	for(int im=-M*Mcut;im<=8*M*Mcut;im++){
		IM_table(im)=((im+100*M)%M);
		//cout << im << " " << IM_table(im) << endl;
	}
	cout << "IM initialized" << endl;
		
	Ih_table.resize(Range(-Ihmax,2*Ihmax));
	for(int ih=-Ihmax;ih<=2*Ihmax;ih++)
	Ih_table(ih)=((ih+Ihmax)%Ihmax);
	cout << "Ih initialized" << endl;
	#endif									
	
	/*-- genera Tri Net, used by transfer matrix minimize --*/
	genera_Net();
	
	#ifdef GRAPHICS	
	// movie	
	if(Ok!=0){
		int sizey = 500; //int(sizex*Ihmax*(1.0/M));	
		int sizex = int(sizey*M*(1.0/Ihmax));
		#ifdef EPS
		//int movie_eps = g2_open_EPSF_CLIP ("movie.ps",sizex,sizey);  
		  //movie = g2_open_EPSF ("movie.ps");  
		  //movie_PNG = g2_open_gd(msg, (int)(sizex), (int)(sizey), g2_gd_png); 
		#else 
		//movie = g2_open_X11(size,size);
		movie=g2_open_vd();
		int movie_X11=g2_open_X11X(sizex,sizey,0,0,(char *)"VMC creep", (char *)"VMC icono",
				(char *)"VMC data",5,5);
		//movie_eps = g2_open_EPSF_CLIP ("movie.eps",sizex,sizey);

		g2_attach(movie,movie_X11);
		//g2_attach(movie,movie_eps);
        	g2_set_coordinate_system(movie, 0, 0, sizex*1.0/M, sizey*1.0/Ihmax);
		#endif
	}
	#endif

	Difi.open("Difi.dat");
	out_conf_critical.open("conf_critical.dat");
	out_conf_biggest_barrier.open("conf_biggest_barrier.dat");
	out_conf_biggest_barrier_excited.open("conf_biggest_barrier_excited.dat");
	out_conf_biggest_barrier_nextmeta.open("conf_biggest_barrier_nextmeta.dat");
	out_conf_ground_state.open("conf_ground_state.dat");
	out_dis.open("disorder_6_4.data");
	out_meta.open("metastable_states.dat");
	
	logfile.open("logfile.dat");
	
	F_fraction = F;
	//F=0.0;
	
	logfile << "M=" << M << endl;
	logfile << "Mcut=" << Mcut << endl;
	logfile << "Ihmax=" << Ihmax << endl;
	logfile << "Adis=" << Adis << endl;
	#ifndef ABSOLUTEDRIVENFORCE
	logfile << "F/Fc=" << F_fraction << endl;
	#else
	logfile << "F [absolute]=" << F_fraction << endl;
	#endif
	logfile << "Ntot=" << Ntot << endl;
	logfile << "Cel=" << Cel << endl; 
	logfile << "beta=" << beta << endl; 
	logfile << "shift=" << shift << endl; 	
	logfile << "Movie=" << Ok << endl; 	
	logfile << "RN seed=" << seed << endl; 	

	#ifdef RANDOM_FIELD
	logfile << "random field" << endl;
	#else
	logfile << "random bond" << endl;
	#endif

}

void sistemavmc::genera_Net()
{
	// initialize Tri-network  
	
	int id; Real dis;
	try{
	Net.resize(INet(Mcut*(Ihmax-1),(Ihmax-1))+1);
	}
	catch(exception& e)
  	{
    		cout << "Error allocating Net: " << e.what() << endl;
  	}
	cout << "size of Net= " << Mcut*(Ihmax-1)<< endl;
	char * nomnet = (char *)malloc(sizeof(char)*100);
	sprintf(nomnet,"Net_Mcut%d_beta%f_Ihmax%d_Cel%f.dat",Mcut,beta,Ihmax,Cel);
	ifstream fNet(nomnet);
	
	if(fNet.good())
	{
		int a1; Real a2, a3;
		for(int i=0;i<INet(Mcut*(Ihmax-1),(Ihmax-1))+1;i++)
		{
			fNet >> Net(i).j >> Net(i).k >> Net(i).E >> Net(i).best >>  Net(i).n;			
			for(int j=0;j<Net(i).n;j++){
				fNet >> a1 >> a2 >> a3;
				(Net(i).list_below).push_back(a1);
				(Net(i).bond_below).push_back(a2);
			        (Net(i).bond_boltz_below).push_back(a3);
			}		
		}	
		cout << "Net initialized from file ..." << endl;	
	}
	
	else
	{
		for (int k=1;k<Ihmax;k++){
		for (int j=-Mcut*k ;j<=Mcut*k;j++){

			int i=INet(j,k);
			Net(i).j=j;
			Net(i).k=k;
			Net(i).E=BIGNUMBER;
			Net(i).best=0;		

			if(k == 1){
				 Net(i).n=1;
				/*(Net(i).list_below).resize(1); 
				(Net(i).bond_below).resize(1); 
				(Net(i).list_below)(0) = 0; 								
				(Net(i).bond_below)(0) = Cel*0.5*pow2((Real)Net(i).j); 								
				*/
				(Net(i).list_below).push_back(0); 								
				(Net(i).bond_below).push_back(Cel*0.5*pow2((Real)Net(i).j)); 								
				(Net(i).bond_boltz_below).push_back(exp(-beta*Cel*0.5*pow2((Real)Net(i).j))); 								
			}
			else{
				/*(Net(i).list_below).resize(2*Mcut+1); 
				(Net(i).bond_below).resize(2*Mcut+1);*/ 
				Net(i).n=0;
				for(int jj=-Mcut*(k-1);jj<=Mcut*(k-1);jj++){					
					id  = INet(jj,k-1); 
					dis = Net(i).j - Net(id).j;
					if(abs(dis)<=Mcut){	
						/*(Net(i).list_below)(Net(i).n) = id; 	
						(Net(i).bond_below)(Net(i).n) = Cel*0.5*pow2(dis);*/ 											
						(Net(i).list_below).push_back(id); 	
						(Net(i).bond_below).push_back(Cel*0.5*pow2(dis)); 											
						(Net(i).bond_boltz_below).push_back(exp(-beta*Cel*0.5*pow2(dis))); 											
						Net(i).n++;
					}
				}				
			}

		}}

		//-- origin
		Net(0).j=0;
		Net(0).k=0;
		Net(0).E=0.0;
		Net(0).best=0;		
		/*(Net(0).list_below).resize(2*Mcut+1); 
		(Net(0).bond_below).resize(2*Mcut+1);*/ 
		Net(0).n=0;  
		for(int jj=-Mcut*(Ihmax-1);jj<Mcut*(Ihmax-1);jj++){					
				id  = INet(jj,Ihmax-1);  
				dis = Net(0).j - Net(id).j;
				if(abs(dis)<=Mcut){	
					/*(Net(0).list_below)(Net(0).n) = id; 	
					(Net(0).bond_below)(Net(0).n) = Cel*0.5*pow2(dis);*/
					(Net(0).list_below).push_back(id); 	
					(Net(0).bond_below).push_back(Cel*0.5*pow2(dis));				 											
					(Net(0).bond_boltz_below).push_back(exp(-beta*Cel*0.5*pow2(dis))); 											
					Net(0).n++;
				}
		}		


		/*
			fNet >> Net(i).j >> Net(i).k >> Net(i).E >> Net(i).best >>  Net(i).n;			
			for(int j=0;j<Net(i).n;j++){
				fNet >> a1 >> a2 >> a3;
				(Net(i).list_below).push_back(a1);
				(Net(i).bond_below).push_back(a2);
			        (Net(i).bond_boltz_below).push_back(a3);
			}		
		*/
		ofstream fNetout(nomnet);
		for(int i=0;i<INet(Mcut*(Ihmax-1),(Ihmax-1))+1;i++)
		{
			fNetout << Net(i).j << " " << Net(i).k << " " << Net(i).E << " " << Net(i).best << " " << Net(i).n << endl;
			for(int j=0;j<Net(i).n;j++){
			fNetout << Net(i).list_below[j] << " " << Net(i).bond_below[j] << " " << Net(i).bond_boltz_below[j] << endl;
			}		
		}
		
		cout << "Net initialized from scratch ..." << endl;
	}

	// testing Tri
/*	int i,j,k,nei;
	cout << "Test Tri: (j,k)=?" << endl;
	j=0;
	while(1){
	cin >> j >> k;
	if(j==111) exit(1);
	i=INet(j,k);
	cout << "chosen: (j,k)=" << "(" << j <<","<< k <<")"<< endl;
	cout << "Net-->(j,k)=" << "(" << Net(i).j <<","<< 
			 Net(i).k <<")"<< endl;
	cout << "Net--> E=" << Net(i).E << endl;	
	cout << "Net--> best=" 	<< Net(i).best << endl;
	cout << "Net--> # neigh_below=" << Net(i).n << endl;
	for(int l=0;l<Net(i).n;l++){
		nei=Net(i).list_below(l);
		cout << "Net--> neigh = (" << Net(nei).j << "," 
		<< Net(nei).k << ")" << " ;  neigh.E =" << Net(nei).E << endl;		
	}
	}	
*/

/*

-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6
	  -4 -3 -2 -1 0 1 2 3 4
	        -2 -1 0 1 2
		          0

-4 -3 -2 -1 0 1 2 3 4
   -3 -2 -1 0 1 2 3
      -2 -1 0 1 2 
         -1 0 1
            0

*/
}

/*-- relatives coordinates (j,k) --> id Net number */
int sistemavmc::INet(int j, int k)
{
	if(j>Mcut*k){cout << "error:  j > Mcut k\n"; exit(1);}
	if(j<-Mcut*k){cout << "error:  j < -Mcut k\n"; exit(1);}
	if(k>=Ihmax){cout << "error:  k >= Ihmax\n"; exit(1);} 
	if(k<0){cout << "error:  k < 0\n"; exit(1);} 
	
	return k*k*Mcut+k+j;
}


//////////////////////////////////////////////////////////////////////////////////////////
/*--- generates disorder --- */

#ifndef PSIX
#define PSIX 1
#endif
#ifndef PSIY
#define PSIY 1
#endif

void sistemavmc::correlate_disorder()
{
	fftw_complex * Din;
	fftw_complex * Cin;
	fftw_complex * Uout;
	
	fftw_plan Dp, Cp, Up;


	Din = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ihmax*M);
	Dp = fftw_plan_dft_2d(Ihmax, M, Din, Din, FFTW_FORWARD, FFTW_MEASURE);

	Cin = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ihmax*M);
	Cp = fftw_plan_dft_2d(Ihmax, M, Cin, Cin, FFTW_FORWARD, FFTW_MEASURE);

	Uout = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ihmax*M);
	Up = fftw_plan_dft_2d(Ihmax, M, Uout, Uout, FFTW_BACKWARD, FFTW_MEASURE);


	int ind, di, dj;

	// cargo desorden en formato de fftw
	for(int i=0; i<Ihmax; i++){
		for(int j=0; j<M; j++){
			ind=M*i+j;
			Din[ind][0]=V(i,j);
			Din[ind][1]=0.0;

			di=(i<=Ihmax*0.5)?(i):(Ihmax-i); 
			dj=(j<=M*0.5)?(j):(M-j); 

			Cin[ind][0]= exp(-di/PSIY)*exp(-dj/PSIX); //short range			
			Cin[ind][1]= 0.0;
	}}
	fftw_execute(Cp);// fft de la funcion de correlacion 
	fftw_execute(Dp);// fft de los numeros descorrelacionados 

	for(int i=0;i<Ihmax;i++){
	for(int j=0;j<M;j++){
		ind=M*i+j;
		Uout[ind][0]=sqrt(fabs(Cin[ind][0]))*Din[ind][0];
		Uout[ind][1]=sqrt(fabs(Cin[ind][0]))*Din[ind][1];
	}}
	fftw_execute(Up);

	for(int i=0;i<Ihmax;i++){
	for(int j=0;j<M;j++){
		ind=M*i+j;
		V(i,j)=Uout[ind][0];
	}}
	V = V/sqrt(mean(V*V));
	
	fftw_destroy_plan(Dp);
	fftw_destroy_plan(Cp);
	fftw_destroy_plan(Up);
    
	fftw_free(Cin); 
	fftw_free(Din); 
	fftw_free(Uout);
	//cout << V << endl;

}

void sistemavmc::set_parabola(double mass2,double w)
{
#ifndef NEWVERSION
	// not optimal: we are generating the same disorder
	init_genrand(initial_seed);
	change_disorder();

	// not optimal: we add a parabola
	#ifndef TIPDRIVENSTATICS
	for(int i=0; i<Ihmax; i++){
	#else
	for(int i=0; i<1; i++){
	#endif
	for(int j=0; j<M; j++){
		double dw=(w-j);
		if(dw>M*0.5) dw-=M;
		if(dw<-M*0.5) dw+=M;
		#ifndef TIPDRIVENSTATICS
		V(i,j)+=mass2*dw*dw*0.5;
		#else
		V(i,j)+=5*dw*dw*0.5;
		#endif
	}}
#else
	parabola_position=w;
	mass_squared=mass2;
#endif
}

void sistemavmc::change_disorder()
{
	out_dis.precision(6);
	out_dis.width(12);
	out_dis.flags(ios_base::fixed);

	//ifstream inp_dis("desorden_JL.dat");

	for(int i=0; i<Ihmax; i++){
#ifdef RANDOM_FIELD
	   V(i,0)=0.; // starting point
	   V(i,M-1)=0.0; // end point
	   double beta=1.0; // gaussian width
	   double zend=V(i,M-1);
	   for( int jj = 1; jj <= M-2; jj++){
		double zinit= V(i,jj-1);
		double tau_1=1./beta;
		double tau_2=(M-1-jj)/beta;
		double zmean=(tau_2* zinit + tau_1 * zend)/(tau_1 + tau_2);
		double sigmasq=1./( 1./tau_1 + 1./tau_2);
		double sigma= sqrt(sigmasq);
		V(i,jj)=zmean + gasdev_real3()*sigma;
	   }
/*
	   V(i,0)=0.0;
	   for(int j=1; j<M; j++){
	  	   V(i,j) = V(i,j-1) + gasdev_real3();
	   }
*/
#elif defined(PERIODIC_DISORDER)
	   int fac=(int)(M*0.25); // divide la muestra en fac bloques
	   int m = (int)(M*1.0000001/fac); // tamanio de bloque
	   for(int j=0; j<m; j++){
		double nro= genrand_real3();
		for(int n=0;n<fac;n++) V(i,(j+n*m)%M)=nro;
	   }	
#elif defined(PERIODIC_CHAIN)
	   int fac=Ihmax; // divide la muestra en L pedazos
	   int m = (int)(M*1.0000001/fac);
	   for(int j=0; j<M; j++){
		if(i==0){ 
		  double nro= genrand_real3();
		  V(i,j)=nro;
		}
		else{
		  int jj = j + (int)(i*M/Ihmax);
		  V(i,j)=V(0,jj%M);		  
		}  
	   }		   
#else
	   for(int j=0; j<M; j++){
		V(i,j)=gasdev_real3();
	   }  
#endif
	}


#ifdef RANDOM_FIELD
	   Array<Real,1> aux(M);
	   int shift=0;//(int)(genrand_real2()*M);
	   for(int i=0; i<Ihmax; i++){
	   aux(Range(0,M-1))=V(i,Range(0,M-1)); 
	   for(int j=0; j<M; j++){
		V(i,j)=aux((j+shift)%M);
	   }}

	  /*for(int j=0; j<M; j++){
		out_dis << V(0,j) << endl;	      
	  }
	  exit(1);*/
#endif

#ifdef STRONG_PINNING
	V=0.0;
	unsigned int Np=(int)(M*Ihmax*0.005);
	int ix, iy;
	for(int n=0;n<Np;n++){
		ix =(int)(genrand_real3()*M);
		iy =(int)(genrand_real3()*Ihmax);
		V(iy,ix)=-1.;
		//V(iy,(ix-1+M)%M)=V(iy,(ix+1)%M)=-2;
	}
#endif

	// print disorder in a file
	for(int i=0; i<Ihmax; i++){
		for(int j=0; j<M; j++){
			//inp_dis >> V(i,j);
			/*#ifdef PERIODIC_DISORDER
			V(i,j)= (1+gasdev_real3())*cos(M_PI*j*0.5);
			//V(i,j)=gasdev_real3();
			#endif*/
			#ifdef PINNING_CENTER
			V(i,j)= ((j-2)*(j-2)+(i-2)*(i-2)<10)?(1.0):(0.0);
			//V(i,j)=gasdev_real3();
			#endif
			//out_dis << V(i,j) << endl;
		}
		//out_dis << endl << endl;
	}
	//exit(1);

/*	another way to write it : used for testing
	for(int i=0; i<Ihmax; i++){
	for(int j=0; j<M; j++){
		V(i,j)=gasdev_real3();
		out_dis << V(i,j) << endl;
	}}
*/

	V=(V-mean(V));
	//V=V/sqrt(mean(V*V));

	V*=Adis;
	V_boltz=exp(-beta*V);	
	Archive_set.clear();
}

void sistemavmc::copy_disorder(sistemavmc &Q, double delta)
{
	Archive_set.clear();

#ifndef RANDOM_FIELD
	for(int i=0; i<Ihmax; i++){
		for(int j=0; j<M; j++){
			V(i,j)= Q.V(i,j) + gasdev_real3()*delta;
		}
	}
	V=(V-mean(V));
	V*=1./sqrt(1.0+delta*delta);
	//V=V/sqrt(mean(V*V));
	cout << "mean=" << mean(Q.V) << " , " << mean(V) << endl;
	cout << "variance=" << mean(Q.V*Q.V) << " , " << mean(V*V) << endl;
#else
	change_disorder();
	for(int i=0; i<Ihmax; i++){
		for(int j=0; j<M; j++){
			V(i,j)= Q.V(i,j) + V(i,j)*delta;
		}
	}
	//V=(V-mean(V));
	V*=1./sqrt(1.0+delta*delta);
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////



/*-- properties functions --*/

//--- energy of Is only
Real sistemavmc::Energy()
{
	Real E=0.0;
	for(int i=0;i<Ihmax-1;i++){ 
		#ifndef NEWVERSION
		E = E + V(i,IM(Is(i))) - F*Is(i) 
		      + Cel*0.5*(Is(Ih(i+1))-Is(i))*(Is(Ih(i+1))-Is(i));	
		#else
		E = E + potV(i,Is(i)) - F*Is(i) 
		      + Cel*0.5*(Is(Ih(i+1))-Is(i))*(Is(Ih(i+1))-Is(i));	
		#endif
	}		

	#ifndef NEWVERSION
	E = E + V(Ihmax-1,IM(Is(Ihmax-1))) - F*Is(Ihmax-1) 
	      + Cel*0.5*pow2(Is(0)-Is(Ihmax-1)+shift*Ihmax);	
	#else
	E = E + potV(Ihmax-1,Is(Ihmax-1)) - F*Is(Ihmax-1) 
	      + Cel*0.5*pow2(Is(0)-Is(Ihmax-1)+shift*Ihmax);	
	#endif
	return E;
}


//--- energy of a given configuration
Real sistemavmc::Energy_conf(Array<int,1> X)
{
	Real E=0.0;
	for(int i=0;i<Ihmax;i++){ 
		#ifndef NEWVERSION
		E = E + V(i,IM(X(i))) - F*X(i) 
		      + Cel*0.5*(X(Ih(i+1))-X(i))*(X(Ih(i+1))-X(i));	
		#else
		E = E + potV(i,X(i)) - F*X(i) 
		      + Cel*0.5*(X(Ih(i+1))-X(i))*(X(Ih(i+1))-X(i));	
		#endif	
	}		
	return E;
}

//--- energy of a given configuration
Real sistemavmc::Width_conf(Array<int,1> X)
{
	Real CM = mean(X);
	Real width = mean((X-CM)*(X-CM)); 
	return width;
}


//--- energy of a given configuration
Real sistemavmc::Energy_elastic(Array<int,1> X)
{
	Real E=0.0;
	for(int i=0;i<Ihmax;i++){ 
		E = E + Cel*0.5*pow2(X(Ih(i+1))-X(i));	
	}		
	return E;
}

//--- energy of a given configuration
Real sistemavmc::Energy_pinning(Array<int,1> X)
{
	Real E=0.0;
	for(int i=0;i<Ihmax;i++){ 
		E = E + V(i,IM(X(i)));	
	}		
	return E;
}


/*---------------------------------*/ 

void sistemavmc::cinema_set_png_name()
{
#ifdef GRAPHICS	
	int sizex = 700;
	int sizey = int(sizex*Ihmax*(1.0/M));	
	char nom[50];
	sprintf(nom,"peli_%d.png",100000000+tiempo);
	printf("%s\n",nom);
	movie_png = g2_open_gd(nom, sizex, sizey, g2_gd_png);
        g2_set_coordinate_system(movie_png, 0, 0, sizex*1.0/M, sizey*1.0/Ihmax);
	g2_attach(movie,movie_png);
#endif	
}

void sistemavmc::cinema_close_png_name()
{
#ifdef GRAPHICS	
	g2_detach(movie,movie_png);
	g2_close(movie_png);
#endif
}


void sistemavmc::cinema(int flag_erase, int color, Array<int,1> X)
{
#ifdef GRAPHICS	
	
	if(Ok!=0){
		int i;
		Real x1,y1,x2,y2;
		struct timespec delay;
		#ifndef MOVIE_PNG
		delay.tv_sec=0;
 	   	delay.tv_nsec=500000;
		#else
		delay.tv_sec=0;
 	   	delay.tv_nsec=0;
		#endif

		if(flag_erase==0) g2_clear(movie);
		for(i=0; i<Ihmax+1; i++){			
			if(flag_erase==1){
			for(int j=0;j<M+1;j++)
			{
				g2_pen(movie,(V(Ih(i),IM(j))>0.0)?(3):(7));
				g2_circle(movie, j, i, fabs(V(Ih(i),IM(j)))/(6.0*Adis));
			}}	

			g2_pen(movie,color);	
			x1= IM(X(Ih(i))); y1=i;		
			x2= IM(X(Ih(i+1))); y2=i+1;

			g2_filled_circle(movie, x1, y1, 0.1);
			if(fabs(x1-x2)<M*0.5)
			g2_line(movie, x1, y1, x2, y2);
		}
		//nanosleep(&delay,NULL);
		//sleep(1);
		tiempo++;
	}
#endif
}

void sistemavmc::cinema(int flag_erase, int color, Array<double,1> X)
{
#ifdef GRAPHICS	
	
	if(Ok!=0){
		int i;
		Real x1,y1,x2,y2;
		struct timespec delay;
		#ifndef MOVIE_PNG
		delay.tv_sec=0;
 	   	delay.tv_nsec=500000;
		#else
		delay.tv_sec=0;
 	   	delay.tv_nsec=0;
		#endif

		if(flag_erase==0) g2_clear(movie);
		for(i=0; i<Ihmax+1; i++){			
			if(flag_erase==1){
			for(int j=0;j<M+1;j++)
			{
				g2_pen(movie,(V(Ih(i),IM(j))>0.0)?(3):(7));
				g2_circle(movie, j, i, fabs(V(Ih(i),IM(j)))/(6.0*Adis));
			}}	

			g2_pen(movie,color);	
			x1= IM((int)X(Ih(i))); y1=i;		
			x2= IM((int)X(Ih(i+1))); y2=i+1;

			g2_filled_circle(movie, x1, y1, 0.1);
			if(fabs(x1-x2)<M*0.5)
			g2_line(movie, x1, y1, x2, y2);
		}
		//nanosleep(&delay,NULL);
		//sleep(1);
		tiempo++;
	}
#endif
}


/*void sistemavmc::fftw(int flag)
{
	Fou.faire(Is);
	if(flag==0){
		Fou.S=Fou.S+Fou.Si;
		Fou.N++;
	}
	else
	{
		Fou.S2=Fou.S2+Fou.Si;
		Fou.N2++;	
	}
}*/

////////////////////////////////////////////////////////////////////////////////////
// Fourier transforms
void Fourier_transform::print(const char * filename)
{	
	int size=S.size();
	ofstream Struc(filename);

	Struc << "#####" << endl;
	for(int i=0;i<size*0.5;i++)
	{
		//Struc << i*2.0*M_PI/size << " " << S(i)/N << " " << S2(i)/N2 << " " << Smaxtot(i)/Ntotal << endl;
		Struc << i*2.0*M_PI/size << " " << S(i)/N << " " << endl;

	}
}

void Fourier_transform::reset()
{	
	Si=0.0;
	N=0;
	S=0.0;
	N2=0;
	S2=0.0;
	Smax=0.0;
	Smaxtot=0.0;	
}

void Fourier_transform::init(int d)
{	
	Sin = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * d);

	Sp = fftw_plan_dft_1d(d, Sin, Sin, FFTW_FORWARD, FFTW_MEASURE);

	Si.resize(d);
	S.resize(d);
	S2.resize(d);
	Smax.resize(d);
	Smaxtot.resize(d);
	Rei.resize(d);
	Imi.resize(d);
}

void Fourier_transform::faire(Array<double,1> Z)
/* structure factor should be:
S(i) = (Sin[i][0]*Sin[i][0]+Sin[i][1]*Sin[i][1]);
*/
{
	for(int i=0;i<Z.size();i++)
	{
		Sin[i][0]=Z(i); Sin[i][1]=0.0;
	}
	//cout << "executing fftw" << endl;
	fftw_execute(Sp);	
	//cout << "ready" << endl;

	for(int i=0;i<Z.size();i++)
	{
		Rei(i)=Sin[i][0];
		Imi(i)=Sin[i][1];
		Si(i)=(Sin[i][0]*Sin[i][0]+Sin[i][1]*Sin[i][1]);
		S(i)=S(i)+Si(i);
	}
	N++;
	
}

void Fourier_transform::faire(Array<double,1> &Z1, Array<double,1> &Z2)
/* structure factor should be:
S(i) = (Sin[i][0]*Sin[i][0]+Sin[i][1]*Sin[i][1]);
*/
{
	for(int i=0;i<Z1.size();i++)
	{
		Sin[i][0]=Z1(i); Sin[i][1]=0.0;
	}
	fftw_execute(Sp); // 1st transform

	for(int i=0;i<Z1.size();i++)
	{
		Rei(i)=Sin[i][0];
		Imi(i)=Sin[i][1];
	}

	for(int i=0;i<Z2.size();i++)
	{
		Sin[i][0]=Z2(i); Sin[i][1]=0.0;
	}
	fftw_execute(Sp); //transform

	for(int i=0;i<Z2.size();i++)
	{
		Si(i)=(Sin[i][0]*Rei(i)+Sin[i][1]*Imi(i));
		S(i)=S(i)+Si(i);
	}
	N++;
}

////////////////////////////////////////////////////////////////////////////////////
// avalanches ...
void avalanche::init(int size)
{
	histo_r.resize(size+1);
	histo_a.resize((size+1)*(size+1));
	histo_ritoz.resize(size+1);
}

void avalanche::acum(int x,int y)
{
	ri=x;
	ai=y;
	r+=x;
	a+=y;
	ritozi=1.0*y;
	
	N++;
	Ntot++;
	histo_r(ri)++;
	histo_a(ai)++;		
	histo_ritoz(ri)+= ritozi;	
}

void avalanche::reset(int flag)
{
	ri=r=ai=a=N=0;	

	if(flag==1){
		Ntot = 0;
		histo_a = 0;
		histo_r = 0;
		histo_ritoz = 0.0;
	}
}

void sistemavmc::print_avalanches(int flag)
{	
	//ofstream Ava("Avalanches.dat");
	
	if(flag==0){	
		FILE * Ava=fopen("Avalanches.dat","a");
		//Ava << r*1.0/N << " " << a*1.0/N << endl;
		fprintf(Ava,"%d %lf %lf %d %lf %lf\n",Ropt.N, Ropt.r*1.0/Ropt.N,Ropt.a*1.0/Ropt.N, 
								   Ropt2.N, Ropt2.r*1.0/Ropt2.N,Ropt2.a*1.0/Ropt2.N);
		fclose(Ava);
	}
	
	if(flag==1){
		ofstream Roptic("Ropt.dat");
		for(int i=0;i<Ropt.histo_r.size();i++)
		{
				if(Ropt.histo_r(i)!=0 || Ropt2.histo_r(i)!=0)
				Roptic << i << " " << (Ropt.histo_r(i)*1.0/Ropt.Ntot) 
							<< " " << (Ropt2.histo_r(i)*1.0/Ropt2.Ntot) <<endl;
		}

		ofstream Aoptic("Aopt.dat");
		for(int i=0;i<Ropt.histo_a.size();i++)
		{
				if(Ropt.histo_a(i)!=0 || Ropt2.histo_a(i)!=0)
				Aoptic << i << " " << (Ropt.histo_a(i)*1.0/Ropt.Ntot) << " " << (Ropt2.histo_a(i)*1.0/Ropt2.Ntot) << endl;
		}

		ofstream Ritoz("Ritoz.dat");
		for(int i=0;i<Ropt.histo_ritoz.size();i++)
		{
				if(Ropt.histo_ritoz(i)!=0)
				Ritoz << i << " " << (Ropt.histo_ritoz(i)*1.0/(i*Ropt.histo_r(i))) << 
							  " " << (Ropt2.histo_ritoz(i)*1.0/(i*Ropt2.histo_r(i))) << endl;
		}
	}

}


void sistemavmc::print_diffusion(int N)
{
	ofstream Dif("Dif.dat");
	for(int i=0;i<Ihmax;i++)
	{
		Dif << i << " " << Is_diffusion(i)*1.0/N << " " << (Is_energy2(i)*1.0/N-pow2(Is_energy(i)*1.0/N)) << 
		" " << Is_thermal_av(i)*1.0/N << " " << -Is_entropy(i)*1.0/N << endl;	
	}
}

void sistemavmc::print_diffusion_inst()
{
	for(int i=0;i<Ihmax;i++)
	{
		Difi << Is_diffusion_inst(i) << endl;
	}
	Difi << endl;	
}

void sistemavmc::Density()
{
	for(int i=0;i<Ihmax;i++) 
	Is_density(i,IM(Is(i)))+=1;
}

void sistemavmc::print_Density()
{
	ofstream Den("Den.dat");
	for(int i=0;i<Ihmax;i++){ 
		Den << Is_density(i,0);
		for(int j=1;j<M;j++){ 
			Den << " " << Is_density(i,j);
		}
		Den << endl;
	}
}

void dprop::acum(Real x)
{
	vi=x;
	v+=x;
	N++;
}

void dprop::reset()
{
	vi=v=0.0; N=0;	
}


void sistemavmc::print_config(Array<int,1> X)
{
	for(int i=0;i<Ihmax;i++){ 
		cout << IM(X(i)) << endl;	
	}
}

void sistemavmc::print_conf(ofstream &out, Array<int,1> X)
{
	//ofstream out(nom);
	out << "#####" << endl;
	for(int i=0;i<Ihmax;i++){ 		
		out << i << " " << IM(X(i)) << " "<< X(i) << endl;	
	}
	out << endl << endl;
	out.flush();
}
void sistemavmc::print_conf(ofstream &out, Array<double,1> X)
{
	//ofstream out(nom);
	out << "#####" << endl;
	for(int i=0;i<Ihmax;i++){ 		
		out << i << " " << IM((int)X(i)) << " "<< X(i) << endl;	
	}
	out << endl << endl;
}

void sistemavmc::print_conf2(Array<int,1> X)
{
	FILE *fmeta=fopen("meta.dat","w");
	fprintf(fmeta,"%12.6f\n",Energy_conf(X));
	for(int i=0;i<Ihmax;i++){ 		
		fprintf(fmeta,"%d\n",X(i));	
	}
	fclose(fmeta);
}

void sistemavmc::check_mcut(Array<int,1> X)
{
	int dif, maxdif=0;
	for(int i=1;i<Ihmax;i++){
		if((dif=abs(X(i-1)-X(i)))>Mcut){
			cout << "error Dh = " << dif << " > " << Mcut << endl;
			exit(1);
		}
		if(maxdif<dif) maxdif=dif;
	}
	cout << "maxdif=" << maxdif << endl;
}



