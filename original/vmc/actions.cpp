#include "vmc.hpp"
//#define DEBUG

/* -- adjust Mcut analyzing elongations in critical configurations ---*/
void sistemavmc::adjust_mcut()
{
	Range I(0,Ihmax-2);
	Range J(1,Ihmax-1);
	int Mcut_new;
	
	find_critical();
	while((Mcut_new=max(abs(Is_critical(I)-Is_critical(J)))) >= Mcut)
	{
		cout << "Mcut_new=" << Mcut_new << endl;
		Mcut=Mcut_new+1;
		find_critical();
	}
	cout << "adapted Mcut=" << Mcut << endl; 
}

/* -- finds critical configuration and force (precision: SMALLNUMBER)---*/
void sistemavmc::find_critical()
{
	Real vmc_delta;
	//Real ff=F;
	//F=0.0;		
	//F=F_fraction;	
	
	for (int j=0;j<Ihmax;j++){ 
		Is(j) = (int)tilted(j);	
	}
	//cout << Is << endl;

	Is_critical = 0;
		
	//cout << "critical force calculation...." << endl;
		
	cout << "going backward..."; cout.flush();
	while((vmc_delta=VMC(BACKWARD))<0.0){ //backward
        }       
	cout << "done\n";

        #ifdef GRAPHICS
        cinema(0,BLUE,Is);
	tiempo++;
        #endif  

	cout << "going forward, F=" << F << "...\n";cout.flush();

	int count=0;
	F-=5.0*PRECISIONFC;
	while(count==0){
        	Is0=Is; 
		while(min(Is-Is0)<2*M)
        	{ //forward
		//cout << mean(Is-Is0) << " ";cout.flush();

			//si sube la energia para cualquier "elemental move"
			// aumento la fuerza impulsora
        		if((vmc_delta=VMC(FORWARD))>=0.0)
        		{
        			F+=(vmc_delta+PRECISIONFC);
                		Is_critical=Is;
		        	#ifdef GRAPHICS
                		cinema(2,RED,Is); // plot metastable
                		tiempo++;
	            		#endif
	                	cout << F << " " << mean(Is-Is0) << endl;
				count++;	 
        		}
            		#ifdef GRAPHICS
            		//cinema(2,BLUE,Is);          
            		#endif  
		}

		// si di la vuelta sin aumentar F, entonces F>Fc
		// es importante empezar con un buen ansatz...
		if(count==0) {
			F-=5.0*PRECISIONFC;cout << "initial F too large\n";
			Is = Is-2*M;
		};
	}
        Fc=F-PRECISIONFC;
	cout << "done\n";

}


/*---- finds the grouns state with periodic b.c. ----*/
void sistemavmc::find_gs(int mode)
{
	Real EG=BIGNUMBER, E; 
	Is_old = Is;
	
	Real ff=F;
	
	F=0;
	//cout << "GS with PBC" << endl;

	int max;
	#ifdef FIXEDTM	
	max=1;
	#else 
	max=M;
	#endif

	// finds sample ground state
	for(int i=0;i<max;i++){	
		Is=i;		

		Minimize_GS(0,mode);

		if(EG>(E=Energy()))
		{
			cout << i << " " << E << endl;
			EG=E;
			Is_gs_pbc=Is;
		}
	}		
	
	Is = Is_old;

	F=ff;
	/*
	#ifdef GRAPHICS
	cinema(2,8);
	#endif
	*/
}

/*---- given a configuration, relax it forward ----*/
void sistemavmc::find_next_meta()
{
	Real delta_min;
	while((delta_min=VMC(FORWARD,0,0))<=0.0){
	//while(VMC(FORWARD)<=0.0){
		#ifdef GRAPHICS 
		//cinema(2,2);          
		#endif		
		cout << delta_min << "\n";
	}
	//cout << "last VMC " <<  VMC(FORWARD)<< endl;
	//#ifdef GRAPHICS 
	//cinema(2,3);          
	//#endif
}

/*---- given a configuration, relax it forward ----*/
Real sistemavmc::find_next_meta(Real w, Real m2)
{
	Real delta_min;
	while((delta_min=VMC(FORWARD,w,m2))<=0.0){
		#ifdef GRAPHICS 
		//cinema(2,2);          
		#endif		
	}
	return delta_min;
	//cout << "last VMC " <<  VMC(FORWARD)<< endl;
	//#ifdef GRAPHICS 
	//cinema(2,3);          
	//#endif
}


/*----- finds escape path for configuration Is_alpha ------*/
/* 
finds: Is_beta, Is_gamma, Lopt, Lrel, Escape Barrier UB. 
arguments: Is_alpha, global biggest barrier (best_barrier), 
maximun avalanche size to limit path search (smax). 
*/
Real sistemavmc::find_Besc2(Real best_barrier, Array<int,1> &Is_alpha, 
	Array<int,1> &Is_beta, Array<int,1> &Is_gamma, int smax, int &Lopt, int &Lrel)
{
	Real Wmax,W;	
	int istart, idist,idist_stop;
	
	// water level imposed to paths. Set to the smallest among all good paths.
	Wmax = BIGNUMBER; 
		
#ifdef ACTIVE_POINTS
	Is_metastable_now = Is; 	
#endif
	
	Is_alpha = Is;
	Is_gamma = Is;
	
#ifndef INSTANTONIC_WALK // to target more rapidly Is_alpha with biggest UB
	idist=2;
	istart=0;	
	idist_stop=Ihmax;			
	do
	{		
		Is_limit = Is_alpha;
		
		// if a forward avalanche of size idist at istart can minimize energy ...
		if(Minimize_local_constrained(istart, idist, Is_limit, Is_limit,Is_limit,0)==1)
		{
			// idist_stop is set to the size of the smallest avalanche that minimize the
			// energy, plus smax. Bounded by Ihmax			
			if(idist_stop==Ihmax) idist_stop=(idist+smax+1>Ihmax)?(Ihmax):(idist+smax+1);			

			// initial metastable (alpha) configuration
			Is = Is_alpha;

			// finds path to escape, gets its minimal barrier W
			// search bounded by Wmax. If no paths available, W=BIGNUMBER 
			W=Besc_new2(Wmax);
			
			// if path barrier is the smallest  
			if(W < Wmax) 
			{
				Wmax=W;  
				Is_gamma = Is; // gamma
				Is_beta = Is_replic2; // beta				
				Lopt = count(Is_alpha != Is_beta); // Lopt
			}
		}	
			
		// next avalanche try 
		istart++;
		if(istart==Ihmax) {idist++; istart=0;}
							
	} while(Wmax >= best_barrier && idist < idist_stop); 	
	// this approximate dynamics is repeated while water 
	// level is higher than best global barrier
	// and forward avalanches minimizing the energy are available ...
#endif
	
#ifdef INSTANTONIC_WALK
idist=idist_stop; 
#endif
	
#ifdef EXACT
	// exact instanton dynamics, no limits. It is done only 
	// if last search can not improve best_barrier.
	if(idist==idist_stop && Wmax >= best_barrier)
	{
		Is = Is_alpha; 
		Is_limit = Is_alpha + M*5; // practically no limits

		W=Besc_new2(Wmax);
		
		if(W <= Wmax) 
		{
			Wmax=W; 
			Is_gamma = Is;    
			Is_beta = Is_replic2;				
			Lopt = count(Is_alpha != Is_beta);
		}	
	} 
#endif
	Is = Is_gamma; // new metastable state
	Lrel =  count(Is_gamma != Is_alpha);
	return Wmax; // new water level
}



/*-- find Besc of the configuration Is: uses set --*/
Real sistemavmc::Besc_new2(Real Wmax)
{
	Real UB, E0, E, E1;
	int inserted;
	set<key_Archive,cmp_Archive>::iterator It;
	
	UB=0.0;
	E0=Energy(); // energy of metastable state alpha: grand mother
		
	do{
		// generates all the configurations dynamically conected
		// with alpha by single VMC steps.		

		//cout << "Is=" << Is << endl;

		if(all(Is<=Is_limit)){		
			cout << "inserting confs" << endl;		
			inserted=grow_Tree_new2(FORWARD, Wmax);
		}
		cout << "inserted = " << inserted << endl;
		cout << "Archive size =" << Archive_set.size() << endl;
		
		//for(int j=1;j<Ihmax;j++) cout << ((It->u)[j-1]).to_ulong()-Mcut << endl;
				
		//-- copy and erase first element --> Is 
		//-- which is the less energetic children 
		//-- Archive is empty if Is can not reproduce
		if(!Archive_set.empty()){
			It = Archive_set.begin();	
			Is(0) = (It->i0); // copy h(0)
			for(int j=1;j<Ihmax;j++)
			{ 
#ifdef USING_CHAR_ARRAY	
				Is(j) = Is(j-1)+(It->u)[j-1]; // copy h(i+1)-h(i)
#endif	
#ifdef USING_BITSET_VECTOR	
				Is(j) = Is(j-1)+((It->u)[j-1]).to_ulong()-Mcut;
				//cout << "Is(j)=" << Is(j) << endl;
#endif	
			}			
			E  = It->E; // copy its energy
			
#ifdef USING_CHAR_ARRAY	
			free(It->u);
#endif	
#ifdef USING_BITSET_VECTOR	
			//cout << "erasing config vector of less energetic copied children" <<  endl;
			//Archive_set.erase(It);
			//for(int j=1;j<Ihmax;j++) ~((It->u)[j-1]);
			//cout << "OK" <<  endl;
#endif	
			cout << "erasing less energetic copied children from set ... " ;
			Archive_set.erase(It); // erase children from Archive
			cout << "OK" <<  endl;
			cout.flush();
		    #ifdef GRAPHICS
            //cinema(2,BLUE,Is);          
            #endif  

		}
		else return BIGNUMBER; // if mother can not reproduce ...

		// highest energy diference between less energetic (l.e.) 
		// children and mother is registered
		if(E-E0>UB) UB=E-E0; 
		
		//cout << "Is=" << Is << endl;
		//cout << "E=" << E << endl;
		Is_replic2 = Is; // save l.e. children before relaxing it to a mother.
		cout << "searching next metastable state ... " ; cout.flush();
		find_next_meta(); // children is now a metastable mother.		
		cout << "OK" << endl ; cout.flush();
	}
	while(Energy()>E0); 
	// Archive grows until a new mother (relaxed less energetic children) 
	// has a lower energy than the grand mother

#ifdef USING_CHAR_ARRAY	
	int memtotal=0;
	for(It = Archive_set.begin();It != Archive_set.end(); ++It)
	{
		free(It->u);
		//memtotal += sizeof(char)*(Ihmax-1)+sizeof(int);
	}
	//cout << memtotal << endl;	
	//system("ps -o %mem,comm | grep VMC");	
#endif	
#ifdef USING_BITSET_VECTOR	
/*	int memtotal=0;
	for(It = Archive_set.begin();It != Archive_set.end(); ++It)
	{
		for(int j=1;j<Ihmax;j++) ~((It->u)[j-1]);		
		//memtotal += sizeof(char)*(Ihmax-1)+sizeof(int);
	}*/
	//cout << memtotal << endl;	
	//system("ps -o %mem,comm | grep VMC");	
#endif	

	Archive_set.clear(); // clear all Archive

	// UB is the minimal barrier, Is is the next metastable state
	return UB;
}



/*---- finds the ground state with periodic b.c. ----*/
void sistemavmc::find_gs(int mode, int i0)
{
	Real EG=BIGNUMBER, E;
	Is_old = Is;

	Real ff=F;

	F=0;
	//cout << "GS with PBC" << endl;

	i0=Is_gs_pbc(0);

	// mass2*jmax^2 ~ 5 -> jmax ~ sqrt(5/mass2)

	// this is the criterion for stopping the search
	// the next metastable state can not be beyond i0+i1
	int i1=Ihmax;	

	// finds sample ground state
	//for(int i=(i0-5+M)%M;i<(i0+5)%M;i++){
	for(int i=i0;i<i0+i1;i++){
		Is=(i%M);

		Minimize_GS(0,mode);

		if(EG>(E=Energy()))
		{
			cout << i << " " << E << endl;
			EG=E;
			Is_gs_pbc=Is;
		}
	}

	Is = Is_old;

	F=ff;
	/*
	#ifdef GRAPHICS
	cinema(2,8);
	#endif
	*/
}



