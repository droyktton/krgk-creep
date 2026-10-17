#include "vmc.hpp"

/* Transfer matrix construction of the ground state 
of a polymer constrained in one point with periodic 
(mode=0) or open (mode=1) boundary conditions
*/
int sistemavmc::Minimize_GS(int measure, int mode)
{
	int iflag;
	int idist, istart, mstart,mend,mb;
	int j,k,id,l,idn,idb,kk, kmin, kmax, kaux;
	Real Emin,Ebond,Edis;		
	Real norm;
	Net(0).w=1.0;

	Is_replic=Is;	
	iflag=0;	
	Is_diffusion_inst=0.0;
	width_zero_T=0;

	idist=Ihmax-1; // size of the avalanche
	istart=0;      // starting fixed point for the avalanche

	//-- end points
	mstart=Is(istart);
	mend=mstart + (int)(shift*Ihmax);


	/* --------------------------------------------------------
	   level by level construction of minimal energy polymers 
	-----------------------------------------------------------*/
	for(k=1;k<=idist;k++)
	{
		if(mode==0){
			//-- "two cones" constraint: PBC
			kmin=( (-k*Mcut) >= (kaux=Mcut*(k-(idist+1))+mend-mstart) )?(-k*Mcut):(kaux);
			kmax=( k*Mcut <= (kaux=Mcut*((idist+1)-k)+mend-mstart) )?(k*Mcut):(kaux);		
		}

		if(mode==1){
			// one cone constraint: OBC
			kmin=(-k*Mcut); 
			kmax=( k*Mcut);				
		}				
		
		norm=0.0;
		for(j=kmin;j<=kmax;j++)
		{
			id=INet(j,k); //-- node (k,j)==(z,x) 
			Net(id).w=0.0;

			Emin=BIGNUMBER; 
			for(l=0;l<Net(id).n;l++){//-- conected nodes below
				idn = Net(id).list_below[l];				
				
				Ebond= Net(id).bond_below[l];
				if(beta>0.0) Net(id).w+= Net(id).bond_boltz_below[l] * Net(idn).w;//exp(-beta*Ebond)*Net(idn).w;

				if((Net(idn).E+Ebond) <Emin){ //-- polymer energy up to node <idn>
					Emin=Net(idn).E+Ebond;
					Net(id).best=idn;
				}
			}

			//-- energy of polymer up to the node <id>=(k,j)										
			#ifndef NEWVERSION
			Net(id).E  = Emin + V(Ih(istart+k),IM(mstart+j));
			#else
			Net(id).E  = Emin + potV(istart+k,mstart+j);						
			#endif						
													    
			if(beta>0.0){
				#ifndef NEWVERSION
				Net(id).w *= V_boltz(Ih(istart+k),IM(mstart+j));//exp(-beta*V(Ih(istart+k),IM(mstart+j)));
				#else
				Net(id).w *= potV_boltz(istart+k,mstart+j);//exp(-beta*V(Ih(istart+k),IM(mstart+j)));
				#endif
				norm += Net(id).w;
			}
		}// transverse at level k


		if(beta>0.0){
			Real u_k=0.0; 
			for(j=kmin;j<=kmax;j++){ 
				id=INet(j,k);
				Net(id).w = Net(id).w/norm; 
				Is_diffusion(k)+=((j*j)*Net(id).w);			
				Is_diffusion_inst(k)+=((j*j)*Net(id).w);			
				u_k += (j*Net(id).w);			
				Is_energy2(k)+=(pow2(Net(id).E)*Net(id).w);			
				Is_energy(k)+=((Net(id).E)*Net(id).w);			
				//Is_entropy(k)+=(Net(id).w*log(Net(id).w));			
			}	
			Is_thermal_av(k)+= (u_k*u_k);
		}
		else{
			Emin=BIGNUMBER;
			//int idmin=-1;
			for(j=kmin;j<=kmax;j++){
				id=INet(j,k);
				if(Net(id).E<Emin){
					//idmin=id;
					Emin=Net(id).E;
					width_zero_T(k)=j;
				}
			}
		}
	}
			
	/* (0,0) points to the starting point for polymer backward reconstruction */
	int idaux; 
	k=Ihmax-1;
	if(mode==0){
		//-- "two cones" constraint: PBC
		kmin=( (-k*Mcut) >= (kaux=Mcut*(k-(idist+1))+mend-mstart) )?(-k*Mcut):(kaux);
		kmax=( k*Mcut <= (kaux=Mcut*((idist+1)-k)+mend-mstart) )?(k*Mcut):(kaux);
	}
	if(mode==1){
		// one cone constraint: OBC
		kmin=(-k*Mcut);
		kmax=( k*Mcut);
	}

	Emin = BIGNUMBER;
	//for(int jj= -Mcut*(Ihmax-1);jj<= Mcut*(Ihmax-1);jj++)
	for(int jj= kmin;jj<= kmax;jj++)
	{
		idaux=INet(jj,k);
		if(Net(idaux).E<Emin){
			Emin=Net(idaux).E;
			Net(0).best=idaux;
		}
	}
	

	/* descends from (0,0) through the minimal energy polymer */
	//cout << "descending" << endl;			
	id=0;
	kaux=0;
	for(kk = Ihmax-1 ;kk >= 1; kk=kk-1)
	{
		idb=Net(id).best;
		id=idb; 
		//if((mb = Net(idb).j+mstart)!=Is(Ih(istart+kk)))
		#ifdef OPTIMIZE3
		mb = Net(idb).j+mstart;
		iflag=1; kaux+=(mb-Is(Ih(istart+kk))); Is(Ih(istart+kk))=mb;
		#else
		mb = Net(idb).j+mstart;
		iflag=1; Is(Ih(istart+kk))=mb;				
		#endif
		if(beta<0.0){ 
			Is_diffusion(kk) += pow2(Net(idb).j); 
			Is_energy2(k)+=pow2(Net(idb).E);
			Is_energy(k)+=(Net(idb).E);
		}
	
	}

	//cout << "change" << endl;
	#ifdef OPTIMIZE3
	if(measure==1) Ropt.acum(idist,kaux);				
	#else
	if(measure==1) Ropt.acum(idist,sum(Is-Is_replic));				
	#endif

	//print_diffusion_inst();

	return iflag;
/* 
iflag=0 --> polymer not changed by minimization
iflag=1 --> polymer changed by minimization
*/
	
}

/* Transfer matrix construction of the ground state 
of a polymer constrained in TWO points 
*/
int sistemavmc::Minimize_local(int start_avalanche, int size_avalanche, Array<int,1> &X)
{
	int iflag;
	int idist, istart, mstart,mend,mb;
	int j,k,id,l,idn,idb,kk, kmin, kmax, kaux;
	Real Emin,Ebond,Edis;		

	iflag=0;	

	idist=size_avalanche; //-- avalanche size
	istart=start_avalanche; //-- avalanche start

	//-- end points
	mstart=X(istart);
	mend=X(Ih(istart+idist));

	/* construction level by level of minimal energy polymers 
	   size polymers = idist 
	   starting point = istart 
	*/
	for(k=1;k<=idist;k++){

		//-- "two cones intersection" limits
		kmin=( (-k*Mcut) >= (kaux=Mcut*(k-idist)+mend-mstart) )?(-k*Mcut):(kaux);
		kmax=( k*Mcut <= (kaux=Mcut*(idist-k)+mend-mstart) )?(k*Mcut):(kaux);		

		for(j=kmin;j<=kmax;j++){
			id=INet(j,k); //-- node (k,j)==(z,x) 

			Emin=BIGNUMBER; 
			for(l=0;l<Net(id).n;l++){//-- conected nodes below
				idn = Net(id).list_below[l];				

				Ebond= Net(id).bond_below[l];

				if((Net(idn).E+Ebond) <Emin){ //-- polymer energy up to node <idn>
					Emin=Net(idn).E+Ebond;
					Net(id).best=idn;
				}
			}
			//-- energy of polymer up to the node <id>=(k,j)
			Net(id).E = Emin + V(Ih(istart+k),IM(mstart+j)) - F*j; 									
		}
	}

	/* descends through minimal energy polymer */
	k=idist; j=mend-mstart;	
	id=INet(j,k);kaux=0;

	idb=Net(id).best;
	if((Net(idb).j+mstart) != X(Ih(istart+idist-1)) ){									
		for(kk=idist-1;kk>=1;kk=kk-1){
			idb=Net(id).best;
			id=idb; 
			#ifdef OPTIMIZE3
			mb = Net(idb).j+mstart;
			iflag=1; kaux+=(mb-X(Ih(istart+kk))); X(Ih(istart+kk))=mb;
			#else
			mb = Net(idb).j+mstart;
			iflag=1; X(Ih(istart+kk))=mb;				
			#endif
		}

		return iflag;			
	}

	return iflag;
}



/* Transfer matrix construction of the ground state 
of a polymer constrained in TWO points with backward 
(default) and forward (controlled by flag) 
constraints Is_kmin, Is_kmax. Solution is X. If 
the solution can be obtained by closer end points 
it is rejected and the initial X is not changed.
*/
int sistemavmc::Minimize_local_constrained(int start_avalanche, int size_avalanche, Array<int,1> &X, 
					Array<int,1> &Is_kmin, Array<int,1> &Is_kmax, int flag_use_limits)
{
	int iflag;
	int idist, istart, mstart,mend,mb;
	int j,k,id,l,idn,idb,kk, kmin, kmax, kaux, kminbelow, kmaxbelow;
	Real Emin,Ebond,Edis;		

	iflag=0;	

	idist=size_avalanche; //-- avalanche size
	istart=start_avalanche; //-- avalanche start

	//-- end points
	mstart=X(istart);
	mend=X(Ih(istart+idist));

	//cout << "minimize local constrained" << endl;
	//cout << mstart << " " << mend << endl;

	/* construction level by level of minimal energy polymers 
	   size polymers = idist 
	   starting point = istart 
	*/
	for(k=1;k<=idist;k++){

		//-- "two cones intersection" limits
		kmin=( (-k*Mcut) >= (kaux=Mcut*(k-idist)+mend-mstart) )?(-k*Mcut):(kaux);
		kmax=( k*Mcut <= (kaux=Mcut*(idist-k)+mend-mstart) )?(k*Mcut):(kaux);		

		// backward constraint	
		if(kmin < Is_kmin(Ih(istart+k))-mstart) kmin=Is_kmin(Ih(istart+k))-mstart;
		
		// forward constraint (on/off by flag)	
		if(flag_use_limits==1){ 
			if(kmax > Is_kmax(Ih(istart+k))-mstart) kmax=Is_kmax(Ih(istart+k))-mstart;
		}
		
		for(j=kmin;j<=kmax;j++)
		{
			id=INet(j,k); //-- node (k,j)==(z,x) 

			kminbelow = Is_kmin(Ih(istart+k-1))-mstart;
			kmaxbelow = Is_kmax(Ih(istart+k-1))-mstart;

			Emin=BIGNUMBER; 
			for(l=0;l<Net(id).n;l++){//-- conected nodes below
				idn = Net(id).list_below[l];				
								
				if((flag_use_limits==0)?(Net(idn).j>=kminbelow):(Net(idn).j>=kminbelow && Net(idn).j<=kmaxbelow))
				{
					Ebond= Net(id).bond_below[l];

					if((Net(idn).E+Ebond) <Emin){ //-- polymer energy up to node <idn>
						Emin=Net(idn).E+Ebond;
						Net(id).best=idn;
					}
				}
			}
			//-- energy of polymer up to the node <id>=(k,j)
			Net(id).E = Emin + V(Ih(istart+k),IM(mstart+j)) - F*j; 									
		}
	}

	/* descends through minimal energy polymer */
	k=idist; j=mend-mstart;	// top point
	id=INet(j,k);kaux=0; 

	idb=Net(id).best; // continuation point below
	
	/*- put solution in X -*/
	// if the solution can be obtained for closer end points
	// the solution is not accepted. 
	if((flag_use_limits==0)?( (Net(idb).j+mstart) != X(Ih(istart+idist-1)) ):(1))
	{										
		for(kk=idist-1;kk>=1;kk=kk-1){
			idb=Net(id).best;
			id=idb; 
			#ifdef OPTIMIZE3
			mb = Net(idb).j+mstart;
			iflag=1; kaux+=(mb-X(Ih(istart+kk))); X(Ih(istart+kk))=mb;
			#else
			mb = Net(idb).j+mstart;
			iflag=1; X(Ih(istart+kk))=mb;				
			#endif
		}
		return 1; // solution 			
	}
		
	return 0; // no solution
}
