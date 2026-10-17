#include "vmc.hpp"

//#define DEBUG

#ifdef DEBUG
ofstream debug("debug.dat");
#endif 

Real sistemavmc::VMC(int idir)
{
	VMC(idir, 0.0, 0.0);
}

//#define EPS
// Arraus: Is, Ih, IM, V
Real sistemavmc::VMC(int idir, float w, float m2)
{
	int iflag, il, i, idelta1, istop, idelta2,j, ipoint, 
	mpoint, mfpoint, maxil, nav, ilmin, imin;

	Real delta, delta_min, daux;
	delta_min=BIGNUMBER;

	Is_replic=Is;
#ifndef LOCAL_MOVES
	maxil=Ihmax+1;
#else
	maxil=2;	
#endif	
	
	nav=0;
	for(il=2;il<=maxil;il++)
	{// size of avalanche
		int max_istart=(il<maxil)?(Ihmax):(1);
		for(i=0;i<max_istart ;i++)
		{// starting point

			// check mcut
			bool mcut_condition;
			//cout << extra << endl;
			istop=Ih(i+il);

			double extra_stop = (istop==0)?(shift*Ihmax):(0.0);
			double extra_start = (i+1==Ihmax)?(shift*Ihmax):(0.0);

			if(il<Ihmax+1){ // tests mcut at proposed avalanche borders
				mcut_condition = abs(Is(Ih(i+1)) + idir - Is(i) + extra_start)<=Mcut 
								  && abs(Is(Ih(istop-1))+idir-Is(istop)-extra_stop)<=Mcut;
			}
			else mcut_condition=1; // only for rigid global displacement 
				
			if(mcut_condition)
			{				
				//cout << "extra = " << extra << endl;		

				// elastic energy in the extremes
				if(il<Ihmax+1) delta=1.0+idir*(Is(Ih(i+1))-Is(i) + extra_start + Is(Ih(istop-1))-Is(istop) - extra_stop);
				else delta=0.0;

				delta*=Cel;

				// energy change for the avalanche: disorder + drive
				for(j=1;j<il;j++){
					ipoint=Ih(i+j);
					mpoint=IM(Is(ipoint));
					mfpoint=IM(Is(ipoint)+idir);
					#ifndef TIPDRIVENDEP
					delta=delta-V(ipoint,mpoint)+V(ipoint,mfpoint)-(idir*F)+m2*0.5*(1+2.0*idir*(Is(ipoint)-w));
					//delta=delta-V(ipoint,mpoint)+V(ipoint,mfpoint)-(idir*F);
					#else
					delta=delta-V(ipoint,mpoint)+V(ipoint,mfpoint)-(idir*F); //+0.1*0.5*(1+2.0*idir*(Is(ipoint)-w))*(ipoint==0);
					#endif	
				} //-- delta: total change of energy for the avalanche
				
				daux=delta/(il-1.0); // change of energy per unit length
				if(daux < delta_min) delta_min=daux;
				
				//-- move if avalanche decreases energy
				if (delta<0){
					//cout << "vmc avalanche dE<0 : " << il << endl;
					for(j=i+1;j<i+il;j++){
						Is(Ih(j))=Is(Ih(j))+idir; // make avalanche
					}
					// return energy of smaller avalanche decreasing Energy 
					return delta_min;
				}
				nav++;
			}
		}
	}
	//cout << "vmc avalanche dE>0 : " << il-1 << endl;
	return delta_min; 
}

/*--- grow the tree of configurations ---
given an initial configuration Is and a 
limiting configuration Is_limit 
*/
int sistemavmc::grow_Tree_new2(int idir, Real Wmax)
{	
	int iflag, il, i, idelta1, istop, idelta2,j, ipoint, 
	mpoint, mfpoint, maxil, stop_flag, inserted;

	Real delta, delta_min, daux, E_fils, E_mere;
	delta_min=BIGNUMBER;

	key_Archive x;
	x.u.resize(Ihmax-1);

#ifndef LOCAL_MOVES
	maxil=Ihmax+1;
#else
	maxil=2;	
#endif	
	
	E_mere = Energy(); // la mere est Is	
	inserted=0;

	for(il=2;il<=maxil;il++)
	{// size of avalanche
		for(i=0;(il<=Ihmax)?(i<Ihmax):(i<1) ;i++)
		{// starting point

			stop_flag=1;
			Is_replic = Is;

			istop=Ih(i+il);
			double extra_stop = (istop==0)?(shift*Ihmax):(0.0);
			double extra_start = (i+1==Ihmax)?(shift*Ihmax):(0.0);

			// check mcut	
			bool mcut_condition;	
			if(il<Ihmax+1){
				mcut_condition = 
				(abs(Is(Ih(i+1))+idir-Is(i)+extra_start)<=Mcut 
				&& abs(Is(Ih(istop-1))+idir-Is(istop)-extra_stop)<=Mcut);
			}
			else mcut_condition=1;
							
			if(mcut_condition)
			{												

				// elastic energy change around the fixed points
				if(il<Ihmax+1) delta=1.0+idir*(Is(Ih(i+1))+Is(Ih(istop-1))-Is(i)-Is(istop)+ extra_start - extra_stop);
				else delta=0.0;
				delta*=Cel; // comment if Cel==1

				for(j=1;j<il && stop_flag==1;j++){
					ipoint=Ih(i+j);
					mpoint=IM(Is(ipoint));
					Is_replic(ipoint)=Is(ipoint)+idir;
					if(Is_replic(ipoint) > Is_limit(ipoint)) // envelope constraint
					{
						stop_flag=0;
					}
					else
					{
						mfpoint=IM(Is(ipoint)+idir);
						delta=delta-V(ipoint,mpoint)+V(ipoint,mfpoint)-(idir*F);
					}
				} //-- delta: total change of energy due to the avalanche (istart,idist)
		
				if(stop_flag==1)
				{
					if(delta>Wmax){ // new configuration penalized: too energetic
						stop_flag=0;
					}			
				}
			
				// if new configuration not rejected: ordered insertion in the list
				if(stop_flag==1)
				{												
					E_fils = E_mere + delta;
					//x=create_fils_set(Is_replic,E_fils,i,il);

					create_fils_set_new(Is_replic,E_fils,i,il,x);

					int s1=Archive_set.size();

					pair<set<key_Archive,cmp_Archive>::iterator, bool> result = Archive_set.insert(x);
					//for(j=0;j<Ihmax-1;j++) cout << "en la cocina " <<  j << " = " << (result.first)->u[j] << endl;
					//cout.flush();
				}					
			}
		}
	}
	return Archive_set.size();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void sistemavmc::create_fils_set_new(Array<int,1> & Y, Real E, int istart, int idist, key_Archive &x)
{
//	key_Archive x;
	//x.u.clear();
	x.E = E;
	x.i0 = Y(0);	

#ifdef USING_BITSET_VECTOR	
//	x.u.reserve(Ihmax-1);
	for(int i=1;i<Ihmax;i++){ 
		if(Y(i)-Y(i-1)+Mcut < 0 || Y(i)-Y(i-1)+Mcut > 9)
		{
			cout << "error: Mcut is too small" << endl;
			exit(1);
		}	
		//else x.u.push_back(int_mcut(Y(i)-Y(i-1)+Mcut));		
		else x.u[i-1]=(int_mcut(Y(i)-Y(i-1)+Mcut));		
		//cout << x.u[i-1] << endl;
	}
#endif		
	//exit(1);
	x.idist = idist; 
	x.istart = istart;
//	return x;
}

//////////////////////////////////////////////////////////////////////////////////////////
key_Archive sistemavmc::create_fils_set(Array<int,1> & Y, Real E, int istart, int idist)
{
	key_Archive x;
	x.E = E;
	x.i0 = Y(0);	

#ifdef USING_CHAR_ARRAY	
	x.u = (char *)malloc(sizeof(char)*(Ihmax-1));	
	for(int i=1;i<Ihmax;i++){
		if(Y(i)-Y(i-1) < -128 || Y(i)-Y(i-1) > 128)
		{
			cout << "error: |elogantion| > 128" << endl;
			exit(1);
		}	
		else		
		x.u[i-1]=char(Y(i)-Y(i-1));
	}
#endif	

#ifdef USING_BITSET_VECTOR	
	x.u.reserve(Ihmax-1);
	for(int i=1;i<Ihmax;i++){ 
		if(Y(i)-Y(i-1)+Mcut < 0 || Y(i)-Y(i-1)+Mcut > 9)
		{
			cout << "error: Mcut is too small" << endl;
			exit(1);
		}	
		else x.u.push_back(int_mcut(Y(i)-Y(i-1)+Mcut));		
	}
#endif		
	x.idist = idist; 
	x.istart = istart;
	return x;
}

