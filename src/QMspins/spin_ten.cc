/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-25-01
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
/*
	spin_ten.cc--> spacial spherical tensors up to rank 2
*/


#include "QMspins/spin_ten.h"

BEGIN_BL_NAMESPACE


bool checksys(SpinSys &A, int on){
	if(on>=A.spins()) return false;
	else return true;
}

void printerr(std::string mess, int on){
	BLEXCEPTION(mess+" spin number" + itost(on)+ "does not exsist")
}

rimatrix T0(SpinSys &A, int on){
 if(!checksys(A,on)){		printerr("T0", on);	}
	return A.Ie(on);
}

rmatrix T1(SpinSys &A, int on1, int m){
	if(!checksys(A,on1)){		printerr("T1", on1);	}

  if(m!=0&&m!=1&&m!=-1){
  	BLEXCEPTION(std::string("tensors cannot m's, ")+itost(m)+", like that!")

  }else{
    if(m==1){
      return -A.Imi(on1);
    }else if(m==0){
      return rmatrix(A.Iz(on1));
    }else if(m==-1){
      return A.Ip(on1);
    }
  }
  return A.F0();
}
rmatrix T11(SpinSys &A, int on)
{
	if(!checksys(A,on)){		printerr("T1", on);	}
	return A.Ip(on);
}

rdmatrix T10(SpinSys &A, int on)
{
	if(!checksys(A,on)){		printerr("T1", on);	}
	return A.Iz(on);
}

rmatrix T1m1(SpinSys &A, int on)
{
	if(!checksys(A,on)){		printerr("T1", on);	}
	return -A.Imi(on);
}


matrix T2(SpinSys &A, int on1, int on2, int m){
  	if(!checksys(A,on1)){		printerr("T2", on1);	}
	if(!checksys(A,on2)){		printerr("T2", on2);	}
	if(m!=0&&m!=1&&m!=-1&&m!=-2&&m!=2){
		BLEXCEPTION(std::string("tensors cannot m's, ")+itost(m)+", like that!")
	}else{
		if(m==-2){
		  return 1./2.*(A.Imi(on1)*A.Imi(on2));
		}else if(m==-1){
		  return 1./2.*(A.Imi(on1)*A.Iz(on2)+A.Iz(on1)*A.Imi(on2));
		}else if(m==0){
		  return 1./pow(6., 0.5)*(2*A.Iz(on1)*A.Iz(on2)-(A.Ix(on1)*A.Ix(on2)+A.Iy(on1)*A.Iy(on2)));
		}else if(m==1){
		  return -1./2.*(A.Ip(on1)*A.Iz(on2)+A.Iz(on1)*A.Ip(on2));
		}else if(m==2){
		  return 1./2.*(A.Ip(on1)*A.Ip(on2));
		}
	}
	return A.F0();
}

matrix T2(SpinSys &A, int on1, int m){ return T2(A, on1, on1, m);	}

rmatrix T22(SpinSys &A, int on1, int on2)
{
	if(!checksys(A,on1)){		printerr("T2", on1);	}
	if(!checksys(A,on2)){		printerr("T2", on2);	}
	return 1./2.*(A.Ip(on1)*A.Ip(on2));
}

rmatrix T21(SpinSys &A, int on1, int on2)
{
	if(!checksys(A,on1)){		printerr("T2", on1);	}
	if(!checksys(A,on2)){		printerr("T2", on2);	}
	return -1./2.*(A.Ip(on1)*A.Iz(on2)+A.Iz(on1)*A.Ip(on2));
}

matrix T20(SpinSys &A, int on1, int on2)
{
	if(!checksys(A,on1)){		printerr("T2", on1);	}
	if(!checksys(A,on2)){		printerr("T2", on2);	}
	return 1./pow(6., 0.5)*(2.*A.Iz(on1)*A.Iz(on2)-(A.Ix(on1)*A.Ix(on2)+A.Iy(on1)*A.Iy(on2)));
}

rmatrix T2m1(SpinSys &A, int on1, int on2)
{
	if(!checksys(A,on1)){		printerr("T2", on1);	}
	if(!checksys(A,on2)){		printerr("T2", on2);	}
	return 1./2.*(A.Imi(on1)*A.Iz(on2)+A.Iz(on1)*A.Imi(on2));
}

rmatrix T2m2(SpinSys &A, int on1, int on2)
{
	if(!checksys(A,on1)){		printerr("T2", on1);	}
	if(!checksys(A,on2)){		printerr("T2", on2);	}
	return 1./2.*(A.Imi(on1)*A.Imi(on2));
}

rmatrix T22(SpinSys &A, int on)
{
	if(!checksys(A,on)){		printerr("T2", on);	}
	return 1./2.*(A.Ip(on)*A.Ip(on));
}

rmatrix T21(SpinSys &A, int on)
{
	if(!checksys(A,on)){		printerr("T2", on);	}
	return -1./2.*(A.Ip(on)*A.Iz(on)+A.Iz(on)*A.Ip(on));
}

matrix T20(SpinSys &A, int on)
{
	if(!checksys(A,on)){		printerr("T2", on);	}
	return 1./pow(6., 0.5)*(2.*A.Iz(on)*A.Iz(on)-(A.Ix(on)*A.Ix(on)+A.Iy(on)*A.Iy(on)));
}

rmatrix T2m1(SpinSys &A, int on)
{
	if(!checksys(A,on)){		printerr("T2", on);	}
	return 1./2.*(A.Imi(on)*A.Iz(on)+A.Iz(on)*A.Imi(on));
}

rmatrix T2m2(SpinSys &A, int on)
{
	if(!checksys(A,on)){		printerr("T2", on);	}
	return 1./2.*(A.Imi(on)*A.Imi(on));
}

matrix T2_rot(SpinSys &A,  int on1, int on2, int m,double alpha, double beta, double gamma){
	if(!checksys(A,on1)){		printerr("T2_rot", on1);	}
	if(!checksys(A,on2)){		printerr("T2_rot", on2);	}

	return T2(A, on1, on2, -2)*wigner2(alpha, beta, gamma, m, -2)+
		T2(A, on1, on2, -1)*wigner2(alpha, beta, gamma, m, -1)+
		T2(A, on1, on2, 0)*wigner2(alpha, beta, gamma, m, 0)+
		T2(A, on1, on2,1)*wigner2(alpha, beta, gamma, m, 1)+
		T2(A, on1, on2, 2)*wigner2(alpha, beta, gamma, m, 2);
}

matrix T2_rot(SpinSys &A, int on1, int m, double alpha, double beta, double gamma)
{
	if(!checksys(A,on1)){		printerr("T2_rot", on1);	}

	return T2(A, on1, -2)*wigner2(alpha, beta, gamma, m, -2)+
		T2(A, on1, -1)*wigner2(alpha, beta, gamma, m, -1)+
		T2(A, on1,  0)*wigner2(alpha, beta, gamma, m, 0)+
		T2(A, on1, 1)*wigner2(alpha, beta, gamma, m, 1)+
		T2(A, on1,  2)*wigner2(alpha, beta, gamma, m, 2);
}



//spin tensors for the CSA.  the l=1 terms are set to 0
//as their spacial tensor comterparts are 0
//these are generated assuming a B field along the z-axis
//these are already quantized along the B field
//so there is no need to rotate them
//
//            [ 1 ]1/2
//     T  =  -| --|   * Iz
//      0,0   [ 3 ]
//
//     T  =    0
//      1,+/-1
//
//                                    [ 1 ]               [ 2 ]1/2
//     T      = 0        T    =   -/+ |---| *I       T   =|---|  * Iz
//      2,+/-2            2,+/-1      [ 2 ]   +/-     2,0 [ 6 ]

rdmatrix T00_CSA(SpinSys &A, int on)
{
	return 1./pow(3., 0.5)*A.Iz(on);
}

rdmatrix T10_CSA(SpinSys &A, int on1)
{
	return A.F0();
}

rdmatrix T11_CSA(SpinSys &A, int on1)
{
	return A.F0();
}

rdmatrix T1m1_CSA(SpinSys &A, int on1)
{
	return A.F0();
}



rdmatrix T20_CSA(SpinSys &A, int on)
{
	return A.Iz(on);
}

rdmatrix T22_CSA(SpinSys &A, int on)
{
	return A.F0();
}

rdmatrix T2m2_CSA(SpinSys &A, int on)
{
	return A.F0();
}

rmatrix T21_CSA(SpinSys &A, int on)
{
	return rmatrix(-1./2.*(A.Ip(on)));
}

rmatrix T2m1_CSA(SpinSys &A, int on)
{
	return rmatrix(1./2.*(A.Imi(on)));
}


rmatrix T_CSA(SpinSys &A, int on, int l, int m){
	if(!checksys(A,on)){		printerr("T_CSA", on);	}
	if(l==0&&m==0){
		return 1./pow(3., 0.5)*A.Iz(on);
	}else if(l==1){
		return rmatrix(A.F0());
	}else if(l==2){
		if(m==2||m==-2){
			return rmatrix(A.F0());
		}else if(m==-1){
			return 1./2.*(A.Imi(on));
		}else if(m==0){
			return A.Iz(on);
		}else if(m==1){
			return -1./2.*(A.Ip(on));
		}
	}
	return A.F0();
}

//here we have the 'perturbted' (high field approx) for the quadropolar
//tensors, both 1st and 2nd order perturbations
//NOTE: the extra 'q' is to avoid a function match with
//one in the gamma library

//
matrix T_Qq(SpinSys &A, int on1,int order, int m){
  if(!checksys(A,on1)){		printerr("T_Qq", on1);	}
/* if(order==1&&(m!=0&&m!=1&&m!=-1&&m!=-2&&m!=2)){
    std::cout<<"*******************"<<std::endl;
    std::cout<<"Your program is screwed"<<std::endl;
    std::cout<<"Quadropolar Spin tensors cannot"<<std::endl;
    std::cout<<" have m's like that!"<<std::endl;
    std::cout<<"*******************"<<std::endl;
        BLEXCEPTION(__FILE__,__LINE__)
  }else if(order==2&&(m!=0&&m!=1&&m!=2&&m!=3&&m!=4)){
    std::cout<<"*******************"<<std::endl;
    std::cout<<"Your program is screwed"<<std::endl;
    std::cout<<"Quadropolar Spin tensors cannot"<<std::endl;
    std::cout<<"Only l=0,1,2,3,4 are allowed (only l=0,2,4 have values!=0)"<<std::endl;
    std::cout<<" have m's like that!"<<std::endl;
    std::cout<<"*******************"<<std::endl;
        BLEXCEPTION(__FILE__,__LINE__)
  }else if(order!=1&&order!=2){
    std::cout<<"*******************"<<std::endl;
    std::cout<<"Your program is screwed"<<std::endl;
    std::cout<<"Our hamiltonian only goes to "<<std::endl;
    std::cout<<"2nd order in perturbation"<<std::endl;
    std::cout<<"*******************"<<std::endl;
        BLEXCEPTION(__FILE__,__LINE__)
  }
  */
	if(order==1){
		if(m==0){
			return (3.*A.Iz(on1)*A.Iz(on1));
		}else{
			return A.F0();
		}

	}else if(order==2){
		if(m==0){
			return A.qn(on1)*A.qn(on1)*A.Ie(on1)-3.*A.Iz(on1)*A.Iz(on1);
		}else if(m==11){
			return  T2m1(A,on1)*T21(A, on1)-T21(A,on1)*T2m1(A,on1);
			//Iz*(4*I*I - 8Iz*Iz - 1)
			//return A.Iz(on1)*(
			//	 4.0*(A.Iz(on1)*A.Iz(on1)+A.Ix(on1)*A.Ix(on1)+A.Iy(on1)*A.Iy(on1))-
			//	 8.0*(A.Iz(on1)*A.Iz(on1))-
			//	 A.Ie(on1));
		}else if(m==22){

			return T2m2(A,on1)*T22(A, on1)-T22(A,on1)*T2m2(A,on1);
			//Iz*(2*I - 2Iz - 1)
			//return A.Iz(on1)*(
			//	 2.0*(A.Iz(on1)*A.Iz(on1)+A.Ix(on1)*A.Ix(on1)+A.Iy(on1)*A.Iy(on1))-
			//	 2.0*(A.Iz(on1)*A.Iz(on1))-
			//	 A.Ie(on1));
		}else{
			return A.F0();
		}
	}

  /*
  //the 'm' here is acctually 'l' not m, these are in the combined space of
  //A2m,T2m-->F0+F2+F4 these are just the spin tensor parts
  //for the second order quadrupolar hamiltonian in the combined space
	}else if(order==2){
    if(m==0){
      return 1./2.*A.Iz(on1)*(A.qn(on1)*A.qn(on1)*A.Ie(on1)-3.*A.Iz(on1)*A.Iz(on1));
    }else if(m==2){
    	return 1./2.*A.Iz(on1)*(8.*A.qn(on1)*A.qn(on1)*A.Ie(on1)-12.*A.Iz(on1)*A.Iz(on1)-3.*A.Ie(on1));
    }else if(m==4){
      return 1./2.*A.Iz(on1)*(18.*A.qn(on1)*A.qn(on1)*A.Ie(on1)-34.*A.Iz(on1)*A.Iz(on1)-5.*A.Ie(on1));
    }else{
      return A.F0();
    }*/
  return A.F0();
}

matrix T_Quad(SpinSys &A, int on1, int m,int order)
{	return T_Qq(A, on1, order,m); }

rdmatrix T_QuadTotal(SpinSys &A, int on, int rank)
{
	switch(rank){
		case 2:
			return A.Iz(on)*2.0*(-8.0*A.Iz(on)*A.Iz(on)-1.0+4.0*A.qn(on)*(A.qn(on)+1));
		case 4:
			return A.Iz(on)*2.0*(-2.0*A.Iz(on)*A.Iz(on)-1.0+2.0*A.qn(on)*(A.qn(on)+1));
		default:
			BLEXCEPTION(" Only rank=2 or =4 are alowed ")
			return A.Ie(on);
	}
}


rdmatrix T00_Dip(SpinSys &A, int on1, int on2, int homo)
{
	return A.F0();
}

rdmatrix T10_Dip(SpinSys &A, int on1, int on2, int homo)
{
	return A.F0();
}
rdmatrix T11_Dip(SpinSys &A, int on1, int on2, int homo)
{
	return A.F0();
}

rdmatrix T1m1_Dip(SpinSys &A, int on1, int on2, int homo)
{
	return A.F0();
}


matrix T20_Dip(SpinSys &A, int on1, int on2, int homo)
{
	if(homo==0){
		return pow(1./6., 0.5)*(2.*A.Iz(on1)*A.Iz(on2)-A.Ix( on1)*A.Ix(on2)-A.Iy(on1)*A.Iy(on2));
	}else if(homo==1){
		return pow(2./3., 0.5)*(A.Iz(on1)*A.Iz(on2));
	}
	return A.F0();
}

rdmatrix T21_Dip(SpinSys &A, int on1, int on2, int homo)
{
	return A.F0();
}

rdmatrix T22_Dip(SpinSys &A, int on1, int on2, int homo)
{
	return A.F0();
}

rdmatrix T2m1_Dip(SpinSys &A, int on1, int on2, int homo)
{
	return A.F0();
}

rdmatrix T2m2_Dip(SpinSys &A, int on1, int on2, int homo)
{
	return A.F0();
}




//here are the 'perturbed' Dipole tensors, only first order is concidered
//for both hetero and homonuclear  couplings
//it returns a zero matrix if  for the system if a 'wrong' m value
hmatrix T_Dip(SpinSys &A, int on1, int on2,  int m){
	if(!checksys(A,on1)){		printerr("T_Dip", on1);	}
	if(!checksys(A,on2)){		printerr("T_Dip", on2);	}
	if(m!=0&&m!=1&&m!=-1&&m!=-2&&m!=2) {
		BLEXCEPTION(std::string("tensors cannot m's, ")+itost(m)+", like that!")
	}
	bool homo=A(on1)==A(on2);
	if(homo){
		if(m==0){
		  return pow(1./6., 0.5)*(2.*A.Iz(on1)*A.Iz(on2)-A.Ix( on1)*A.Ix(on2)-A.Iy(on1)*A.Iy(on2));
		}else{
		  return A.F0();
		}
	}else{
		if(m==0){
		  return pow(4./6., 0.5)*(A.Iz(on1)*A.Iz(on2));
		}else{
		  return A.F0();
		}
	}
	return A.F0();
}


hmatrix T_J(SpinSys &A, int on1, int on2, int m, int strong)
{
	if(!checksys(A,on1)){		printerr("T_J", on1);	}
	if(!checksys(A,on2)){		printerr("T_J", on2);	}
	if(m==0){

		if(!strong)	return 1./pow(3., 0.5)*A.Iz(on1)*A.Iz(on2);
		else return 1./pow(3., 0.5)*(A.Iz(on1)*A.Iz(on2)+A.Ix(on1)*A.Ix(on2)+A.Iy(on1)*A.Iy(on2));
	}
	return A.F0();
}


END_BL_NAMESPACE

