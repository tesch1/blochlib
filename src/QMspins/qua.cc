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
	qua.cc--> methods spin interaction quadrupole
*/
#include "QMspins/qua.h"
#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE


void Qua::reset(){
	Q_=0.;
	eta_=0.;
	si_=0.;
	chi_=0.;
	psi_=0.;
	on_=0;
	order=0;
	crystalAs.resize(5, ZeroType<complex>::zero());
}

Qua::Qua(){	reset();	}



void Qua::setParam(const std::string param, double value){
	if(param=="Q" || param=="del") {	Q_=value; setCrystalAs();}
	if(param=="eta") { eta_=value; setCrystalAs();}
	if(param=="si") {si_=value; setCrystalAs();}
	if(param=="chi") {chi_=value; setCrystalAs();}
	if(param=="psi") {psi_=value; setCrystalAs();}
	if(param=="on" || param=="on1") {on_=int(value);}
	if(param=="Bfield") {	Bfield_=value;	}
	if(param=="order") {	order=int(value);	}
}

double Qua::getParam(const std::string param){
	if(param=="iso") {	return 0.;}
	if(param=="del" || param=="Q") { return Q_;}
	if(param=="eta") { return eta_;}
	if(param=="si" || param=="alpha") { return si_;}
	if(param=="chi" || param=="beta") { return chi_;}
	if(param=="psi" || param=="gamma") { return psi_;}
	if(param=="chain" || param=="on") {return double(on_);}
	if(param=="Bfield") {	return Bfield_;	}
	if(param=="order") {	return double(order);	}
	return 0.;
}

void Qua::setSpinMats(SpinSys &A)
{
	T20=T_Quad(A,  on_, 0,1);
	T22T2m2=T_Quad(A,  on_, 22,2);
	T21T2m1=T_Quad(A,  on_, 11,2);
	c2=T_QuadTotal(A, on_, 2);
	c4=T_QuadTotal(A, on_,4);
}

//returns a std::string like Q(spin)
std::string Qua::name()
{	return "Q("+itost(on_)+")";	}


/*set two Qua equal to each other
 */

Qua Qua::operator=(const Qua &another){
	if(this==&another) return *this;
	Q_=another.Q_;
	eta_=another.eta_;
	si_=another.si_;
	chi_=another.chi_;
	psi_=another.psi_;
	on_=another.on_;
	crystalAs=another.crystalAs;
	order=another.order;
	T20=another.T20;
	T22T2m2=another.T22T2m2;
	T21T2m1=another.T21T2m1;
	return *this;
}

bool Qua::operator==(const Qua &another)
{
	if(Q_!=another.Q_) return false;
	if(eta_!=another.eta_)return false;;
	if(si_!=another.si_)return false;;
	if(chi_!=another.chi_)return false;;
	if(psi_!=another.psi_)return false;;
	if(on_!=another.on_)return false;;
	if(Bfield_!=another.Bfield_)return false;;
	return true;
}

bool Qua::operator!=(const Qua &another)
{
	if(Q_==another.Q_) return false;
	if(eta_==another.eta_)return false;;
	if(si_==another.si_)return false;;
	if(chi_==another.chi_)return false;;
	if(psi_==another.psi_)return false;;
	if(on_==another.on_)return false;;
	if(Bfield_==another.Bfield_)return false;;
	return true;
}



std::ostream& operator<<(std::ostream &otr,const Qua &qua){
	otr<<"-------------Quad--------------"<<std::endl;
  otr<<"spin\tQ\teta\talpha\tbeta\tgamma";
  if(qua.Bfield_>0){otr<<"\tBo(MHz)"<<std::endl;}
  else{	otr<<std::endl; }

	otr<<qua.on_<<"\t"<<qua.Q_<<"\t"<<qua.eta_<<"\t"<<
			qua.si_<<"\t"<<qua.chi_<<"\t"<<qua.psi_;
	if(qua.Bfield_>0){ otr<<"\t"<<qua.Bfield_<<std::endl;}
  	else{	otr<<std::endl; }
	return otr;
}

void Qua::display(){
  std::cout<<*this;
}

void Qua::write(std::ostream &otr){	//writes a line corresponding to the input
	otr<<"Q "<<Q_<<" "<<eta_<<" "<<on_<<" "<<si_<<" "<<chi_<<" "<<psi_<<" "<<order<<std::endl;
}


/*double Qua::freq(double a, double b, double g, double p, double t, double s){
  complex w(0.,0.);

  for(int m=-2;m<=2;m++){
		w+=wigner2(a,b,g, 0, m)*wigner2(p,t,s, m, 0)*A_Q_pas(eta, 0);
		w+=wigner2(a,b,g, 0, m)*wigner2(p,t,s, m, 2)*A_Q_pas(eta, 2);
		w+=wigner2(a,b,g, 0, m)*wigner2(p,t,s, m, -2)*A_Q_pas(eta, -2);
	}
	return pow(3./2., 0.5)*Q*Re(w);
}*/


dmatrix Qua::get_H(SpinSys &A,double alpha, double beta,double theta, double phi, double gam){
	Rotations rots(phi, theta, gam, alpha, beta);
	return get_H(A, rots);
}


void Qua::setCrystalAs(){
	if(si_==0&&chi_==0&&psi_==0){
		crystalAs[0]=A_Q_pas(eta_,2); //A-2
		crystalAs[1]=complex(0.,0.); //A-1
		crystalAs[2]=A_Q_pas(eta_,0); //A0
		crystalAs[3]=complex(0.,0.); //A1
		crystalAs[4]=A_Q_pas(eta_,-2); //A2
	}else{
		double prad=si_*DEG2RAD;
		double thrad=chi_*DEG2RAD;
		double grad=psi_*DEG2RAD;

//		double Co=1./8.*sqrt(3.)/sqrt(2.)*(3.*cos(thetap)*cos(thetap)-1.+etap*sin(thetap)*sin(thetap)*cos(2.*phip));
//		double C1=1./4.*sqrt(3./2.)*sin(2.*thetap)*(1.-1./3.*etap*cos(2.*phip));
//		double C2=1./8.*sqrt(3./2.)*(sin(thetap)*sin(thetap)+etap/3.*(cos(thetap)*cos(thetap)+1.)*cos(2.*phip));
//		double S1=1./2./sqrt(6.)*etap*sin(thetap)*sin(2.*phip);
//		double S2=1./4./sqrt(6.)*etap*cos(thetap)*sin(2.*phip);

		//A2-2
		crystalAs[0]=wigner2(prad,thrad,grad, -2,2)*A_Q_pas(eta_,2)+
					wigner2(prad,thrad,grad, -2,0)*A_Q_pas(eta_,0)+
					wigner2(prad,thrad,grad, -2,-2)*A_Q_pas(eta_,-2);
		//A2-1
		crystalAs[1]=wigner2(prad,thrad,grad, -1,2)*A_Q_pas(eta_,2)+
					wigner2(prad,thrad,grad, -1,0)*A_Q_pas(eta_,0)+
					wigner2(prad,thrad,grad, -1,-2)*A_Q_pas(eta_,-2);
		//A20
		crystalAs[2]=wigner2(prad,thrad,grad, 0,2)*A_Q_pas(eta_,2)+
					wigner2(prad,thrad,grad, 0,0)*A_Q_pas(eta_,0)+
					wigner2(prad,thrad,grad, 0,-2)*A_Q_pas(eta_,-2);
		//A21
		crystalAs[3]=wigner2(prad,thrad,grad, 1,2)*A_Q_pas(eta_,2)+
					wigner2(prad,thrad,grad, 1,0)*A_Q_pas(eta_,0)+
					wigner2(prad,thrad,grad, 1,-2)*A_Q_pas(eta_,-2);
		//A22
		crystalAs[4]=wigner2(prad,thrad,grad, 2,2)*A_Q_pas(eta_,2)+
					wigner2(prad,thrad,grad, 2,0)*A_Q_pas(eta_,0)+
					wigner2(prad,thrad,grad, 2,-2)*A_Q_pas(eta_,-2);
	}
	rmatrix rot=rotationMatrix3D(si_*DEG2RAD, chi_*DEG2RAD, psi_*DEG2RAD);
	PAS_.resize(3,3);
	PAS_.fill(0.0);
	PAS_(0,0)=(eta_-1.0)/2.0;
	PAS_(1,1)= -(eta_+1.0)/2.0;
	PAS_(2,2)=1.0;
	PAS_=rot*
	     PAS_*
	     transpose(rot);
}

dmatrix Qua::get_H(SpinSys &A, Rotations &rr){
	dmatrix temp=A.F0();
//	complex rot(0.,0.),tmp(0.,0.);
//	int m,mp;

	double spin=A[on_].qn();
	double Qp_=Q_/(4.0*spin*(2.0*spin-1.0));
	//static rmatrix secO=PAS_;
	static rmatrix rotm(3,3,0.0);

	if( rr.beta==0.0){
		if(eta_==0.0 && Qp_==0.0){	return temp; }
		else {
			rotm=rr.cartPowder*PAS_*transpose(rr.cartPowder);
			if(order==1 || abs(Bfield_)<=1e-12){
				return Qp_*rotm(2,2)*T20;
			}else if(abs(Bfield_)>0.0 && order==2){
				return Qp_*Qp_/(Bfield_*A[on_].relativeFrequency())*
						(c4*(pow(rotm(0,0)-rotm(1,1), 2.0)/4.0+rotm(0,1)*rotm(0,1))
						-c2*(rotm(0,2)*rotm(0,2)+rotm(1,2)*rotm(1,2)));
				//temp+= Qp_*Qp_/(Bfield_*A[on_].relativeFrequency())*
				//	(c4*(pow((rotm(0,0)-rotm(1,1)), 2.0)/4.0+ rotm(0,1)*rotm(0,1))
				//	-c2*(rotm(0,2)*rotm(0,2)+rotm(1,2)*rotm(1,2)));
			}else if(abs(Bfield_)>0.0){
				return Qp_*rotm(2,2)*T20+(Qp_*Qp_/(Bfield_*A[on_].relativeFrequency())*
						(c4*(pow(rotm(0,0)-rotm(1,1), 2.0)/4.0+rotm(0,1)*rotm(0,1))
						-c2*(rotm(0,2)*rotm(0,2)+rotm(1,2)*rotm(1,2))));
			}
		}
	}else {
		rotm=rr.cartSpin*rr.cartPowder;
		rotm=rotm*PAS_*transpose(rotm);
		if(order==1 || abs(Bfield_)<=1e-12){
			return Qp_*rotm(2,2)*T20;
		}else if(abs(Bfield_)>0.0 && order==2){
			return Qp_*Qp_/(Bfield_*A[on_].relativeFrequency())*
					(c4*(pow(rotm(0,0)-rotm(1,1), 2.0)/4.0+rotm(0,1)*rotm(0,1))
					-c2*(rotm(0,2)*rotm(0,2)+rotm(1,2)*rotm(1,2)));
			//temp+= Qp_*Qp_/(Bfield_*A[on_].relativeFrequency())*
			//	(c4*(pow((rotm(0,0)-rotm(1,1)), 2.0)/4.0+ rotm(0,1)*rotm(0,1))
			//	-c2*(rotm(0,2)*rotm(0,2)+rotm(1,2)*rotm(1,2)));
		}else if(abs(Bfield_)>0.0){
			return Qp_*rotm(2,2)*T20+(Qp_*Qp_/(Bfield_*A[on_].relativeFrequency())*
					(c4*(pow(rotm(0,0)-rotm(1,1), 2.0)/4.0+rotm(0,1)*rotm(0,1))
					-c2*(rotm(0,2)*rotm(0,2)+rotm(1,2)*rotm(1,2))));
		}
	}
	return temp;
}

/*** First ORder Hamiltonians ***/
dmatrix Qua::get_H1(SpinSys &A, Rotations &rr){
	dmatrix temp=A.F0();
	complex rot(0.,0.),tmp(0.,0.);
	int m,mp;

	if(rr.alpha==0.0 || rr.beta==0.0){
		complex a20;
		if(eta_==0.0 && Q_==0.0){	return temp; }
		else if(eta_==0.0 && Q_!=0.0){
				a20=rr.powderAs(2,0)*crystalAs[0]
					+rr.powderAs(2,1)*crystalAs[1]
					+rr.powderAs(2,2)*crystalAs[2]
					+rr.powderAs(2,3)*crystalAs[3]
					+rr.powderAs(2,4)*crystalAs[4];
				temp+=Q_*a20*T20;
			return temp;
		}else if(Q_!=0.0){
			complex a20(0,0);
			a20=rr.powderAs(2,2)*crystalAs[0]
				+rr.powderAs(2,1)*crystalAs[1]
				+rr.powderAs(2,2)*crystalAs[2]
				+rr.powderAs(2,3)*crystalAs[3]
				+rr.powderAs(2,4)*crystalAs[4];
				temp+=Q_*a20*T20;
			return temp;
		}
	}else {
		for(m=0;m<=4;m++){
			tmp=complex(0.,0.);
			for(mp=0;mp<=4;mp++){
				tmp+=rr.powderAs(mp,m)*crystalAs[mp];
			}
			tmp *=rr.spinnerAs[m];

			rot+=tmp;
		}
			//temp+=rot*T20;
		return Q_*rot*T20;
	}
	return Q_*temp;
}

/*** second ORder Hamil ***/
dmatrix Qua::get_H2(SpinSys &A, Rotations &rr){
	dmatrix temp=A.F0();
	//complex rot(0.,0.),tmp(0.,0.);
	//int m,mp;

	if(rr.alpha==0.0 || rr.beta==0.0){
		complex a21,a22,a2m1,a2m2; //a20
		if(eta_==0.0 && Q_==0.0){	return temp; }
		else if(eta_==0.0&&Q_!=0.0){
			if(Bfield_>0.0){
				a21=rr.powderAs(3,0)*crystalAs[0]
						+rr.powderAs(3,1)*crystalAs[1]
						+rr.powderAs(3,2)*crystalAs[2]
						+rr.powderAs(3,3)*crystalAs[3]
						+rr.powderAs(3,4)*crystalAs[4];
				a2m1=rr.powderAs(1,0)*crystalAs[0]
						+rr.powderAs(1,1)*crystalAs[1]
						+rr.powderAs(1,2)*crystalAs[2]
						+rr.powderAs(1,3)*crystalAs[3]
						+rr.powderAs(1,4)*crystalAs[4];
				a22=rr.powderAs(4,0)*crystalAs[0]
						+rr.powderAs(4,1)*crystalAs[1]
						+rr.powderAs(4,2)*crystalAs[2]
						+rr.powderAs(4,3)*crystalAs[3]
						+rr.powderAs(4,4)*crystalAs[4];
				a2m2=rr.powderAs(0,0)*crystalAs[0]
						+rr.powderAs(0,1)*crystalAs[1]
						+rr.powderAs(0,2)*crystalAs[2]
						+rr.powderAs(0,3)*crystalAs[3]
						+rr.powderAs(0,4)*crystalAs[4];
				temp-=Q_*Q_/Bfield_*(a22*a2m2*T22T2m2+a21*a2m1*T21T2m1);
			}
			return temp;
		}else if(Q_!=0.0){
			//complex a20(0,0);
			if(Bfield_>0.0){
				a21=rr.powderAs(3,0)*crystalAs[0]
						+rr.powderAs(3,1)*crystalAs[1]
						+rr.powderAs(3,2)*crystalAs[2]
						+rr.powderAs(3,3)*crystalAs[3]
						+rr.powderAs(3,4)*crystalAs[4];
				a2m1=rr.powderAs(1,0)*crystalAs[0]
						+rr.powderAs(1,1)*crystalAs[1]
						+rr.powderAs(1,2)*crystalAs[2]
						+rr.powderAs(1,3)*crystalAs[3]
						+rr.powderAs(1,4)*crystalAs[4];
				a22=rr.powderAs(4,0)*crystalAs[0]
						+rr.powderAs(4,1)*crystalAs[1]
						+rr.powderAs(4,2)*crystalAs[2]
						+rr.powderAs(4,3)*crystalAs[3]
						+rr.powderAs(4,4)*crystalAs[4];
				a2m2=rr.powderAs(0,0)*crystalAs[0]
						+rr.powderAs(0,1)*crystalAs[1]
						+rr.powderAs(0,2)*crystalAs[2]
						+rr.powderAs(0,3)*crystalAs[3]
						+rr.powderAs(0,4)*crystalAs[4];

				temp-=Q_*Q_/Bfield_*(a22*a2m2*T22T2m2+a21*a2m1*T21T2m1);
			}
			return temp;
		}
	}else {
		if(Bfield_>0.0){
			Vector<complex> as(5,0);
			rr.set_spinnerAsMat(rr.alpha, rr.beta);

			as=rr.powderAs*crystalAs;
			as=rr.spinnerAsMat*as;
			if(order!=2){
				return Q_*as[2]*T20
				   -Q_*(Q_/Bfield_)*(
					   as[1]*as[3]*T21T2m1+
					   as[0]*as[4]*T22T2m2
					);
			}else{
				return  -Q_*(Q_/Bfield_)*(
					   as[1]*as[3]*T21T2m1+
					   as[0]*as[4]*T22T2m2
					);
			}

		}
	}
	return Q_*temp;
}


/*
hmatrix Qua::get_H(SpinSys &A,double alpha, double beta, double gam,double theta, double phi, double B_field){
	hmatrix temp=A.F0();                 //first order perturbation piece
	hmatrix temp2=A.F0();
  	static double rot;
  	static complex koo;
 //ignore second order for these cases...it will be added below
  	if(si_==0.0 && chi_==0.0 && psi_==0.0 && alpha==0.0 && beta==0.0 && eta_==0.){  //no eta_ and not other rotations (i.e. no spin)
    	rot=Q_*(pow(1.5, 0.5)*(1. + 3.*cos(2.*theta)))/4.;
    	temp+=T20*rot;
  	}else if(si_==0.0 && chi_==0.0 && psi_==0.0 && alpha==0.0 && beta==0.0){  //w/ eta_ and no spinning
    	rot=Q_*(pow(1.5, 0.5)*(1. + 3.*cos(2.*theta) + 2.*eta_*cos(2*phi)*sin(theta)*sin(theta)))/4.;
    	temp+=T20*rot;
    }else if(si_==0.0 && chi_==0.0 && psi_==0.0 && eta_==0.){  //no eta_, but spinning
    	rot=Q_*(pow(1.5, 0.5)*(1. + 3.*cos(2.*theta) + 3.*cos(2.*beta)*(1 + 3.*cos(2.*theta)) +
       12.*cos(2.*(alpha))*sin(beta)*sin(beta)*sin(theta)*sin(theta) +
       12.*cos(alpha)*sin(2.*beta)*sin(2.*theta)))/16.;
    	temp+=Q_*rot*T20;
    }else if(si_==0.0 && chi_==0.0 && psi_==0.0){ //spinning with eta_ (messy, but as simple as it gets
    	koo=Q_*0.688919*(1./9. + 1./3.*cos(2.*theta) +
     4./9.*eta_*cos(4.*phi)*cos(2.*alpha - 2.*gam + 2.*phi)*pow(cos(0.5*theta),4.)*
      pow(sin(beta),double(2.)) + 4./9.*eta_*cos(2.*alpha - 2.*(gam + phi))*pow(cos(0.5*theta),double(4.))*
      pow(sin(beta),double(2.)) + complex(0.,4./9.)*eta_*cos(2.*alpha - 2.*gam + 2.*phi)*
      pow(cos(0.5*theta),double(4.))*pow(sin(beta),double(2.))*sin(4.*phi) -
     complex(0.,4./9.)*eta_*cos(4.*phi)*pow(cos(0.5*theta),4.)*pow(sin(beta),double(2.))*
      sin(2.*alpha - 2.*gam + 2.*phi) + 4./9.*eta_*pow(cos(0.5*theta),double(4.))*pow(sin(beta),double(2.))*
      sin(4.*phi)*sin(2.*alpha - 2.*gam + 2.*phi) +
     complex(0.,4./9.)*eta_*pow(cos(0.5*theta),double(4.))*pow(sin(beta),double(2.))*
      sin(2.*alpha - 2.*(gam + phi)) + 4./9.*eta_*cos(2.*alpha - 2.*gam + 2.*phi)*
      pow(sin(beta),double(2.))*pow(sin(0.5*theta),double(4.)) +
     4./9.*eta_*cos(4.*phi)*cos(2.*alpha - 2.*(gam + phi))*pow(sin(beta),double(2.))*
      pow(sin(0.5*theta),double(4.)) + complex(0.,4./9.)*eta_*cos(2.*alpha - 2.*(gam + phi))*
      pow(sin(beta),double(2.))*sin(4.*phi)*pow(sin(0.5*theta),double(4.)) -
     complex(0.,4./9.)*eta_*pow(sin(beta),double(2.))*sin(2.*alpha - 2.*gam + 2.*phi)*
      pow(sin(0.5*theta),double(4.)) + complex(0.,4./9.)*eta_*cos(4.*phi)*pow(sin(beta),double(2.))*
      sin(2.*alpha - 2.*(gam + phi))*pow(sin(0.5*theta),double(4.)) -
     4./9.*eta_*pow(sin(beta),double(2.))*sin(4.*phi)*sin(2.*alpha - 2.*(gam + phi))*
      pow(sin(0.5*theta),double(4.)) + 16./3.*cos(beta)*cos(alpha - gam)*cos(theta)*sin(beta)*
      sin(theta) - 16./9.*eta_*cos(beta)*cos(alpha - gam)*cos(2.*phi)*cos(theta)*sin(beta)*
      sin(theta) - 16./9.*eta_*cos(beta)*sin(beta)*sin(alpha - gam)*sin(2.*phi)*sin(theta) +
     2./9.*eta_*cos(2.*phi)*pow(sin(theta),double(2.)) +
     2./3.*cos(2.*phi)*cos(2.*alpha - 2.*gam + 2.*phi)*pow(sin(beta),double(2.))*pow(sin(theta),double(2.)) +
     2./3.*cos(2.*phi)*cos(2.*alpha - 2.*(gam + phi))*pow(sin(beta),double(2.))*pow(sin(theta),double(2.)) +
     complex(0.,2./3.)*cos(2.*alpha - 2.*gam + 2.*phi)*pow(sin(beta),double(2.))*sin(2.*phi)*
      pow(sin(theta),double(2.)) + complex(0.,2./3.)*cos(2.*alpha - 2.*(gam + phi))*pow(sin(beta),double(2.))*
      sin(2.*phi)*pow(sin(theta),double(2.)) - complex(0.,2./3.)*cos(2.*phi)*pow(sin(beta),double(2.))*
      sin(2.*alpha - 2.*gam + 2.*phi)*pow(sin(theta),double(2.)) +
     2./3.*pow(sin(beta),double(2.))*sin(2.*phi)*sin(2.*alpha - 2.*gam + 2.*phi)*pow(sin(theta),double(2.)) +
     complex(0.,2./3.)*cos(2.*phi)*pow(sin(beta),double(2.))*sin(2.*alpha - 2.*(gam + phi))*
      pow(sin(theta),double(2.)) - 2./3.*pow(sin(beta),double(2.))*sin(2.*phi)*sin(2.*alpha - 2.*(gam + phi))*
      pow(sin(theta),double(2.)) + cos(2.*beta)*(1./3. + cos(2.*theta) +
        2./3.*eta_*cos(2.*phi)*pow(sin(theta),double(2.))));
    	temp+=T20*koo;
//FIX ME
  }else{}

  return (temp);

}*/

dmatrix Qua::get_HINT(SpinSys &A, double wr, double alpha1,double alpha2,
		                double beta, double theta, double phi,double B_field){
	dmatrix temp=A.F0();
  return temp;
}

/*hthis reads in the various peices of needed for a quadrpole
 *from any old file
 *for quadropoles, Q and eta and the spin label
	also we need the spin qunatum number and orientation
with respect to the PAS.  if no angles are given
values of 0,0,0 are assumed
if no spin quanutm number is given, m=1 is assumed
the rest must be put down
eg
coupling type	Q	  eta	 on	 alpha 	beta	gamma
	Q	          12	43	 0	 	60		40		3
*/

Vector<Qua> read_qua(std::string in){
	std::ifstream infile(in.c_str());
	if(infile.fail()){
		BLEXCEPTION(std::string("QUA:: no spin file name like ")+in)
	}
 	return read_qua(infile);
}

Vector<Qua> read_qua(const char * in){
	std::ifstream infile(in);
	if(infile.fail()){
		BLEXCEPTION(std::string("QUA:: no spin file name like ")+in)
	}
 	return read_qua(infile);
}

Vector<Qua> read_qua(std::ifstream & infile){
	char liner[1000];
	Vector<std::string> ii;
	int temp=0;
	if(infile.fail()){
		BLEXCEPTION("given Spin system file cannot be opened ")
  	}

	while((temp=infile.peek())!=EOF){
    	infile.getline(liner,1000,'\n');
    	ii.push_back(std::string(liner));
	}
	infile.close();
	return read_qua(ii);
}

Vector<Qua> read_qua(const Vector<std::string> &ss){
	Vector<Qua> qs;
	int qs_adv=-1;			//counter for num csa
	Vector<std::string> tmm;
	std::string tmS;
	int ll=0, i, len=ss.size();
	bool got=false;
	for(i=0;i<len;i++){
		got=false;
		tmm=parse_param(ss[i]);
		ll=tmm.size();
		if(!tmm.empty()){
			if(ll>=1){
				if(tmm[0]=="Q" || tmm[0]=="QUAD" || tmm[0]=="Quad"){
					if(ll>=4){
						qs_adv++;
						qs.push_back(Qua());
						tmS=removeWhite(tmm[1]);
						if(tmS.size()>0)	qs[qs_adv].setParam("del",std::atof(tmS.c_str()));

						tmS=removeWhite(tmm[2]);
						if(tmS.size()>0)	qs[qs_adv].setParam("eta",std::atof(tmS.c_str()));

						tmS=removeWhite(tmm[3]);
						if(tmS.size()>0){
							qs[qs_adv].setParam("on",std::atoi(tmS.c_str()));
							got=true;
						}else{
							got=false;
						}
					}
					if(ll>=5){
						tmS=removeWhite(tmm[4]);
						if(tmS.size()>0)
						{	qs[qs_adv].setParam("si",std::atof(tmS.c_str()));	}
					}
					if(ll>=6){
						tmS=removeWhite(tmm[5]);
						if(tmS.size()>0)
						{	qs[qs_adv].setParam("chi",std::atof(tmS.c_str()));	}
					}
					if(ll>=7){
						tmS=removeWhite(tmm[6]);
						if(tmS.size()>0)
						{	qs[qs_adv].setParam("psi",std::atof(tmS.c_str()));	}
					}
					if(ll>=8){
						tmS=removeWhite(tmm[7]);
						if(tmS.size()>0)
						{	qs[qs_adv].setParam("order",std::atof(tmS.c_str()));	}
					}
					if(!got){
						BLEXCEPTION(std::string(" A Quadrpole must have at least 3 parameters..")+
							std::string("\n Q <Q> <eta> <spin> <alpha> <beta> <gamma> <order>"))
					}
					qs[qs_adv].setCrystalAs();
				}
			}
		}
	}
  	return qs;
}

END_BL_NAMESPACE

