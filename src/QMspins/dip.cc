/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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
 dip.h --> spin interactions dip
*/


#include "QMspins/dip.h"

BEGIN_BL_NAMESPACE



void Dip::reset(){
	del_=0.;
	si_=0.;
	chi_=0.;
	psi_=0.;
	on1_=0;
	on2_=0;
	crystalAs.resize(5, ZeroType<complex>::zero());
}

Dip::Dip()	{	reset();	}

//this function sets the parameters inside the class
//
void Dip::setParam(const std::string param,double value){
	if(param=="del") {	del_=value; setCrystalAs();}
	else if(param=="si"|| param=="alpha") {si_=value; setCrystalAs(); }
	else if(param=="chi" || param=="beta") {chi_=value; setCrystalAs();}
	else if(param=="psi" || param=="gamma") {psi_=value; setCrystalAs();}
	else if(param=="on1" || param=="chain1") {on1_=int(value);}
	else if(param=="on2" || param=="chain2") {on2_=int(value);}
}

double Dip::getParam(const std::string param){
	if(param=="iso") {	return 0.;}
	if(param=="del") { return del_;}
	if(param=="eta") { return 0.;}
	if(param=="si" || param=="alpha") { return si_;}
	if(param=="chi" || param=="beta") { return chi_;}
	if(param=="psi" || param=="gamma") { return psi_;}
	if(param=="chain1" || param=="on1") {return double(on1_);}
	if(param=="chain2" || param=="on2") {return double(on2_);}
	return 0.;
}


//returns a std::string like 'J(spin1, spin2)'
std::string Dip::name()
{	return "D("+itost(on1_)+","+itost(on2_)+")";	}

void Dip::setSpinMats(SpinSys &A)
{
	T20=T_Dip(A, spinON1(), spinON2(),0);
}

void Dip::display(){
	std::cout<<"-------------Dip--------------"<<std::endl;
  std::cout<<"spin\tdel\talpha\tbeta\tgamma"<<std::endl;
	std::cout<<"D"<<on1_<<","<<on2_<<"\t"<<del_<<"\t"<<
			si_<<"\t"<<chi_<<"\t"<<psi_<<"\t"<<std::endl;
}

void Dip::write(std::ostream &otr){	//writes a line corresponding to the input
	otr<<"D "<<del_<<" "<<on1_<<" "<<on2_<<" "<<" "<<si_<<" "<<chi_<<" "<<psi_<<std::endl;
}


Dip Dip::operator=(const Dip &another){
	if(this==&another) return *this;
	del_=another.del_;
	si_=another.si_;
	chi_=another.chi_;
	psi_=another.psi_;
	on1_=another.on1_;
	on2_=another.on2_;
	crystalAs=another.crystalAs;
	T20=another.T20;
	return *this;
}

bool Dip::operator==(const Dip &another){
	if(this==&another) return true;
	if(del_!=another.del_) return false;
	if(si_!=another.si_) return false;
	if(chi_!=another.chi_) return false;
	if(psi_!=another.psi_) return false;
	if(on1_!=another.on1_) return false;
	if(on2_!=another.on2_) return false;
	return true;
}

bool Dip::operator!=(const Dip &another){
	if(this==&another) return false;
	if(del_==another.del_) return false;
	if(si_==another.si_) return false;
	if(chi_==another.chi_) return false;
	if(psi_==another.psi_) return false;
	if(on1_==another.on1_) return false;
	if(on2_==another.on2_) return false;
	return true;
}



std::ostream& operator<<(std::ostream &otr,const Dip &dip){
	std::cout<<"-------------Dip--------------"<<std::endl;
  std::cout<<"spin\tdel\talpha\tbeta\tgamma"<<std::endl;
	std::cout<<"D"<<dip.on1_<<","<<dip.on2_<<"\t"<<dip.del_<<"\t"<<
			dip.si_<<"\t"<<dip.chi_<<"\t"<<dip.psi_<<std::endl;
	return otr;
}

void Dip::setCrystalAs(){
	if(si_==0&&chi_==0&&psi_==0){
		crystalAs[0]=ZeroType<complex>::zero();
		crystalAs[1]=ZeroType<complex>::zero();
		crystalAs[2]=A_Dip_pas(del_,0);
		crystalAs[3]=ZeroType<complex>::zero();
		crystalAs[4]=ZeroType<complex>::zero();
	}else{
		double prad=si_*DEG2RAD;
		double thrad=chi_*DEG2RAD;
		double grad=psi_*DEG2RAD;
		//A2-2
		crystalAs[0]=wigner2(prad,thrad,grad, -2,0)*A_Dip_pas(del_,0);

		//A2-1
		crystalAs[1]=wigner2(prad,thrad,grad, -1,0)*A_Dip_pas(del_,0);
		//A20
		crystalAs[2]=wigner2(prad,thrad,grad, 0,0)*A_Dip_pas(del_,0);
		//A21
		crystalAs[3]=wigner2(prad,thrad,grad, 1,0)*A_Dip_pas(del_,0);
		//A22
		crystalAs[4]=wigner2(prad,thrad,grad, 2,0)*A_Dip_pas(del_,0);
	}
	rmatrix rot=rotationMatrix3D(si_*DEG2RAD, chi_*DEG2RAD, psi_*DEG2RAD);
	PAS_.resize(3,3);
	PAS_.fill(0.0);
	PAS_(0,0)=-0.5*del_*pow(6.0, 0.5);
	PAS_(1,1)=-0.5*del_*pow(6.0, 0.5);
	PAS_(2,2)=1.0*del_*pow(6.0, 0.5);
	PAS_=rot*
		 PAS_*
		 transpose(rot);

}

hmatrix Dip::get_H(SpinSys &A, Rotations &rr){
	static rmatrix rotm(3,3);
	if(rr.beta==0.0){

		if(del_!=0.0){
			/*	return (rr.powderAs(2,0)*crystalAs[0]+
				rr.powderAs(2,1)*crystalAs[1]+
				rr.powderAs(2,2)*crystalAs[2]+
				rr.powderAs(2,3)*crystalAs[3]+
				rr.powderAs(2,4)*crystalAs[4])*T20;
			*/
			rotm=rr.cartPowder*PAS_*transpose(rr.cartPowder);
			return rotm(2,2)*T20;
		}
	}else {
	/*	complex rot(0.0,0.0);
		int m, mp;
		if(si_==0.0 && chi_==0.0 && psi_==0.0){
			for(m=0;m<=4;m++){
				rot+=rr.spinnerAs[m]*rr.powderAs(m,2)*crystalAs[2];
			}
		}else{
			int m,mp;
			for(m=0;m<=4;m++){
				for(mp=0;mp<=4;mp++){
					rot+=rr.spinnerAs[m]*rr.powderAs(m,mp)*crystalAs[mp];
				}
			}
		}*/
		rotm=rr.cartSpin*rr.cartPowder;
		rotm=rotm*PAS_*transpose(rotm);
		return rotm(2,2)*T20;

		//return rot*T20;
	}
	return A.F0();
}

hmatrix Dip::get_H(SpinSys &A,double alpha, double beta, double theta, double phi, double gam)
{
	Rotations rots(phi, theta, gam, alpha, beta);
	return get_H(A, rots);
}


hmatrix Dip::get_HINT(SpinSys &A, double wr, double tau1,double tau2,
		                double beta, double theta, double phi){
	hmatrix temp=A.F0();



	complex rot(0.0,0.0);
	if(wr==0.0 || (tau2-tau1)==0.0){ return get_H(A, 0,0,theta, phi); }

	else if(si_==0.0 && chi_==0.0 && psi_==0.0){
  	rot= -((pi*tau1*(wr) - pi*tau2*( wr) +
        3.*pi*(tau1 - tau2)*( wr)*cos(2.*beta) +
        3.*pi*(tau1 - tau2)*( wr)*(1. + 3.*cos(2.*beta))*
         cos(2.*theta) + 6.*sin(2.*beta)*sin(2.*theta)*
         sin(phi + 2.*pi*tau1*( wr)) +
        3.*pow(sin(beta),double(2.))*pow(sin(theta),double(2.))*
         sin(2.*(phi + 2.*pi*tau1*( wr))) -
        6.*sin(2.*beta)*sin(2.*theta)*sin(phi + 2.*pi*tau2*( wr)) -
        3.*pow(sin(beta),double(2.))*pow(sin(theta),double(2.))*
         sin(2.*(phi + 2.*pi*tau2*( wr)))))/(8.*( wr));
  	rot*=A_Dip_pas(del_, 0);
  	return T20*rot;
  }
  return temp;
}

/*given a file name, this little bit will dynamically add dipoles
 *to a list.
 *the systax for the text file is as follows.
 *
	option: homo_(0) or hetero(1)
	the oriientaion angles are optional the rest should be there
	eg
	  flag	coupling  onspin tospin	 alpha beta gamma
		 D   700	   0	   1       0   0    0

*/

Vector<Dip> read_dip(std::string in){
	std::ifstream infile(in.c_str());
	 if(infile.fail()){
		BLEXCEPTION(std::string("DIP:: no spin file name like ")+in)
  }
 	return read_dip(infile);
}

Vector<Dip> read_dip(const char * in){
	std::ifstream infile(in);
	 if(infile.fail()){
    	BLEXCEPTION(std::string("DIP:: no spin file name like ")+in)
  }
 	return read_dip(infile);
}

Vector<Dip> read_dip(std::ifstream &infile){
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
	return read_dip(ii);
}

Vector<Dip> read_dip(const Vector<std::string> &ss){
	Vector<Dip> dip;

	int dip_adv=-1;			//counter for num dip
	Vector<std::string> tmm;
	int i, len=ss.size(),ll=0;
	for(i=0;i<len;i++){
		tmm=parse_param(ss[i]);
		ll=tmm.size();
		if(!tmm.empty()){
			switch(tmm[0][0]){
			case '#': case '\0': case '\n':   //ignore lines that start with these
			  break;
			default:
				if(ll>=1){
					if(tmm[0]=="D" || tmm[0]=="Dip" || tmm[0]=="dip" || tmm[0]=="DIP"){
						if(ll>=4){
							dip_adv++;
							dip.push_back(Dip());
							dip[dip_adv].setParam("del",std::atof(tmm[1].c_str()));
							dip[dip_adv].setParam("on1",std::atof(tmm[2].c_str()));
							dip[dip_adv].setParam("on2",std::atof(tmm[3].c_str()));
						}
						if(ll>=5) if(tmm[4].size()>0) {	dip[dip_adv].setParam("si",std::atof(tmm[4].c_str()));	}
						if(ll>=6) if(tmm[5].size()>0) {	dip[dip_adv].setParam("chi",std::atof(tmm[5].c_str()));	}
						if(ll>=7) if(tmm[6].size()>0) {	dip[dip_adv].setParam("psi",std::atof(tmm[6].c_str()));	}
						if(ll<4){
							BLEXCEPTION(std::string(" A dipole must have at least 4 parameters..")+
									std::string("\n D <coupling> <spin 1> <spin 2> <alpha> <beta> <gamma>"))
						}
						dip[dip_adv].setCrystalAs();
					}
				}
				break;
			}
		}
  }
  return dip;
}

END_BL_NAMESPACE


