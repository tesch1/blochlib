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
	csa.cc --> methods for spin interaction CSA
*/

#include "QMspins/csa.h"

BEGIN_BL_NAMESPACE


void Csa::reset(){
	iso_=0.;
	del_=0.;
	eta_=0.;
	si_=0.;
	chi_=0.;
	psi_=0.;
	chain_=0;
	sig1_=sig2_=sig3_=0;
	crystalAs.resize(5,0);
}

Csa::Csa(){	reset();	}

void Csa::setParam(const std::string param, double value){
	if(param=="iso") {
		iso_=value;
		setSigma();
		setCrystalAs();
	}
	if(param=="del") {
		del_=value;
		setSigma();
		setCrystalAs();
	}
	if(param=="eta") {
		eta_=value;
		setSigma();
		if(eta_>1. || eta_<0){
			std::cerr<<std::endl<<*this<<std::endl;
			BLEXCEPTION(" 'eta' set to a negative or value >1...")
		}
		setCrystalAs();
	}
	if(param=="sig1"){
		sig1_=value;
		iso_=(sig1_+sig2_+sig3_)/3.;
		if(abs(sig3_)<=abs(sig1_)) swap_(sig3_, sig1_);
		if(abs(sig3_)<=abs(sig2_)) swap_(sig3_, sig2_);
		if(abs(sig1_)<=abs(sig2_)) swap_(sig2_, sig1_);
		del_=sig3_-iso_;
		eta_=(sig1_-sig2_)/del_;
		if(eta_<0) eta_=-eta_;
		//sig2=iso+del/2.*(1-eta);
		if(eta_>1. || eta_<0){
			std::cerr<<std::endl<<*this<<std::endl;
			BLEXCEPTION(" 'eta' set to a negative or value >1...")
		}
		setCrystalAs();
	}
	if(param=="sig2") {
		sig2_=value;
		iso_=(sig1_+sig2_+sig3_)/3.;
		if(abs(sig3_)<=abs(sig1_)) swap_(sig3_, sig1_);
		if(abs(sig3_)<=abs(sig2_)) swap_(sig3_, sig2_);
		if(abs(sig1_)<=abs(sig2_)) swap_(sig2_, sig1_);
		del_=sig3_-iso_;
		eta_=(sig1_-sig2_)/del_;
		if(eta_<0) eta_=-eta_;
		//sig1=iso+del/2.*(1+eta);
		if(eta_>1. || eta_<0){
			std::cerr<<std::endl<<*this<<std::endl;
			BLEXCEPTION(" 'eta' set to a negative or value >1...")
		}
		setCrystalAs();
	}
	if(param=="sig3") {
		sig3_=value;
		iso_=(sig1_+sig2_+sig3_)/3.;
		if(abs(sig3_)<=abs(sig1_)) swap_(sig3_, sig1_);
		if(abs(sig3_)<=abs(sig2_)) swap_(sig3_, sig2_);
		if(abs(sig1_)<=abs(sig2_)) swap_(sig2_, sig1_);
		del_=sig3_-iso_;
		eta_=(sig1_-sig2_)/del_;
		if(eta_<0) eta_=-eta_;
		if(eta_>1. || eta_<0){
			std::cerr<<std::endl<<*this<<std::endl;
			BLEXCEPTION(" 'eta' set to a negative or value >1...")
		}
		setCrystalAs();
	}
	if(param=="si" || param=="alpha") {si_=value; setCrystalAs();}
	if(param=="chi" || param=="beta") {chi_=value; setCrystalAs();}
	if(param=="psi" || param=="gamma") {psi_=value; setCrystalAs();}
	if(param=="chain") {chain_=int(value);}
}

void Csa::setSigma()
{
	sig1_=iso_+del_/2.*(1+eta_);
	sig2_=iso_+del_/2.*(1-eta_);
	sig3_=iso_-del_;
	if(abs(sig3_)<=abs(sig1_)) swap_(sig3_, sig1_);
	if(abs(sig3_)<=abs(sig2_)) swap_(sig3_, sig2_);
	if(abs(sig1_)<=abs(sig2_)) swap_(sig2_, sig1_);
}


double Csa::getParam(const std::string param){
	if(param=="iso") {	return iso_;}
	if(param=="del") { return del_;}
	if(param=="eta") { return eta_;}
	if(param=="sig1"){ return sig1_;}
	if(param=="sig2"){ return sig2_;}
	if(param=="sig3"){ return sig3_;}
	if(param=="si" || param=="alpha") { return si_;}
	if(param=="chi" || param=="beta") { return chi_;}
	if(param=="psi" || param=="gamma") { return psi_;}
	if(param=="chain" || param=="on") {return double(chain_);}
	return 0.;
}

//returns a std::string like 'C(spin1)'
std::string Csa::name()
{	return "C("+itost(chain_)+")";	}

void Csa::setSpinMats(SpinSys &A)
{
	T20=T20_CSA(A, spinON());
	T00=T00_CSA(A, spinON());
}

//display the bits in our linked list

void Csa::display(){

	std::cout<<"#-------------CSA--------------"<<std::endl;
 	std::cout<<"#spin\tiso\tdel\teta\talpha\tbeta\tgamma"<<std::endl;
	std::cout<<chain_<<"\t"<<iso_<<"\t"<<del_<<"\t"<<eta_<<"\t"<<
			si_<<"\t"<<chi_<<"\t"<<psi_<<std::endl;

}

void Csa::write(std::ostream &otr){	//writes a line corresponding to the input
	otr<<"C "<<iso_<<" "<<del_<<" "<<eta_<<" "<<chain_<<" "<<si_<<" "<<chi_<<" "<<psi_<<std::endl;
}


Csa Csa::operator=(const Csa &another){
	if(this==&another) return *this;
	iso_=another.iso_;
	del_=another.del_;
	eta_=another.eta_;
	si_=another.si_;
	chi_=another.chi_;
	psi_=another.psi_;
	sig1_=another.sig1_;
	sig2_=another.sig2_;
	sig3_=another.sig3_;
	chain_=another.chain_;
	crystalAs=another.crystalAs;
	T20=another.T20;
	T00=another.T00;
	return *this;
}

bool Csa::operator==(const Csa &another)
{
	if(iso_!=another.iso_) return false;
	if(del_!=another.del_)return false;
	if(eta_!=another.eta_)return false;
	if(si_!=another.si_)return false;
	if(chi_!=another.chi_)return false;
	if(psi_!=another.psi_)return false;
	if(chain_!=another.chain_)return false;
	return true;
}

bool Csa::operator!=(const Csa &another)
{
	if(iso_!=another.iso_) return true;
	if(del_!=another.del_)return true;
	if(eta_!=another.eta_)return true;
	if(si_!=another.si_)return true;
	if(chi_!=another.chi_)return true;
	if(psi_!=another.psi_)return true;
	if(chain_!=another.chain_)return true;
	return false;
}


std::ostream& operator<<(std::ostream &otr, const Csa &csa){
	otr<<"#-------------CSA--------------"<<std::endl;
 	otr<<"#spin\tiso\tdel\teta\talpha\tbeta\tgamma"<<std::endl;
	otr<<csa.chain_<<"\t"<<csa.iso_<<"\t"<<csa.del_<<"\t"<<csa.eta_<<"\t"<<
			csa.si_<<"\t"<<csa.chi_<<"\t"<<csa.psi_<<std::endl;
	return otr;

}


void Csa::setCrystalAs(){
	if(si_==0&&chi_==0&&psi_==0){
		crystalAs[0]=A_CSA_pas(iso_,del_,eta_,2,-2);
		crystalAs[1]=complex(0.,0.);
		crystalAs[2]=A_CSA_pas(iso_,del_,eta_,2,0);
		crystalAs[3]=complex(0.,0.);
		crystalAs[4]=A_CSA_pas(iso_,del_,eta_,2,2);
	}else{
		double prad=si_*DEG2RAD;
		double thrad=chi_*DEG2RAD;
		double grad=psi_*DEG2RAD;
		//A2-2
		crystalAs[0]=wigner2(prad,thrad,grad, -2,2)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad, -2,0)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad, -2,-2)*A_CSA_pas(iso_,del_,eta_,2,-2);
		//A2-1
		crystalAs[1]=wigner2(prad,thrad,grad, -1,2)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad, -1,0)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad, -1,-2)*A_CSA_pas(iso_,del_,eta_,2,-2);
		//A20
		crystalAs[2]=wigner2(prad,thrad,grad, 0,2)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad, 0,0)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad, 0,-2)*A_CSA_pas(iso_,del_,eta_,2,-2);
		//A21
		crystalAs[3]=wigner2(prad,thrad,grad, 1,2)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad, 1,0)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad, 1,-2)*A_CSA_pas(iso_,del_,eta_,2,-2);
		//A22
		crystalAs[4]=wigner2(prad,thrad,grad, 2,2)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad, 2,0)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad, 2,-2)*A_CSA_pas(iso_,del_,eta_,2,-2);

	/*	crystalAs[0]=wigner2(prad,thrad,grad, 2,-2)*A_CSA_pas(iso_,del_,eta_,2,2)+
						wigner2(prad,thrad,grad,0, -2)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad, -2,-2)*A_CSA_pas(iso_,del_,eta_,2,-2);
		//A2-1
		crystalAs[1]=wigner2(prad,thrad,grad, 2,-1)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad, 0,-1)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad, -2,-1)*A_CSA_pas(iso_,del_,eta_,2,-2);
		//A20
		crystalAs[2]=wigner2(prad,thrad,grad,2, 0)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad,0, 0)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad,-2, 0)*A_CSA_pas(iso_,del_,eta_,2,-2);
		//A21
		crystalAs[3]=wigner2(prad,thrad,grad,2, 1)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad,0, 1)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad,-2, 1)*A_CSA_pas(iso_,del_,eta_,2,-2);
		//A22
		crystalAs[4]=wigner2(prad,thrad,grad,2, 2)*A_CSA_pas(iso_,del_,eta_,2,2)+
					wigner2(prad,thrad,grad,0, 2)*A_CSA_pas(iso_,del_,eta_,2,0)+
					wigner2(prad,thrad,grad,-2, 2)*A_CSA_pas(iso_,del_,eta_,2,-2);
*/
	}
	rmatrix rot=rotationMatrix3D(si_*DEG2RAD, chi_*DEG2RAD, psi_*DEG2RAD);
	PAS_.resize(3,3);
	PAS_.fill(0.0);
	PAS_(0,0)=del_*(eta_-1.0)/2.0;
	PAS_(1,1)= -del_*(eta_+1.0)/2.0;
	PAS_(2,2)=1.0*del_;
	PAS_=rot*
		 PAS_*
		 transpose(rot);

}


dmatrix Csa::get_H(SpinSys &A, Rotations &rr){
	dmatrix temp=A.F0();
	complex rot(0.,0.);
	static rmatrix rotm(3,3);
	//static double curPhi=-1010101.0, curTheta=0.0, CurGamma=0.0, curA=0.0, curB=0.0;
	//int m,mp;

	if(iso_!=0){
		rot=A_CSA_pas(iso_, 0., 0., 0,0);
		temp=T00*rot;
		if(del_==0.0){ return temp; }
	}

	if(rr.beta==0.0){
		if(eta_==0.0 && del_==0.0){	return temp; }
		else{
			/*rot=rr.powderAs(2,0)*crystalAs[0]+
				rr.powderAs(2,1)*crystalAs[1]+
				rr.powderAs(2,2)*crystalAs[2]+
				rr.powderAs(2,3)*crystalAs[3]+
				rr.powderAs(2,4)*crystalAs[4];
			temp+=T20*rot;
			return temp;*/
			rotm=rr.cartPowder*PAS_*transpose(rr.cartPowder);
			return rotm(2,2)*T20+temp;
		}
	}else {
		/*rot=ZeroType<complex>::zero();
		for(m=0;m<=4;m++){
			for(mp=0;mp<=4;mp++){
				rot+=rr.spinnerAs[m]*rr.powderAs(mp,m)*crystalAs[mp];
			}
		}
		temp+=rot*T20;
		return temp;*/
		rotm=rr.cartSpin*rr.cartPowder;
		rotm=rotm*PAS_*transpose(rotm);
		return rotm(2,2)*T20+temp;
	}
	return temp;
}


dmatrix Csa::get_H(SpinSys &A,double alpha, double beta, double theta, double phi, double gam){
	Rotations rots(phi,theta, gam, alpha, beta);
	return get_H(A,rots);
}


dmatrix Csa::get_HINT(SpinSys &A,double wr, double tau1, double tau2, double beta, double theta, double phi){
	dmatrix temp=A.F0();
	complex rot=complex(0,0);
  if(wr==0.0 || (tau2-tau1)==0.0){ return get_H(A, 0,0,theta, phi); }
	if(iso_!=0){
		rot=A_CSA_pas(iso_, 0., 0., 0,0)*2.*Pi*(tau2-tau1);
 		temp+=T00*rot;
 		if(del_==0.0){ return temp; }
	}
	if(si_==0.0 && chi_==0.0 && psi_==0.0){
		if(eta_==0.0){

	  			rot=-((pi*tau1*(wr) - pi*tau2*( wr) +
        3.*pi*(tau1 - tau2)*( wr)*cos(2.*beta) +
        3.*pi*(tau1 - tau2)*( wr)*(1. + 3.*cos(2.*beta))*
         cos(2.*theta) + 6.*sin(2.*beta)*sin(2.*theta)*
         sin(phi + 2.*pi*tau1*( wr)) +
        3.*pow(sin(beta),double(2.))*pow(sin(theta),double(2.))*
         sin(2.*(phi + 2.*pi*tau1*( wr))) -
        6.*sin(2.*beta)*sin(2.*theta)*sin(phi + 2.*pi*tau2*( wr)) -
        3.*pow(sin(beta),double(2.))*pow(sin(theta),double(2.))*
         sin(2.*(phi + 2.*pi*tau2*( wr)))))/(8.*( wr));

  		temp+=A_CSA_pas(0, del_, 0, 2,0)*rot*T20;
  		return temp;
  	}
}

  return temp;
}

Vector<Csa> read_csa(std::string in){
	std::ifstream infile(in.c_str());
	 if(infile.fail()){
		BLEXCEPTION(std::string("CSA:: no spin file name like ")+in)
  }
 	return read_csa(infile);
}

Vector<Csa> read_csa(const char * in){
	std::ifstream infile(in);
	 if(infile.fail()){
		BLEXCEPTION(std::string("CSA:: no spin file name like ")+in)
  }
 	return read_csa(infile);
}

Vector<Csa> read_csa(std::ifstream & infile){
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
	return read_csa(ii);
}

Vector<Csa> read_csa(const Vector<std::string> &ss){

	Vector<Csa> csa_in;
	int csa_adv=-1;			//counter for num csa
	Vector<std::string> tmm;
	int i, len=ss.size(), ll=0;
	for(i=0;i<len;i++){
		tmm=parse_param(ss[i]);
		ll=tmm.size();
		if(!tmm.empty()){
			switch(tmm[0][0]){
			case '#': case '\0': case '\n':   //ignore lines that start with these
			  break;
			default:
				if(ll>=1){
					if(tmm[0]=="C" || tmm[0]=="CSA" || tmm[0]=="csa"){
						if(ll>=5){
							csa_adv++;
							csa_in.push_back(Csa());
							csa_in[csa_adv].setParam("iso",std::atof(tmm[1].c_str()));
							csa_in[csa_adv].setParam("del",std::atof(tmm[2].c_str()));
							csa_in[csa_adv].setParam("eta",std::atof(tmm[3].c_str()));
							csa_in[csa_adv].setParam("chain",std::atoi(tmm[4].c_str()));
						}
						if(ll>=6) if(tmm[5].size()>0) {	csa_in[csa_adv].setParam("si",std::atof(tmm[5].c_str()));	}
						if(ll>=7) if(tmm[6].size()>0)  {	csa_in[csa_adv].setParam("chi",std::atof(tmm[6].c_str()));	}
						if(ll>=8) if(tmm[7].size()>0) {	csa_in[csa_adv].setParam("psi",std::atof(tmm[7].c_str()));	}
						if(ll<=4){
							BLEXCEPTION(std::string(" A CSA must have at least 4 parameters..")+
									" C <iso> <del> <eta> <spin> <alpha> <beta> <gamma>")
						}
						csa_in[csa_adv].setCrystalAs();
					}

				}
				break;
			}
		}
  }
  return csa_in;
}

END_BL_NAMESPACE

