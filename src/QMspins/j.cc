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
 j.h --> spin interactions scaler coupling
*/


#include "QMspins/j.h"

BEGIN_BL_NAMESPACE


void J::reset(){
	iso_=0.;
	on1_=0;
	on2_=0;
	strong_=0;
}

J::J(){	reset();	}

void J::setParam(const std::string param, double value){
	if(param=="iso") 	{	iso_=value;}
	if(param=="on1") 	{on1_=int(value);}
	if(param=="on2") 	{on2_=int(value);}
	if(param=="strong") {strong_=int(value);}
	if(param=="1")		{on1_=int(value);}
	if(param=="2")		{on2_=int(value);}
}

double J::getParam(const std::string param){
	if(param=="iso") {	return iso_;}
	if(param=="del") { return 0.;}
	if(param=="eta") { return 0.;}
	if(param=="si" || param=="alpha") { return 0.;}
	if(param=="chi" || param=="beta") { return 0.;}
	if(param=="psi" || param=="gamma") { return 0.;}
	if(param=="chain1" || param=="on1") {return double(on1_);}
	if(param=="chain2" || param=="on2") {return double(on2_);}
	if(param=="strong") {return double(strong_);}
	return 0.;
}
/*set two J equal to each other
 */

//returns a std::string like Q(spin)
std::string J::name()
{	return "J("+itost(on1_)+","+itost(on2_)+")";	}


J J::operator=(const J &another){
	if(this==&another) return *this;
	iso_=another.iso_;
	on1_=another.on1_;
	on2_=another.on2_;
	strong_=another.strong_;
	T00=another.T00;
	return *this;
}

bool J::operator==(const J &another){
	if(this==&another) return true;
	if(iso_!=another.iso_) return false;
	if(on1_!=another.on1_) return false;
	if(on2_!=another.on2_) return false;
	if(strong_!=another.strong_) return false;
	return true;
}

bool J::operator!=(const J &another){
	if(this==&another) return true;
	if(iso_==another.iso_) return false;
	if(on1_==another.on1_) return false;
	if(on2_==another.on2_) return false;
	if(strong_==another.strong_) return false;
	return true;
}


std::ostream& operator<<(std::ostream &otr,const J &j){
	std::cout<<"#--------------J's--------------"<<std::endl;
  std::cout<<"#spin\tiso\tstrong\t"<<std::endl;
	std::cout<<"J"<<j.on1_<<","<<j.on2_<<"\t"<<j.iso_<<"\t"<<j.strong_<<std::endl;
	return otr;
}

void J::display(){
  std::cout<<"#--------------J's--------------"<<std::endl;
  std::cout<<"#spin\tiso\tstrong"<<std::endl;
	std::cout<<"J"<<on1_<<","<<on2_<<"\t"<<iso_<<"\t"<<strong_<<std::endl;
}

void J::write(std::ostream &otr){	//writes a line corresponding to the input
	otr<<"J "<<iso_<<" "<<on1_<<" "<<on2_<<" "<<strong_<<std::endl;
}


//*****************SCALER COUPLINGS********************
//generally in solids one can rarely see this anyway
//the cahnces of decifering between a I.I term versus
//Iz.Iz are quite small, but i will include the Iz.Iz
//terms incase somebody simulates a 50000Hz spinning speed
//or something as silly as that
//
//   H  =  J   Iz  *  Iz
//    J     i,j  i      j
//

void J::setSpinMats(SpinSys &A)
{
	T00=T_J(A, on1_, on2_,0,strong_);
}


hmatrix J::get_H(SpinSys &A){
 return iso_*sqrt(3.)*T00;
}

Vector<J> read_js(std::string in){
	std::ifstream infile(in.c_str());
	 if(infile.fail()){
    	BLEXCEPTION(std::string("J:: no spin file name like ")+in)
  }
 	return read_js(infile);
}

Vector<J> read_js(const char * in){
	std::ifstream infile(in);
	 if(infile.fail()){
		BLEXCEPTION(std::string("J:: no spin file name like ")+in)
  }
 	return read_js(infile);
}

Vector<J> read_js(std::ifstream &infile){
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
	return read_js(ii);
}

Vector<J> read_js(const Vector<std::string> &ss){
	Vector<J> js;

	int j_adv=-1;			//counter for num csa
	Vector<std::string> tmm;
	int i, len=ss.size(),ll=0;
	for(i=0;i<len;i++){
		tmm=parse_param(ss[i]);
		ll=tmm.size();
		if(!tmm.empty()){
			if(ll>=1){
				if(tmm[0]=="J" || tmm[0]=="j" || tmm[0]=="jcop" || tmm[0]=="Jcop"){
					if(ll>=4){
						j_adv++;
						js.push_back(J());
						js[j_adv].setParam("iso",std::atof(tmm[1].c_str()));
						js[j_adv].setParam("on1",std::atoi(tmm[2].c_str()));
						js[j_adv].setParam("on2",std::atoi(tmm[3].c_str()));
					}
					if(ll>=5){
						if(tmm[4].size()>0)  js[j_adv].setParam("strong",std::atoi(tmm[4].c_str()));
					}
					if(ll<=3){
							BLEXCEPTION(std::string(" A j-coupling must have at least 3 parameters..")+
									std::string(" J <coupling> <spin 1> <spin 2>"))
					}
				}
			}
		}
	}
  return js;
}


END_BL_NAMESPACE

