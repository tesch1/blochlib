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
 sys.cc --> methods manages the spin system, AND spin interactions (csa, j, dip, qua)
*/

#include "QMspins/sys.h"

BEGIN_BL_NAMESPACE


const void SolidSys::SizeErr()
{
	BLEXCEPTION(" Desired element NOT in System...")
}

SolidSys::SolidSys(const std::string in){ fname=in;	read(in); }
SolidSys::SolidSys(const char *in){	fname=std::string(in); read(in);	}

SolidSys::SolidSys(const Vector<std::string> &in){
	read(in);
}

SolidSys::SolidSys(const SolidSys &cp): SpinSys(cp),
	csa(cp.csa), dip(cp.dip), jcop(cp.jcop), qua(cp.qua)
{
	theRotations=cp.theRotations;
	fname=cp.fname;
	to_gamma=cp.to_gamma;  //true=do a third power angle
	roeq=cp.roeq;		//eq density matrix
	ZZZ=cp.ZZZ;
	III=cp.III;
	Bfield_=(cp.Bfield_);
	setBfield(Bfield_);
	setSpinMats();
	setCrystalAs();
	isdiag=cp.isdiag;

};

SolidSys SolidSys::operator=(const SolidSys &rhs)
{
	if(&rhs==this) return *this;
	SpinSys::operator=(rhs);
	theRotations=rhs.theRotations;
	Bfield_=rhs.Bfield_;
	curhamil=rhs.curhamil;
	csa=rhs.csa;
	dip=rhs.dip;
	jcop=rhs.jcop;
	qua=rhs.qua;
	fname=rhs.fname;
	//A=rhs.A;
	to_gamma=rhs.to_gamma;  //true=do a third power angle
	roeq=rhs.roeq;		//eq density matrix
 	ZZZ=rhs.ZZZ;
 	III=rhs.III;
 	setBfield(Bfield_);
	setSpinMats();
 	setCrystalAs();
 	isdiag=rhs.isdiag;
 	return *this;
}


void SolidSys::read(const std::string in){
	fname=std::string(in);
 	read(in.c_str());
}

void SolidSys::read(const char *in){
	fname=std::string(in);
	std::ifstream infile(in);
	//make sure you've given me the correct file name
 	if(infile.fail()){
    BLEXCEPTION(std::string("SolidSys:: no spin file name like ")+in)
  }
 	read(in, infile);

}

void SolidSys::read(const char * in, std::ifstream &infile){
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
	read(ii);
}

//rimatrix SolidSys::III=rimatrix();
//rdmatrix SolidSys::ZZZ=rdmatrix();

void SolidSys::read(const Vector<std::string> &ss){

/*here is am seting up the data input structures
 *they will assigned to the clsses as we go
 *the last step is to set up the spin system
 */
	int numspin=1, i, len=ss.size();

  //this is the little structure that reads in the type
  //to be later inputed to the spin system
  Vector<std::string> name;
  Vector<int> on1;
	bool gotone=false;
	name.push_back("1H");
	on1.push_back(0);

	 Vector<std::string> tmm;
	int ll;
	for(i=0;i<len;i++){
		tmm=parse_param(ss[i]);
		ll=tmm.size();

		if(!tmm.empty() && ll>=2){
			if(tmm[0]=="numspin" || tmm[0]=="spins"){
				if(ll>=2)	numspin=std::atoi(tmm[1].c_str());
			}else if(tmm[0]=="T" || tmm[0]=="type" || tmm[0]=="Type"){
				if(ll>=3){
      		   		if(gotone==false){
						name[0]=tmm[1];
				    	on1[0]=(std::atoi(tmm[2].c_str()));
				    	gotone=true;
					}else{
						name.push_back(tmm[1]);
						on1.push_back(std::atoi(tmm[2].c_str()));
					}

				}else{
					BLEXCEPTION(std::string("\n You have a 'T' specifier, but with the wrong number of inputs")+
					"\n there should be::  \"Type 1H 2\" ( {Type flag} {nucleous} {on which spin})")
				}
			}
		}
	}


/*now that we have our data, we need to set up the spin system
 *which all the hamiltonian types are pointing to, do some garbage collecting
 *and be done with the read in process
 **/

 //here is our base spin_system
//setup so that the '*A' in the hamitonian types
//have something to point to

 	//SpinSys B(numspin);
 	//A.resize(numspin);
	SpinSys::resize(numspin);
	len=name.size();

	for(i=0;i<len;i++){
		//A.isotope(on1[i], name[i]);
		isotope(on1[i], name[i]);
	}
	//A=B;


	//A.GenSpinOps();
	theRotations.RotationType=Rotations::Cartesian;
	GenSpinOps();
	qua=read_qua(ss);
	csa=read_csa(ss);
 	dip=read_dip(ss);
 	jcop=read_js(ss);
 	//the final initializer as all 'read' arguments are redirected here
 	//need to set up the matrices inside the spin types and an few other
 	//needed initial parameters
 	setBfield(Bfield_);
 	setRoEQ();
	setMats(*this);
	to_gamma=check_gamma();
 	setSpinMats();
 	setCrystalAs();


}

void SolidSys::setBfield(double Bf)
{
	Bfield_=Bf;
	for(int i=0;i<qua.size();i++)	qua[i].Bfield(Bf);
}

//given an input std::string will see if that
//type of interaction and spins
//are acctually available (i.e. weather or
//not a D01 (dipole between spin 0 and spin 1)
//exsits
bool SolidSys::check_spin_avail(std::string in){
	int i, len;
	int it1, it2;
	std::string moo="";
	if(in[0]=='C'){
		if(csa.empty())		return false;
		if(!isdigit(in[1]))	return false;
		else{
			len=csa.size();
			moo+=in[1];
			it1=std::atoi(moo.c_str());
			moo="";
			for(i=0;i<len;i++){
				if(it1==csa[i].spinON()) return true;
			}
			return false;
		}
	}else if(in[0]=='D'){
		if(dip.empty())	 	return false;
		if(!isdigit(in[1])||!isdigit(in[2])) return false;
		else{
			len=dip.size();
			moo+=in[1];
			it1=std::atoi(moo.c_str());
			moo="";
			moo+=in[2];
			it2=std::atoi(moo.c_str());
			moo="";
			for(i=0;i<len;i++){
				if((it1==dip[i].spinON1() && it2==dip[i].spinON2()) ||(it2==dip[i].spinON1() && it1==dip[i].spinON2())  ) return true;
			}
			return false;
		}
	}else if(in[0]=='Q'){
		if(qua.empty())		return false;
		if(!isdigit(in[1]))	return false;
		else{
			len=qua.size();
			moo+=in[1];
			it1=std::atoi(moo.c_str());
			moo="";
			for(i=0;i<len;i++){
				if(it1==qua[i].spinON()) return true;
			}
			return false;
		}

	}else if(in[0]=='J'){
		if(jcop.empty())	return false;
		if(!isdigit(in[1])||!isdigit(in[2]))	return false;
		else{
			len=jcop.size();
			moo+=in[1];
			it1=std::atoi(moo.c_str());
			moo="";
			moo+=in[2];
			it2=std::atoi(moo.c_str());
			moo="";
			for(i=0;i<len;i++){
				if((it1==jcop[i].spinON1() && it2==jcop[i].spinON2()) ||(it2==jcop[i].spinON1() && it1==jcop[i].spinON2())  ) return true;
			}
			return false;
		}
	}
	return false;
}

double SolidSys::getSpinParam(std::string in){
	bool tt=check_spin_avail(in);
	if(tt!=false){
		int le=in.length(); int on1=0, on2=0;
		std::string bit=""; char tmp[2]; tmp[1]='\0';
		if(le>=3 && in[0]=='C' || in[0]=='Q'){
			if(in[2]=='e'){
				bit="eta";
				tmp[0]=in[1];
				on1=std::atoi(tmp);
				return getSpinParam(in[0], bit, on1, on2);
			}else if(in[2]=='d'){
				bit="del";
				tmp[0]=in[1];
				on1=std::atoi(tmp);
				return getSpinParam(in[0], bit, on1, on2);
			}else if(in[2]=='i' && in[0]!='Q'){
				bit="iso";
				tmp[0]=in[1];
				on1=std::atoi(tmp);
				return getSpinParam(in[0], bit, on1, on2);
			}else if(in[2]=='a' ){		//alpha
				bit="alpha";
				tmp[0]=in[1];
				on1=std::atoi(tmp);
				return getSpinParam(in[0], bit, on1, on2);
			}else if(in[2]=='b' ){		//beta
				bit="beta";
				tmp[0]=in[1];
				on1=std::atoi(tmp);
				return getSpinParam(in[0], bit, on1, on2);
			}else if(in[2]=='g' ){		//gamma
				bit="gamma";
				tmp[0]=in[1];
				on1=std::atoi(tmp);
				return getSpinParam(in[0], bit, on1, on2);
			}else if(le>=6){
				if(in[2]=='s' && in[3]=='i' && in[4]=='g'){
					if(in[5]=='1'){
						bit="sig1";
						tmp[0]=in[1];
						on1=std::atoi(tmp);
						return getSpinParam(in[0], bit, on1, on2);
					}else if(in[5]=='2'){
						bit="sig2";
						tmp[0]=in[1];
						on1=std::atoi(tmp);
						return getSpinParam(in[0], bit, on1, on2);
					}else if(in[5]=='3'){
						bit="sig3";
						tmp[0]=in[1];
						on1=std::atoi(tmp);
						return getSpinParam(in[0], bit, on1, on2);
					}
				}
			}
		//these are defaults values i.e. a "C1" would simply set the isotropic shift
		}else if(in[0]=='Q'){
			bit=="del";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			return getSpinParam(in[0], bit, on1, on2);
		}else if(in[0]=='C'){
			bit="iso";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			return getSpinParam(in[0], bit, on1, on2);
		}else if(in[0]=='J'){
			bit="iso";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			tmp[0]=in[2];
			on2=std::atoi(tmp);
			return getSpinParam(in[0], bit, on1, on2);
		}else if(in[0]=='D'){
			bit="del";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			tmp[0]=in[2];
			on2=std::atoi(tmp);
			return getSpinParam(in[0], bit, on1, on2);
		}
	}else{
		BLEXCEPTION(std::string(" Requested spin for parameter:")+in+" Does NOT exsist...")
	}
	return 0.;
}


double SolidSys::getSpinParam(char type, std::string bit,int on1, int on2){
	int i=0;
	if(type=='J'){
		for(i=0;i<jcop.size();i++){
			if((on1==jcop[i].spinON1()&&on2==jcop[i].spinON2())||(on1==jcop[i].spinON2()&&on2==jcop[i].spinON1())){
				 return jcop[i].getParam(bit);
			}else if((i+1)==jcop.size()){
				 BLEXCEPTION(std::string("no ")+type+" coupling exsits between spin "+itost(on1)+" and "+itost(on2))
			}
		}
	}else if(type=='D'){
		for(i=0;i<dip.size();i++){
			if((on1==dip[i].spinON1()&&on2==dip[i].spinON2())||(on1==dip[i].spinON2()&&on2==dip[i].spinON1())){
				return dip[i].getParam(bit);
				//std::cerr<<my_func.dip[i].couple()<<std::endl;
			}else if((i+1)==dip.size()){
				BLEXCEPTION(std::string("no Dipole coupling exsits between spin ")+itost(on1)+" and "+itost(on2))
			}
		}
	}else if(type=='Q'){
		return getSpinParam(type, bit,on1);
	}else if(type=='C'){
		return getSpinParam(type, bit, on1);
	}
	return 0.;
}


double SolidSys::getSpinParam(char type, std::string bit,int on1){
	int i=0;
	if(type=='C'){
		for(i=0;i<csa.size();i++){
			if(on1==csa[i].spinON()){
				return csa[i].getParam(bit);
			}
		}
		if((i+1)==csa.size()){
			BLEXCEPTION(std::string("no CSA coupling exsits on spin ")+itost(on1))
		}
	}else if(type=='Q'){
		for(i=0;i<qua.size();i++){
			if(on1==qua[i].spinON()){
				return qua[i].getParam(bit);
			}else if((i+1)==qua.size()){
				BLEXCEPTION(std::string("no Quad coupling exsits on spin ")+itost(on1))
			}
		}
	}
	return 0.;
}

//these guys set a spin parameter
//given a char type (D, Q, J, C)
//a bit (iso, eta, del)
//a spin number (0...numspins)
//and the num to set the coupling to
////D and J require two spin numbers
//
//this one parses an input std::string and checks to see if
//i can acctullay set the parameter
void SolidSys::setSpinParam(std::string in, double nm){
	bool tt=check_spin_avail(in);
	if(tt==false){
		BLEXCEPTION(std::string(" Your spin parameter ")+in+" Ain't no good")
	}
	int le=in.length(); int on1=0, on2=0;
	std::string bit=""; char tmp[2]; tmp[1]='\0';
	if(le>=3 && in[0]=='C' || in[0]=='Q'){
		if(in[2]=='e'){
			bit="eta";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			setSpinParam(in[0], bit, on1, on2, nm);
		}else if(in[2]=='d'){
			bit="del";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			setSpinParam(in[0], bit, on1, on2, nm);
		}else if(in[2]=='i' && in[0]!='Q'){
			bit="iso";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			setSpinParam(in[0], bit, on1, on2, nm);
		}else if(in[2]=='a' ){		//alpha
			bit="alpha";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			setSpinParam(in[0], bit, on1, on2, nm);
		}else if(in[2]=='b' ){		//beta
			bit="beta";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			setSpinParam(in[0], bit, on1, on2, nm);
		}else if(in[2]=='g' ){		//gamma
			bit="gamma";
			tmp[0]=in[1];
			on1=std::atoi(tmp);
			setSpinParam(in[0], bit, on1, on2, nm);
		}else if(le>=6){
			if(in[2]=='s' && in[3]=='i' && in[4]=='g'){
				if(in[5]=='1'){
					bit="sig1";
					tmp[0]=in[1];
					on1=std::atoi(tmp);
					setSpinParam(in[0], bit, on1, on2, nm);
				}else if(in[5]=='2'){
					bit="sig2";
					tmp[0]=in[1];
					on1=std::atoi(tmp);
					setSpinParam(in[0], bit, on1, on2, nm);
				}else if(in[5]=='3'){
					bit="sig3";
					tmp[0]=in[1];
					on1=std::atoi(tmp);
					setSpinParam(in[0], bit, on1, on2, nm);
				}
			}
		}
	//these are defaults values i.e. a "C1" would simply set the isotropic shift
	}else if(in[0]=='Q'){
		bit=="del";
		tmp[0]=in[1];
		on1=std::atoi(tmp);
		setSpinParam(in[0], bit, on1, on2, nm);
	}else if(in[0]=='C'){
		bit="iso";
		tmp[0]=in[1];
		on1=std::atoi(tmp);
		setSpinParam(in[0], bit, on1, on2, nm);
	}else if(in[0]=='J'){
		bit="iso";
		tmp[0]=in[1];
		on1=std::atoi(tmp);
		tmp[0]=in[2];
		on2=std::atoi(tmp);
		setSpinParam(in[0], bit, on1, on2, nm);
	}else if(in[0]=='D'){
		bit="del";
		tmp[0]=in[1];
		on1=std::atoi(tmp);
		tmp[0]=in[2];
		on2=std::atoi(tmp);
		setSpinParam(in[0], bit, on1, on2, nm);
	}

}

void SolidSys::setSpinParam(char type, std::string bit,int on1, int on2, double num){
	int i=0;
	if(type=='J'){
		for(i=0;i<jcop.size();i++){
			if((on1==jcop[i].spinON1()&&on2==jcop[i].spinON2())||(on1==jcop[i].spinON2()&&on2==jcop[i].spinON1())){
				 jcop[i].setParam(bit,num);
			}else if((i+1)==jcop.size()){
				BLEXCEPTION(std::string("no ")+type+" coupling exsits between spin "+itost(on1)+" and "+itost(on2))
			}
		}
	}else if(type=='D'){
		for(i=0;i<dip.size();i++){
			if((on1==dip[i].spinON1()&&on2==dip[i].spinON2())||(on1==dip[i].spinON2()&&on2==dip[i].spinON1())){
				dip[i].setParam(bit,num);
				//std::cerr<<my_func.dip[i].couple()<<std::endl;
			}else if((i+1)==dip.size()){
				BLEXCEPTION(std::string("no Dipole coupling exsits between spin ")+itost(on1)+" and "+itost(on2))
			}
		}
	}else if(type=='Q'){
		setSpinParam(type, bit,on1, num);
	}else if(type=='C'){
		setSpinParam(type, bit, on1, num);
	}
}


void SolidSys::setSpinParam(char type, std::string bit,int on1, double num){
	int i=0;
	if(type=='C'){
		for(i=0;i<csa.size();i++){
			if(on1==csa[i].spinON()){
				csa[i].setParam(bit, num);
				return;
			}
		}
		if((i+1)==csa.size()){
			BLEXCEPTION(std::string("no CSA coupling exsits on spin ")+itost(on1))
		}
	}else if(type=='Q'){
		for(i=0;i<qua.size();i++){
			if(on1==qua[i].spinON()){
				qua[i].setParam(bit,num);
				return;
			}else if((i+1)==qua.size()){
				BLEXCEPTION(std::string("no Quad coupling exsits on spin ")+itost(on1))
			}
		}
	}
}

void SolidSys::addDip(Dip lhs)
{
	if(lhs.on1()>size() || lhs.on2()>size())
	//{	A=SpinSys(max(lhs.on1(), lhs.on2()));	}
	{	resize(std::max(lhs.on1(), lhs.on2()));	}

	dip.push_back(lhs);
	setRoEQ();
	setMats(*this);
	setSpinMats();
	setCrystalAs();
}

void SolidSys::addCsa(Csa lhs)
{
	if(lhs.on()>size())
	//{	A=SpinSys(lhs.on());	}
	{	resize(lhs.on());	}

	csa.push_back(lhs);
	setRoEQ();
	setMats(*this);
	setSpinMats();
	setCrystalAs();
}

void SolidSys::addJ(J lhs)
{
	if(lhs.on1()>size() || lhs.on2()>size())
	//{	A=SpinSys(max(lhs.on1(), lhs.on2()));	}
	{	resize(std::max(lhs.on1(), lhs.on2()));	}

	jcop.push_back(lhs);
	setRoEQ();
	setMats(*this);
	setSpinMats();
	setCrystalAs();
}

void SolidSys::addQ(Qua lhs)
{
	if(lhs.on()>size())
	//{	A=SpinSys(lhs.on());	}
	{	resize(lhs.on());	}

	qua.push_back(lhs);
	setBfield(Bfield_);
	setRoEQ();
	setMats(*this);
	setSpinMats();
	setCrystalAs();
}

void SolidSys::setRoEQ(){
	//if(A.homonuclear()){	roeq=A.Fz();	return;	}
	if(homonuclear()){	roeq=Fz();	return;	}


	//double g0=A.gamma(0), g1;
	double g0=gamma(0), g1;
	//roeq=A.Iz(0);
	roeq=Iz(0);
	int i=0;
	//for(i=1;i<A.spins();i++){
	for(i=1;i<spins();i++){
		//g1=A.gamma(i);
		g1=gamma(i);
	//	roeq+=(g1/g0)*A.Iz(i);
		roeq+=(g1/g0)*Iz(i);
	}
}

void SolidSys::setMats(const SpinSys &As){
	//III=Fe(As)/As.spins();
//	ZZZ=A.F0();
//	III=A.Fe();
	ZZZ=F0();
	III=Fe();
}

void SolidSys::setSpinMats(){
	 int i=0;
	isdiag=true;
	for(i=0;i<(csa.size());i++){
		csa[i].setSpinMats(*this);
	}
	for(i=0;i<(dip.size());i++){
		dip[i].setSpinMats(*this);
		if(!dip[i].is_hetero(*this)) isdiag=false;
	}
	for(i=0;i<(jcop.size());i++){
		jcop[i].setSpinMats(*this);
		isdiag=false;
	}
	for(i=0;i<(qua.size());i++){
		qua[i].setSpinMats(*this);
		//qua[i].F4=T_Qq(A, qua[i].spinON(),2,4);

	}

	setMats(*this);
}


void SolidSys::setCrystalAs(){
	int i=0;
	for(i=0;i<(csa.size());i++){
		csa[i].setCrystalAs();
	}
	for(i=0;i<(dip.size());i++){
		dip[i].setCrystalAs();
	}

	for(i=0;i<(qua.size());i++){
		qua[i].setCrystalAs();
	}
}


bool SolidSys::check_gamma(){
	int i=0;
	for(i=0;i<(qua.size());i++){
		if(qua[i].eta_cop()!=0) return true;
	}
	for(i=0;i<(csa.size());i++){
		if(csa[i].eta_cop()!=0) return true;
	}
	return false;
}

void SolidSys::setPowderAngles(double theta, double phi, double gam)
{	theRotations.setPowderAngles(phi, theta, gam);	}

void SolidSys::setRotorAngles(double wr, double be, double chi)
{	theRotations.setRotorAngles(wr, be, chi);	}

void SolidSys::setAngles(double wr, double be, double theta, double phi)
{
	setPowderAngles(theta, phi);
	setRotorAngles(wr, be);
}

//HAMILTONINAS
hmatrix &SolidSys::Hamiltonian(double alpha, double beta, double theta, double phi, double t1, double t2)
{
	curhamil=ZZZ;
	//first we set the wigner elements way back in the 'space_ten.h' section
	//std::cout<<powderAs;

	setAngles(alpha, beta, theta, phi);
	int i=0;

	//std::cout<<std::endl<<"SIZE:"<<csa.size()<<std::endl;
	for(i=0;i<csa.size();++i){
		curhamil+=csa[i].Hamiltonian(*this,theRotations);
		//std::cout<<curhamil<<std::endl;
	}
	for(i=0;i<dip.size();++i){
		curhamil+=dip[i].Hamiltonian(*this,theRotations);
	}
	for(i=0;i<jcop.size();++i){
		curhamil+=jcop[i].Hamiltonian(*this);
	}
	for(i=0;i<qua.size();++i){
		curhamil+=qua[i].Hamiltonian(*this,theRotations);
	}

	return curhamil;
}

//HAMILTONINAS
hmatrix &SolidSys::H()
{	return Hamiltonian(0.0, 1.0);	}

hmatrix &SolidSys::Hamiltonian()
{	return Hamiltonian(0.0, 1.0);	}

//HAMILTONINAS
hmatrix &SolidSys::Hamiltonian(double t1, double t2, double wr)
{
	curhamil=ZZZ;

	int i=0;
	double tmid=(t2-t1)/2.0 + t1;
	setRotorAngles(tmid*wr*PI2, theRotations.beta);
	for(i=0;i<csa.size();++i){
		curhamil+=csa[i].H(*this,theRotations);
	}
	for(i=0;i<dip.size();++i){
		curhamil+=dip[i].H(*this,theRotations);
	}
	for(i=0;i<jcop.size();++i){
		curhamil+=jcop[i].H(*this);
	}
	for(i=0;i<qua.size();++i){
		curhamil+=qua[i].H(*this,theRotations);
	}

	return curhamil;
}


//calls up the hamils display functions
//and displays the info about the spin system


void SolidSys::print(std::ostream &oo)
{

	oo<<"#--------------------------------"<<std::endl;
	oo<<"#--------spin parameters---------"<<std::endl;
	oo<<"#--------------------------------"<<std::endl;
	//std::cout<<"you have "<<A.spins()<<" spins in your system"<<std::endl;
	//	std::cout<<A;
	oo<<"#you have "<<spins()<<" spins in your system"<<std::endl;
	oo<<*spinSys()<<std::endl;
	int i=0;
	for(i=0;i<jcop.size();i++){ oo<<jcop[i]; }
	for(i=0;i<dip.size();i++){ oo<<dip[i]; }
	for(i=0;i<qua.size();i++){ oo<<qua[i]; }
	for(i=0;i<csa.size();i++){ oo<<csa[i]; }
}

void SolidSys::display()
{	print(std::cout);	}

std::ostream &operator<<(std::ostream &otr, SolidSys &out){
  out.print(otr);
  return otr;
}

void SolidSys::write(std::string fname){
	 std::ofstream out(fname.c_str());
	write(out);
}

void SolidSys::write(std::ostream &otr){
	 int i;
	int j;
	//otr<<"numspin "<<A.spins()<<std::endl;
	//for(j=0;j<A.spins();j++){
	//	otr<<"T "<<A.symbol(j)<<" "<<j<<std::endl;
	//}
	otr<<"numspin "<<spins()<<std::endl;
	for(j=0;j<spins();j++){
		otr<<"T "<<symbol(j)<<" "<<j<<std::endl;
	}
	for(i=0;i<(csa.size());i++){
		csa[i].write(otr);
	}
	for(i=0;i<(dip.size());i++){
		dip[i].write(otr);
	}

	for(i=0;i<(qua.size());i++){
		qua[i].write(otr);
	}
	for(i=0;i<(jcop.size());i++){
		jcop[i].write(otr);
	}
}

END_BL_NAMESPACE

