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
	powder.cc--> generates the powder angles from various methods
	like, zcw (the best one), alderman, sophie, planar, spherical,
	and 'none' for liquid type sims
*/


#ifndef _powder_cc_
#define _powder_cc_ 1

#include "QMspins/driver/powder.h"

BEGIN_BL_NAMESPACE



void powder::__reset(){
	_tstep=1;
	_pstep=1;
	_gstep=1;
	_method=None;
	_fname="NULL";
	_underlimit=true;
	_ct=0;
}

powder::pTypes powder::AssignM(std::string method)
{
	if(method=="alderman")	return alderman;
	else if(method=="sphere")	return sphere;
	else if(method=="rect")	return rect;
	else if(method=="zcw" || method=="zcw2d" || method=="zcw2D")	return zcw;
	else if(method=="sophie")	return sophie;
	else if(method=="liquid")	return liquid;
	else	return None;
}

powder::pTypes powder::canCalculate(std::string method)
{	return AssignM(method);	}

powder::powder(){
	__reset();
}

powder::powder(std::string method){
	//bool got=false;
	reset();
	_method=AssignM(method);
	_fname="NULL";
	if(_method==None){
		std::cout<<std::endl<<" powder::powder(method)"<<std::endl
			<<" method not found...assuming a file read in"<<std::endl
			<<std::endl;
		_method=None;
		_fname=method;
	}
}

powder::powder(pTypes method){
	//bool got=false;
	reset();
	_method=method;
	_fname="NULL";
	if(_method==None){
		std::cout<<std::endl<<" powder::powder(method)"<<std::endl
			<<" method not found...assuming a file read in"<<std::endl
			<<std::endl;
		_method=None;
		_fname=method;
	}
}

powder::powder(std::string method, int tstep, int pstep, int gstep){
	//bool got=false;
	reset();
	_method=AssignM(method);
	_fname="NULL";
	_tstep=tstep;
	_pstep=pstep;
	_gstep=gstep;

	if(_method==None){
		std::ifstream ftest(method.c_str());
		if(ftest.fail()){
			BLEXCEPTION("method not found...  OR the file was NOT found ")
		}else{
			_fname=method;
		}
		ftest.close();
	}
	calcPow();
}

powder::powder(pTypes method, int tstep, int pstep, int gstep){
	//bool got=false;
	reset();
	_method=method;
	_fname="NULL";
	_tstep=tstep;
	_pstep=pstep;
	_gstep=gstep;

	if(_method==None){
		BLEXCEPTION(" method not found... Using wrong constructor for a non file reading")
	}
	calcPow();
}


powder::powder(double oneth, double oneph, double oneg){
	reset();
	_theta.resize(0);
	_phi.resize(0);
	_weight.resize(0);
	_gamma.resize(0);
	_theta.push_back(oneth);
	_phi.push_back(oneph);
	_gamma.push_back(oneg);
	_method=None;
	_weight.push_back(1.);
}

powder::powder(double oneth){
	reset();
	_theta.push_back(oneth);
	_phi.push_back(oneth);
	_gamma.push_back(oneth);
	_method=None;
	_weight.push_back(1.);
}

powder::powder(int size){
	reset();
	_theta.resize(size, 0);
	_phi.resize(size, 0);
	_gamma.resize(size, 0);
	_method=None;
	_weight.resize(size, 1);
}


void powder::operator=(const powder &rhs){
	if(this==&rhs) return;
	_theta=rhs._theta;
	_phi=rhs._phi;
	_gamma=rhs._gamma;
	_weight=rhs._weight;
	_tstep=rhs._tstep;
	_pstep=rhs._pstep;
	_method=rhs._method;
	_fname=rhs._fname;
	_ct=rhs._ct;
	_underlimit=rhs._underlimit;
}

void powder::calcPow(){	calcPow(_method, _tstep, _pstep, _gstep); }

void powder::calcPow(powder::pTypes method, int tstep, int pstep, int gammast){

	if(method==None && _fname!="NULL"){
		read(_fname);
		return;
	}

	if(method==None) return;

	int k=0,m=0,i=0,j=0;
	_theta.resize(0);
	_phi.resize(0);
	_weight.resize(0);
	_gamma.resize(0);
	double theta=0.;
//	double phi=0.;
//vars for alderman
	double R=0.;
//vars for zcw
	int gott=0, gotp=0;
	const int vatstep[10]={21, 34, 55, 89, 144, 233, 377, 616, 987, 4184};
	const int vapstep[10]={8,  13, 21, 34,  55,  89, 144, 233, 377, 2584};
	double f1=0., f0=0., fac=0.;
	int gstep=gammast;
	if(gstep<=0) gstep=1;
	double gammaNOW=0;

	for(int gstepper=0; gstepper<gstep;++gstepper)
	{
		gammaNOW=double(gstepper)/double(_gstep)*PI;
		if(_gstep==0) gammaNOW=0;
		switch(method)
		{
			case sphere:
				for(k=1;k<=tstep;k++){
					for(m=1;m<=pstep;m++){
						theta=(pi/2.0)*double((k/double((tstep))));
						_theta.push_back(theta);
						_phi.push_back((2*pi)*double(((m)/double(pstep+1))));
						_weight.push_back(sin(theta));
						_gamma.push_back(gammaNOW);
					}
				}
				break;

			case rect:
				for(theta=pi/(2.*tstep); theta<=pi/2.; theta+=(pi/(2.*tstep))){
					double phi_step_max=2.*sin(theta)*_pstep;
					if(phi_step_max==0&&theta==0)		 phi_step_max=1;
					for(int q=0; q<phi_step_max; q++){
						_phi.push_back(double(2*pi*q)/double(phi_step_max));
						_theta.push_back(theta);
						_weight.push_back(sin(theta));
						_gamma.push_back(gammaNOW);
					}
				}
				break;

			case sophie:
				for( k=1;k<tstep;k++){
					for( m=1;m<=k;m++){
						_theta.push_back((pi/2)*double((k/double((tstep-1)))));
						_phi.push_back((pi/2)*double(((k-m+1)/double(k))));
						_weight.push_back(1.);
						_gamma.push_back(gammaNOW);
					}
				}
				break;

			case zcw:
				for(i=0;i<10;i++){
					if(tstep==vatstep[i]) 	gott=1;
					if(pstep==vapstep[i])	gotp=1;
				}

				if(gott==0 || gotp==0 || tstep<=pstep){
					std::cerr<<std::endl<<"Error: Powder.calPow()"<<std::endl;
					std::cerr<<" For using the ZCW method "<<std::endl;
					std::cerr<<" 'theta_step' and 'phi_step'  must be fibonacci numbers"<<std::endl;
					std::cerr<<" where 'theta_step'  > 'phi_step' "<<std::endl;
					std::cerr<<" valid values are"<<std::endl;
					std::cerr<<" theta_step: ";
					for(i=0;i<10;i++) std::cerr<<vatstep[i]<<" ";
					std::cerr<<std::endl<<" phi_step: ";
					for(i=0;i<10;i++) std::cerr<<vapstep[i]<<" ";
					std::cerr<<std::endl<<" A total of 'theta_step-1' are generated"<<std::endl;
					BLEXCEPTION("")
				}
				for(k=1;k<tstep;k++){

					fac=double(k)/double(tstep);
					f0=fac;
					f0=f0-floor(f0);
					f1=double(fac*pstep);
					f1 =f1-floor(f1);
					_phi.push_back(2.*pi*f0);
					_theta.push_back(pi*f1);
					double ww=1.0/double(tstep)*sin(pi*f1);
					_weight.push_back(ww);
					_gamma.push_back(gammaNOW);
				}
				break;

			case alderman:
					for(k=0;k<=1;k++){
					for( i=1; i<=tstep;i++){
						for( j=1;j<=tstep-i;j++){
							R=pow(double(i*i+j*j +(tstep-i-j)*(tstep-i-j)), double(0.5));
							_weight.push_back(1./pow(R, double(3.)));
							_theta.push_back(acos(double(tstep-i-j)/R));
							_phi.push_back(atan(double(j)/double(i)) +k*pi/2);
							_gamma.push_back(gammaNOW);
						}
					}
				}
				break;

			case liquid:
				_weight.push_back(1.);
				_theta.push_back(0.);
				_phi.push_back(0.);
				_gamma.push_back(0.);
				break;
			default:
				std::cout<<std::endl<<" powder::calcPow()"<<std::endl
					<<" a correct method was not found...nothing calculated"<<std::endl<<std::endl;

				break;
		}
	}
	//now normalize the weights to add to 1'
	_weight/=sum(_weight);
}

void powder::read(std::string file){
	std::ifstream infile(file.c_str());	//open the file for reading
  //make sure you've given me the correct file name
  if(infile.fail()){
	BLEXCEPTION(std::string(" no POWDER file name like ")+file)
  }
  _fname=std::string(file);
  read(infile);
}

void powder::read(const char *in){
	std::ifstream infile(in);	//open the file for reading
  //make sure you've given me the correct file name
  if(infile.fail()){
	BLEXCEPTION(std::string("no POWDER file name like ")+in)
  }
  _fname=std::string(in);
  read(infile);
}

void powder::read(std::ifstream &infile){

	Vector<double> blank;
	_theta=blank;
	_phi=blank;
	_weight=blank;
	char liner[500];   			//will be a std::string of each line in the file
  	int temp;
 		//make sure you've given me the correct file name
  	if(infile.fail()){
		BLEXCEPTION(std::string("no POWDER file \"")+_fname+std::string("\" found"))
  	}

	double tmwe=1.;		//temporary weight
	double tmth=0.;		//temp theta
	double tmph=0.;		//tmp phi
  //start reading the file
	bool gotth=false;	//lets me know if we had a valid theta
	int len=0;
	while((temp=infile.peek())!=EOF){
		infile.getline(liner,1000,'\n');  //get a line o' text

		gotth=false;
		Vector<std::string> tmm=parse_param(liner);
		len=tmm.size();
		if(len>0){
			if(tmm[0][0]!='#' && tmm[0][0]!='%'){
				if(len>=1){
					tmth=std::atof(tmm[0].c_str());
					_phi.push_back(tmth);
					gotth=true;
				}

				if(len>=2 && gotth){
					tmph=std::atof(tmm[1].c_str());
					_theta.push_back(tmph);
				}else if(gotth){
					_theta.push_back(0.);
				}

				if(len==3){
					tmwe=std::atof(tmm[2].c_str());
					_weight.push_back(tmwe);
					_gamma.push_back(0.0);
				}else if(len>=4){
					tmwe=std::atof(tmm[2].c_str());
					_gamma.push_back(tmwe);
					tmwe=std::atof(tmm[3].c_str());
					_weight.push_back(tmwe);
				}else if(gotth){
					_weight.push_back(1.0);
					_gamma.push_back(0.0);
				}
			}
		}
	}
}


END_BL_NAMESPACE


#endif

