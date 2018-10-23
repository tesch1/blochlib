


/* biot.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-28-01
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
	biot.cc--> calculates magnetic fields along grid point...implimnetation

	two classes here::
		1) the 'shape' generator--> generate the Coil shapes or reads them from a file
		2) the biot --> calculates the magnetic field along the grid points.....
*/

#ifndef _biot_cc_
#define _biot_cc_ 1


#include "bloch/biot.h"
#include "utils/matlab5.h"
#include "utils/utils.h"


BEGIN_BL_NAMESPACE


//
// this reads in the Parameter set "coil" from a file
// determins the shape, the number of points desired and
// caluclates the shape (or reads in a file)
//


void BiotCoil::operator=(const BiotCoil &rhs)
{
	if(&rhs==this) return;
	section=rhs.section;
	Coil=rhs.Coil;
	amps=rhs.amps;
	nloops=rhs.nloops;
	current=rhs.current;
	Min=rhs.Min;
	Max=rhs.Max;
	type_=rhs.type_;
	constField_=rhs.constField_;
	Symmetry=rhs.Symmetry;
	Direction=rhs.Direction;
	ShapeGenerate=rhs.ShapeGenerate;
	pset=rhs.pset;
}


void BiotCoil::read(){	read(pset);	}


void BiotCoil::read(std::string fname,std::string sec)
{
	pset.read(fname);
	section=sec;
	read(pset,sec);
}

void BiotCoil::read(const Vector<std::string> &pseti,std::string sec)
{
	Parameters pset(pseti);
	read(pset, sec);
}

void BiotCoil::read(Parameters &pseti,std::string sec)
{

	//section=sec;
	if(sec!=""){
		pseti.addSection(sec);
		pset=Parameters(pseti.section(sec));
	}else{
		pset=pseti;
	}

	section="";
	type_=pset.getParamS("type", section);
	if(type_=="constant"){
		constField_=pset.getParamCoordD("constField", section, ',');
	}else if(type_!="file"){
		amps=pset.getParamD("amps",section);
		nloops=pset.getParamD("loops",section);
		current=nloops*amps;

		/*std::string logg=pset.getParamS("logfile", "coil", false);
			if(logg=="" || logg=="std::cout" || logg=="stdout" || logg=="stderr"){
			}else{
				filebuf kl(logg.c_str(), std::ios::out);
				logfile=(kl);
				logfile<<std::endl<<"Begin:: BiotCoil logfile"<<std::endl;
		}*/

	/*	std::string sym=pset.getParamS("symmetry", sec, false, "");
		if(sym!=""){
			bool fw=true;
			if(sym=="C1"){	Symmetry=BiotSym::C1;	}
			else if(sym=="ignore"){	Symmetry=BiotSym::none;	}
			else if(sym=="none"){	Symmetry=BiotSym::none;	}
			else if(sym=="C2"){	Symmetry=BiotSym::C2;	}
			else if(sym=="C3"){	Symmetry=BiotSym::C3;	}
			else if(sym=="C4"){	Symmetry=BiotSym::C4;	}
			else if(sym=="Cinfh"){	Symmetry=BiotSym::Cinfh;	}
			else if(sym=="Dinfh"){	Symmetry=BiotSym::Dinfh;	}
			else if(sym=="D2h"){	Symmetry=BiotSym::D2h;	}
			else if(sym=="D4h"){	Symmetry=BiotSym::D4h;	}
			else{
				std::cerr<<std::endl<<"Warning: BiotCoil::read()"<<std::endl
						 <<" Symmetry '"<<sym<<"' is not valid only "<<std::endl
						 <<"   C1, C2, C3, C4, Cinfh, Dinfh, D2h, D4h" <<std::endl
						 <<" are allowed at this time..."<<std::endl
						 <<" setting value to 'C1' "<<std::endl;
				Symmetry=BiotSym::C1;
				fw=false;
			}

			if(fw && sym != "C1"){
				std::string dir=pset.getParamS("direction", sec, true); //direction is required
			}
		}else{
			Symmetry=BiotSym::C1;
		}*/
		calcShape();

	}else if(type_=="file"){
		fname=pset.getParamS("filename", section);
		readShape(fname);
		double tmamps=pset.getParamD("amps",section, false, 1);
		if(amps==0.0){ amps=tmamps;	}
		tmamps=pset.getParamD("loops",section, false, 1);
		if(nloops==0.0){ nloops=tmamps;	}
	}
	FindMinMax();
}



// ** AUX functions ** //
void BiotCoil::Translate(double x, double y, double z)
{	Translate(coord<>(x,y,z));	}

void BiotCoil::Translate(const coord<> &dir)
{
	for(int i=0;i<Coil.size();++i){
		for(int j=0;j<Coil[i].size();++j){	Coil[i][j]+=dir;	}
	}
	Min+=dir;
	Max+=dir;
}

coord<> BiotCoil::min() const
{	return Min;	}

coord<> BiotCoil::max() const
{	return Max;	}


void BiotCoil::FindMinMax()
{
	static double NonCalc=1e30;
	double xmax=-NonCalc, xmin=NonCalc;
	double ymax=-NonCalc, ymin=NonCalc;
	double zmax=-NonCalc, zmin=NonCalc;

	for(int i=0;i<Coil.size();++i){
		for(int j=0;j<Coil[i].size();++j){
			xmax=std::max(Coil[i][j].x(), xmax);
			ymax=std::max(Coil[i][j].y(), ymax);
			zmax=std::max(Coil[i][j].z(), zmax);

			xmin=std::min(Coil[i][j].x(), xmin);
			ymin=std::min(Coil[i][j].y(), ymin);
			zmin=std::min(Coil[i][j].z(), zmin);
		}
	}
	Min=coord<>(xmin, ymin, zmin);
	Max=coord<>(xmax, ymax, zmax);
}

void BiotCoil::calcShape()
{
	//genShapeCont_t::iterator myit;
	//myit=BiotFunctions.find(type_);
	//if(myit!=BiotFunctions.end()){
	//	ShapeGenerate=myit->second;
	//}else{	readShape(type_);	}
	ShapeGenerate=BiotFunctions.find(type_);
	if(ShapeGenerate) ShapeGenerate(pset, Coil);
	else readShape(type_);
}


//prints out the info...readable by this class
void BiotCoil::print(std::ostream &out)
{
	out<<section<<"{"<<std::endl;
	out<<std::endl<<"#BiotCoil.."<<std::endl;
	if(!pset.empty()){
		Vector<std::string> boitv=pset.section(section);
		for(int i=0;i<boitv.size();++i)
		{	out<<boitv(i)<<std::endl;	}
	}else{
		out<<"\ttype "<<type_<<std::endl;
		if(type_=="file") out<<"\tfilename "<<fname<<std::endl;
		out<<"\tamps "<<amps<<std::endl;
		out<<"\tloops "<<nloops<<std::endl;
	}
	out<<"}"<<std::endl;
}

void BiotCoil::readShape(std::string in) 
{
	std::ifstream oo(in.c_str(), std::ios::in);
	try{
		readShape(oo);
	}catch(BL_exception e){
		oo.close();
		e.print();
	}
}

void BiotCoil::readShape(std::ifstream & in) 
{
	if(in.fail()){	BLEXCEPTION("Shape file name cannot be read...")	}
	int counter=0;
	Coil.resize(1, Vector<coord<> >(MAXPOINTS,0.0));
	char liner[1000];
	Vector<std::string> tmm;
	int  master=0;
	int curp=in.tellg();
	while(!in.eof()){
		in.getline(liner,1000,'\n');
		if(std::string(liner)=="START BiotCoil")	break;
	}
	if(in.eof()){ in.seekg(curp); return;	}

	while(!in.eof()){
		in.getline(liner,1000,'\n');
		if(std::string(liner)=="END BiotCoil" || std::string(liner)=="END MultiBiotCoil")	break;
		tmm=parse_param(std::string(liner));
		if(counter>MAXPOINTS){
			std::cout<<std::endl<<"Error: BiotCoil::readShape(std::string)"<<std::endl;
			std::cout<<" Too many points in the file...stopping at "<<MAXPOINTS<<std::endl;
			return;
		}
		if(tmm.size()>0){
			if(tmm[0]=="section" && counter!=0){
				Coil[master].resizeAndPreserve(counter);
				Coil.push_back(Vector<coord<> >(MAXPOINTS,0.0));
				counter=0;
				master++;
			}
			if(tmm[0]=="amps" || tmm[0]=="amps:") amps=std::atof(tmm[1].c_str());
			if(tmm[0]=="loops" || tmm[0]=="loops") nloops=std::atof(tmm[1].c_str());

			if(tmm.size()>=4){
				int ll=int(std::atof(tmm[0].c_str()))-1;
				Coil[master][ll][0]=std::atof(tmm[1].c_str());
				Coil[master][ll][1]=std::atof(tmm[2].c_str());
				Coil[master][ll][2]=std::atof(tmm[3].c_str());
				counter++;
			}else if(tmm.size()>=3){
				Coil[master][counter][0]=std::atof(tmm[0].c_str());
				Coil[master][counter][1]=std::atof(tmm[1].c_str());
				Coil[master][counter][2]=std::atof(tmm[2].c_str());
				//Coil.push_back(coord<>(std::atof(tmm[0].c_str()),std::atof(tmm[1].c_str()),std::atof(tmm[2].c_str())));
				counter++;
			}
		}
		//else{
		//	std::cout<<std::endl<<"Warning: BiotCoil::readShape(std::string)"<<std::endl;
		//	std::cout<<" line "<<temp<<" in file "<<pfname<<" does not have 3 or 4 columns...it will ignored"<<std::endl;
		//}
	}
	if(counter==0){
		std::cout<<std::endl<<"Warning: BiotCoil::readShape(std::string)"<<std::endl;
		std::cout<<" nothing could be read from file "<<std::endl;
	}
	Coil[master].resizeAndPreserve(counter);

	numpts=counter;
	FindMinMax();
}

/*** SHAPE WRITERS **/
void BiotCoil::writeShape(std::string &out) 
{
	std::fstream oo(out.c_str(), std::ios::out);
	try{
		writeShape(oo);
	}catch(BL_exception e){
		oo.close();
		e.print();
	}
}

void BiotCoil::writeShape(std::fstream &out) 
{
	if(out.fail()){ BLEXCEPTION("out file cannot be opened.")	}
	
	out<<std::endl<<"START BiotCoil"<<std::endl;
	out<<"amps "<<amps<<std::endl;
	out<<"loops "<<nloops<<std::endl;
	for(int i=0;i<Coil.size();++i){
		out<<std::endl<<"section"<<std::endl;
		for(int j=0;j<Coil[i].size();++j){
			out<<Coil[i][j]<<std::endl;
		}
	}
	out<<std::endl<<"END BiotCoil"<<std::endl;
}


/********************************************************************/
/***  BIOTSYMBIT CLASS ***/
/********************************************************************/

std::ostream &operator<<(std::ostream &out, const BiotSymBit &oo)
{
	out<<"Index: "<<oo.index<<std::endl;
	out<<"Sym Points: "<<oo.symindex<<std::endl;
	out<<"Signs :"<<oo.sign<<std::endl;
	return out;
}



/********************************************************************/
/***  BIOTSYMCALC CLASS ***/
/********************************************************************/


BiotSymBit &BiotSymCalc::operator[](int i)
{	return syms(i);	}

BiotSymBit &BiotSymCalc::operator()(int i)
{	return syms(i);	}

BiotSymBit BiotSymCalc::operator[](int i) const
{	return syms(i);	}

BiotSymBit BiotSymCalc::operator()(int i) const
{	return syms(i);	}

END_BL_NAMESPACE


#endif
