

/* biot_meth.h ********/


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
	biot_meth.h--> calculates magnetic fields along grid point...implimnetation

	the templated functions
*/

#ifndef _biot_meth_h_
#define _biot_meth_h_ 1


//#include "bloch/biot.h"
#include "utils/matlab5.h"
#include "mpi/mpi_config.h"


BEGIN_BL_NAMESPACE




/********************************************************************/
/***  BIOTSYMCALC CLASS ***/
/********************************************************************/



template<class Grid_t>
void BiotSymCalc::calcSym(const Grid_t &grid, const BiotCoil &coil)
{
	//fill the syms with the total indexes...for 'C1' symmtry
	if(coil.Symmetry == BiotSym::Dinfh){
		calcDinfh(grid, coil);
		return;
	}

	//Last possible Symmetry..i.e. NO symmetry
	Grid_t &gridr=const_cast<Grid_t &>(grid);
	syms.resize(grid.size());
	typename Grid_t::iterator myit(gridr);
	while(myit){	syms(myit.curpos()).index=myit.curpos(); ++myit;	}
	return;
}


//Dinfh Symmetry
// becuase we do not know if there is any grid points along the mirror axis
// we must chop the grid into octants
// for this we use the GRID not the shape...we find the large octant
// in the grid from the coils center point
template<class Grid_t>
void BiotSymCalc::calcDinfh(const Grid_t &grid, const BiotCoil &coil)
{
		//no saving time if less then 8 grid points
		Grid_t &gridr=const_cast<Grid_t &>(grid);
		if(grid.size()<8){
			syms.resize(grid.size());
			typename Grid_t::iterator myit(gridr);
			while(myit){	syms(myit.curpos()).index=myit.curpos(); ++myit;	}
			return;
		}

		coord<> center=((coil.max()-coil.min())/2.0)+coil.min();
		//Grid<UniformGrid> gridcp(grid);

		//first see if the center point of the coil is also the center point
		//of the grid (then we have equal sized octants)
		coord<> gcenter=((gridr.Max()-gridr.Min())/2.0)+gridr.Max();
		double Cutoff=1e-10;
		if(square_norm(gcenter-center)<Cutoff){
			if(gridr.size()%2 !=0) { //odd count include the points on the line
				int dims=int(ceil(double(gridr.size())/double(8)));
				syms.resize(dims);
			}else{
				syms.resize(gridr.size()/8);
			}

			coord<> rectMin(center), rectMax(gridr.Max());
			if(center.x() > gridr.Max().x() || center.y() > gridr.Max().y()||center.z() > gridr.Max().z()){
				rectMin=gridr.Max();
				rectMax=center;
			}
			XYZrect rect(rectMin, rectMax); //take the top most octant
			typename Grid_t::iterator myit(gridr);
			int ct=0;
			while(myit){ //get all the points in the octant
				if(rect.ShapeFunc(myit.Point())){	syms[ct].index=myit.curpos(); ++ct;	}
				++myit;
			}

			//now get the symmetry index points
			coord<> altx, alty, altz, altxy, altxz, altyz, altxyz;
			coord<> curpt;
			for(int i=0;i<ct;++i){
				curpt=gridr.Point(syms[i].index);
				typename Grid_t::iterator myitinner(gridr);
				altx(-curpt.x(), curpt.y(), curpt.z());
				alty( curpt.x(),-curpt.y(), curpt.z());
				altz( curpt.x(), curpt.y(),-curpt.z());
				altxy(-curpt.x(),-curpt.y(), curpt.z());
				altxz(-curpt.x(), curpt.y(),-curpt.z());
				altyz( curpt.x(),-curpt.y(),-curpt.z());
				altxyz( -curpt.x(),-curpt.y(),-curpt.z());

				//the sign changes based on the axis or rotation...it defaults to the Z axis
				// The sign of the field changes in the direction of the 'alter' coordinate
				while(myitinner){
					switch(coil.Direction){
				// if the rotation axis is X then
				//  an alteration of -x changes the sign of the By and Bz...
				//  an alteration of -y changes the sign of the By
				//  an alteration of -z changes sign of Bz
					case BiotSym::X:
						if(square_norm(altx-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,-1,-1));
						}else if(square_norm(alty-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,-1,1));
						}else if(square_norm(altxy-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,1,-1));
						}else if(square_norm(altz-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,1,-1));
						}else if(square_norm(altxz-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,-1,1));
						}else if(square_norm(altyz-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,-1,-1));
						}else if(square_norm(altxyz-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,1,1));
						}
						break;
				// if the rotation axis is Y then
				//  an alteration of -x changes the sign of the Bx...
				//  an alteration of -y changes the sign of the Bz and Bx
				//  an alteration of -z changes sign of Bz
					case BiotSym::Y:
						if(square_norm(altx-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(-1,1,1));
						}else if(square_norm(alty-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(-1,1,-1));
						}else if(square_norm(altxy-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,1,-1));
						}else if(square_norm(altz-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,1,-1));
						}else if(square_norm(altxz-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(-1,1,-1));
						}else if(square_norm(altyz-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(-1,1,1));
						}else if(square_norm(altxyz-myitinner.Point())<Cutoff){
							syms[i].symindex.push_back(myitinner.curpos());
							syms[i].sign.push_back(coord<>(1,1,1));
						}
						break;

				// if the rotation axis is Z then
				//  an alteration of -x changes the sign of the Bx...
				//  an alteration of -y changes the sign of the By
				//  an alteration of -z does Bx and By change sign
						case BiotSym::Z:
						default:
							if(square_norm(altx-myitinner.Point())<Cutoff){
								syms[i].symindex.push_back(myitinner.curpos());
								syms[i].sign.push_back(coord<>(-1,1,1));
							}else if(square_norm(alty-myitinner.Point())<Cutoff){
								syms[i].symindex.push_back(myitinner.curpos());
								syms[i].sign.push_back(coord<>(1,-1,1));
							}else if(square_norm(altxy-myitinner.Point())<Cutoff){
								syms[i].symindex.push_back(myitinner.curpos());
								syms[i].sign.push_back(coord<>(-1,-1,1));
							}else if(square_norm(altz-myitinner.Point())<Cutoff){
								syms[i].symindex.push_back(myitinner.curpos());
								syms[i].sign.push_back(coord<>(-1,-1,1));
							}else if(square_norm(altxz-myitinner.Point())<Cutoff){
								syms[i].symindex.push_back(myitinner.curpos());
								syms[i].sign.push_back(coord<>(1,-1,1));
							}else if(square_norm(altyz-myitinner.Point())<Cutoff){
								syms[i].symindex.push_back(myitinner.curpos());
								syms[i].sign.push_back(coord<>(-1,1,1));
							}else if(square_norm(altxyz-myitinner.Point())<Cutoff){
								syms[i].symindex.push_back(myitinner.curpos());
								syms[i].sign.push_back(coord<>(1,1,1));
							}

							break;
					}
					++myitinner;
				}
			}
		}

/*

		//the 8 octants..
		coord<> oct1=gridcp.Max(); double v1=prod(abs(oct1- center));
		coord<> oct2=gridcp.Min(); double v2=prod(abs(oct2- center));
		coord<> oct3(oct1.x(), oct1.y(), oct2.z()); double v3=prod(abs(oct3- center));
		coord<> oct4(oct1.x(), oct2.y(), oct1.z()); double v4=prod(abs(oct4- center));
		coord<> oct5(oct2.x(), oct1.y(), oct1.z()); double v5=prod(abs(oct5- center));
		coord<> oct6(oct2.x(), oct2.y(), oct1.z()); double v6=prod(abs(oct6- center));
		coord<> oct7(oct2.x(), oct1.y(), oct2.z()); double v7=prod(abs(oct7- center));
		coord<> oct8(oct1.x(), oct2.y(), oct2.z()); double v8=prod(abs(oct8- center));
*/
}

/********************************************************************/
/***  BIOT CLASS ***/
/********************************************************************/

template<class Grid_t>
void Biot<Grid_t>::read(std::string &in, bool readgrid)
{
	BiotCoil::read(in);
	if(readgrid) pset.addSection("grid");
	Biot::read(pset);
}

template<class Grid_t>
void Biot<Grid_t>::read(Parameters &in, bool readgrid)
{
	if(readgrid)
	{
		pset.addSection("grid");
		coord<> mymax=pset.getParamCoordD("max", "grid");
		coord<> mymin=pset.getParamCoordD("min", "grid");
		coord<int> mydim=pset.getParamCoordI("dim", "grid");

		Grid<UniformGrid> tmps(mymin, mymax, mydim);
		grid=XYZshape<XYZfull>(tmps, XYZfull());
	}
	symcalc.calcSym(tmps, *this);
	Bfield.resize(grid.size());
}

template<class Grid_t>
void Biot<Grid_t>::read(){	BiotCoil::read(); read(BiotCoil::pset, true);	}

template<class Grid_t>
bool Biot<Grid_t>::calculateField()
{

	//Grid<UniformGrid>::iterator myit(grid);
	if(BiotCoil::type()=="constant"){
		Bfield.resize(grid.size());
		Bfield.fill(BiotCoil::constField());
		return 	true;
	}
	current=nloops*amps;
	if(current==0.0){
		Bfield.resize(grid.size(), 0.0);
	 	return true;
	}
	Bfield.resize(grid.size());
	Bfield.fill(0.0);
	//double dist=0.0, Cutoff=1.e-3;
	int i=0, j=0, ct=0;
	int begin=0, end=symcalc.size(),  done=-1;
//first make sure we are Parrellel
	if(Controller.parallel()){
	//Perform the "master/slave" loops
		if(Controller.master()){
			int *procIdx; procIdx=new int[Controller.size()];

		//First we see if we send an initial request to each
			for(int rank=1;rank<Controller.size();++rank){
				procIdx[rank]=ct;
				Controller.put(ct, rank); ++ct;
				if(ct>=end) break;
			}

		//keep sending to any proc that is willing to get one
			coord<> curBf=0.0;
			while(ct<end){
				int get=Controller.getAny(curBf); //get the calculated Bfield point
				Bfield(symcalc[procIdx[get]].index)+=curBf;//put it in the correct place
				procIdx[get]=ct;
				Controller.put(ct,get); //send the proc a new point to do
				++ct;
			}

		//get the last returns
			for(int qq=1;qq<Controller.size();++qq){
				Controller.get(curBf,qq);
				Bfield(symcalc[procIdx[qq]].index)+=curBf;//put it in the correct place
			}

		//put the final kills
			for(int qq=1;qq<Controller.size();++qq)
				Controller.put(done, qq);

			delete [] procIdx;

	//the slave procs...
		}else{
			coord<> curB=0.0;
			int cur;
			while(1)
			{
				Controller.get(cur,0);
				if(cur==-1) break;
				curB=calculateField(cur); //calculate the field at that point
				Controller.put(curB,0);
			}
		}

//here is the serial Implimentation
	}else{
		for(int ct=begin;ct<end;++ct){
			Bfield(symcalc[ct].index)+=calculateField(ct);
		}
	}

//fix up the 'rest' of the Bfield for symmetry sake
	if(symcalc.size() != grid.size())
	{
		for(i=begin;i<end;++i){
			for(j=0;j<symcalc[i].symindex.size();++j){
				Bfield(symcalc[i].symindex[j])=symcalc[i].sign[j]*Bfield(symcalc[i].index);
			}
		}
	}

//finally let everyone have the total Bfield
	Controller.scatter(Bfield);
	return true;
}




//calculates the field at the index 'ct'
template<class Grid_t>
coord<> Biot<Grid_t>::calculateField(int ct)
{
	if(type() == "constant"){
		return constField();
	}
	int i,j;
	coord<> r;
	double factor(0.0);
	static double Cutoff=1.e-3;
	coord<> curB=0.0;
	for(i=0;i<Coil.size();++i)
	{
		Vector<coord<> > ds(Coil[i].size()+1,0.0);
		Range I(1,Coil[i].size()), J(I-1);
		ds(J)=Coil[i](I)-Coil[i](J);

		ds[ds.size()-1]=ds[ds.size()-2];

		for(j=0;j<Coil[i].size()-1; ++j){
			//distance from current point to the coil point
			r=(Coil[i][j]+Coil[i][j+1])/2.0-grid.Point(symcalc[ct].index);

			//constant factor
			factor=cube(norm(r));
			if(factor>Cutoff){
				factor=current*1.0e-7/(factor*1.0e-6);
				curB+=factor*cross(ds(j),r);
			}else{
				curB=10.0;
			}
		}
	}
	return curB;
}

template<class Grid_t>
void Biot<Grid_t>::writeMatlab(std::string fname) 
{
	matstream matout(fname, std::ios::out| std::ios::binary);
	try{
		writeMatlab(matout);
	}catch(BL_exception e){
		matout.close();
		e.print();
	}
}

template<class Grid_t>
void Biot<Grid_t>::writeMatlab(matstream &matout)  
{

	if(matout.fail()){ BLEXCEPTION("out file cannot be opened.")	}

	matout.put("B", Bfield);
	Vector<coord<> > Points(BiotCoil::size(), 0.0);
	Vector<double> pdiv;
	int ct=0;
	if(Coil.type()!="constant"){
		for(int i=0;i<Coil.size();++i){
			pdiv.push_back(ct);
			for(int j=0;j<Coil[i].size();++j){
				Points[ct]=Coil[i][j];
				ct++;
			}
		}
	}
	pdiv.push_back(ct);
	matout.put("pdiv", pdiv);
	matout.put("P", Points);
	matout.put("grid", grid.data());
	matout.put("N", grid.dim());
	coord<> tmp=grid.Max()-grid.Min();
	matout.put("r", tmp);
	matout.put("res",grid.dr());
	matout.close();
}

//writes a text file with both the
//Shape and the Bfield
template<class Grid_t>
void Biot<Grid_t>::write(std::string fname) 
{
	std::fstream matout(fname, std::ios::out| std::ios::binary);
	try{
		write(matout);
	}catch(BL_exception e){
		matout.close();
		e.print();
	}
}

template<class Grid_t>
void Biot<Grid_t>::write(std::fstream &out) 
{
	if(out.fail()){ BLEXCEPTION("out file cannot be opened.")	}
	writeShape(out);
	out<<std::endl<<"START BiotCoilBfield"<<std::endl;
	out<<"size "<<Bfield.size()<<std::endl;
	for(int i=0;i<Bfield.size();++i){
		out<<Bfield[i]<<std::endl;
	}
	out<<std::endl<<"END BiotCoilBfield"<<std::endl;
}

//reads a text file with both the
//Shape and the Bfield
template<class Grid_t>
void Biot<Grid_t>::read(std::ifstream &in) 
{
	if(in.fail()){	BLEXCEPTION(" Shape file name cannot be read...")	}
	
	BiotCoil::readShape(in);
	int counter=0;
	char liner[1000];
	Vector<std::string> tmm;

	int curp=in.tellg();
	while(!in.eof()){
		in.getline(liner,1000,'\n');
		if(std::string(liner)=="START BiotCoilBfield")	break;
	}
	if(in.eof()){ in.seekg(curp); return;	}

	while(!in.eof()){
		in.getline(liner,1000,'\n');
		if(std::string(liner)=="END BiotCoilBfield" || std::string(liner)=="END MultiBiotCoil")	break;
		tmm=parse_param(std::string(liner));


		if(tmm.size()>0)
			if(tmm[0]=="size")
				Bfield.resize(std::atoi(tmm[1].c_str()),0.0);

		if(tmm.size()>=3){
			if(counter>Bfield.size()){
				std::cout<<std::endl<<"Error: Biot<Grid>::read(std::fstream)"<<std::endl;
				std::cout<<" Too many Magnetic Field points in the file...stopping at "<<Bfield.size()<<std::endl;
				return;
			}
			Bfield[counter][0]=std::atof(tmm[0].c_str());
			Bfield[counter][1]=std::atof(tmm[1].c_str());
			Bfield[counter][2]=std::atof(tmm[2].c_str());
			counter++;
		}
	}
	if(counter==0){
		std::cout<<std::endl<<"Warning: Biot<Grid>::read(std::fstream)"<<std::endl;
		std::cout<<" nothing could be read from file "<<std::endl;
	}
}


/********************************************************************/
/*** Multi  BIOT CLASS ***/
/********************************************************************/

template<class Grid_t>
MultiBiot<Grid_t>::MultiBiot(Parameters &pset, std::string sec)
{	read(pset,sec);	}

template<class Grid_t>
MultiBiot<Grid_t>::MultiBiot(Grid_t &ingrid, Parameters &pset, std::string sec):
	grid(ingrid)
{	read(pset,sec, false);	}

template<class Grid_t>
void MultiBiot<Grid_t>::read(Parameters &pset,std::string sec, bool readgrid)
{
	Parameters subp;
	if(sec!=""){
		pset.addSection(sec);
		subp=Parameters(pset.section(sec));
	}else{
		subp=pset;
	}
//signifies the entre MULTI BIOT is in one file
	std::string toF=pset.getParamS("type", sec, false, "");
	if(toF=="file"){
		std::string loo=pset.getParamS("filename", sec);
		std::ifstream fin(loo.c_str(), std::ios::in);
		read(fin);
		fin.close();
	}else{
		std::string subbase=subp.getParamS("base","",false, "subcoil");
		int maxFit=subp.getParamI("numpars","",false, 1000000);
		int numCoils=0;
	//count the number of params present
		while(subp.addSection(subbase+itost(numCoils+1)) && numCoils<=maxFit )
		{	numCoils++;		}
	//add a parameter to our master list
		Coils.resize(numCoils);
		int i=0;
		while(i<Coils.size())
		{
			Coils[i]=BiotCoil(subp, subbase+itost(i+1));
			++i;
		}
	}
	if(readgrid){
		pset.addSection("grid");
		grid=Grid_t(Grid<UniformGrid>(pset.getParamCoordD("min", "grid"),
					pset.getParamCoordD("max", "grid"),
					pset.getParamCoordI("dim", "grid")));
	}
}

template<class Grid_t>
int MultiBiot<Grid_t>::CoilSize() const
{	int ll=0; for(int i=0;i<Coils.size();++i){	ll+=Coils[i].size();	}	return ll;	}

template<class Grid_t>
void MultiBiot<Grid_t>::calculateField()
{
	Bfield.resize(grid.size());
	Bfield.fill(0.0);
	for(int i=0;i<Coils.size();++i){
		Biot<Grid_t> tmB(Coils[i], grid);
		tmB.Controller=Controller;
		tmB.calculateField();
		Bfield+=tmB.Bfield;
	}
}

template<class Grid_t>
coord<> MultiBiot<Grid_t>::MaxField(coord<int> rotframe){
	double sigz=1.0, sigx=1.0, sigy=1.0,
	Bz=-1e30, By=Bz, Bx=Bz,
	oldBz=Bz, oldBy=By, oldBx=Bx;

	for(int i=0;i<Bfield.size();++i){
		if(rotframe.z()){
			Bz=max(Bz, abs(Bfield(i).z()));
			if(oldBz!=Bz){ sigz=sign(Bfield(i).z()); oldBz=Bz;	}
		}else{
			Bz=0.0;
		}
		if(rotframe.x()){
			Bx=max(Bx, abs(Bfield(i).x()));
			if(oldBx!=Bx){ sigx=sign(Bfield(i).x()); oldBx=Bx;	}
		}else{
			Bx=0.0;
		}
		if(rotframe.y()){
			By=max(By, abs(Bfield(i).y()));
			if(oldBy!=By){ sigy=sign(Bfield(i).y()); oldBy=By;	}
		}else{
			By=0.0;
		}
	}
	return coord<>(sigx*Bx, sigy*By, sigz*Bz);
}

template<class Grid_t>
coord<> MultiBiot<Grid_t>::AverageField(coord<int> rotframe)
{
	coord<> suml=sum(Bfield);

	if(!rotframe.x()) suml.x()=0;
	if(!rotframe.y()) suml.y()=0;
	if(!rotframe.z()) suml.z()=0;

	return suml/Bfield.size();
}

template<class Grid_t>
void MultiBiot<Grid_t>::writeMatlab(std::string fname)
{
	matstream matout(fname, std::ios::out| std::ios::binary);
	writeMatlab(matout);
}

template<class Grid_t>
void MultiBiot<Grid_t>::writeMatlab(matstream &matout)
{
	matout.put("B", Bfield);
	Vector<coord<> > Points(CoilSize(), 0.0);
	Vector<double> pdiv;
	int ct=0;
	for(int j=0;j<Coils.size();++j){
		if(Coils[j].type()!="constant"){
			for(int i=0;i<Coils[j].Coil.size();++i){
				pdiv.push_back(ct);
				for(int k=0;k<Coils[j].Coil[i].size();++k){
					Points[ct]=Coils[j].Coil[i][k];
					ct++;
				}
			}
		}
	}
	pdiv.push_back(ct);
	matout.put("pdiv", pdiv);
	matout.put("P", Points);
	matout.put("grid", grid.data());
	matout.put("N", grid.dim());
	coord<> tmp=grid.Max()-grid.Min();
	matout.put("r", tmp);
	matout.put("res",grid.dr());
	matout.close();
}

//reads coil points from files...
template<class Grid_t>
void MultiBiot<Grid_t>::readShape(std::string fname) 
{
	std::ifstream ii(fname.c_str(), std::ios::in);
	try{
		readShape(ii);
	}catch(BL_exception e){
		ii.close();
		e.print();
	}
}

template<class Grid_t>
void MultiBiot<Grid_t>::readShape(std::ifstream &in) 
{
	if(in.fail()){ BLEXCEPTION(" Shape File cannot be read ") }
	char liner[1000];
	int curp=in.tellg();
	while(!in.eof()){
		in.getline(liner,1000,'\n');
		if(std::string(liner)=="START MultiBiotCoil")
			break;
	}
	if(in.eof()){ in.seekg(curp); return;	}
	BiotCoil tmpC;
	while(!in.eof()){
		curp=in.tellg();
		in.getline(liner,1000,'\n');
		//cout<<liner<<endl;
		if(std::string(liner)=="END MultiBiotCoil")	break;
		if(std::string(liner)=="START BiotCoil"){
			in.seekg(curp);
			tmpC.readShape(in);
			//tmpC.Coil.print(cout, "\n");
			Coils.push_back(tmpC);
		}
	}
}

//writes coil points from files...
template<class Grid_t>
void MultiBiot<Grid_t>::writeShape(std::string out)
{
	std::fstream oo(out.c_str(), std::ios::out);
	writeShape(oo);
}

template<class Grid_t>
void MultiBiot<Grid_t>::writeShape(std::fstream &oo)	
{
	if(oo.fail()){  BLEXCEPTION(" Shape File cannot be read ")	}
	
	oo<<std::endl<<"START MultiBiotCoil"<<std::endl;
	for(int i=0;i<Coils.size();++i){
		Coils[i].writeShape(oo);
	}
	oo<<std::endl<<"END MultiBiotCoil"<<std::endl;
}

//writes a text file with both the
//Shape and the Bfield
template<class Grid_t>
void MultiBiot<Grid_t>::write(std::string fname) 
{
	std::fstream matout(fname.c_str(), std::ios::out| std::ios::binary);
	try{
		write(matout);
	}catch(BL_exception e){
		matout.close();
		e.print();
	}
}

template<class Grid_t>
void MultiBiot<Grid_t>::write(std::fstream &oo)	
{
	if(oo.fail()){ BLEXCEPTION(" Shape File cannot be written to ")	}
	
	writeShape(oo);
	oo<<std::endl<<"START MultiBiotCoilBfield"<<std::endl;
	oo<<"size "<<Bfield.size()<<std::endl;
	for(int i=0;i<Bfield.size();++i){
		oo<<Bfield[i]<<std::endl;
	}
	oo<<std::endl<<"END MultiBiotCoilBfield"<<std::endl;

}

//reads a text file with both the
//Shape and the Bfield
template<class Grid_t>
void MultiBiot<Grid_t>::read(std::string in) 
{
	std::ifstream matout(fname, std::ios::out| std::ios::binary);
	try{
		read(matout);
	}catch(BL_exception e){
		matout.close();
		e.print();
	}
}

template<class Grid_t>
void MultiBiot<Grid_t>::read(std::ifstream &in) 
{
	if(in.fail()){    BLEXCEPTION(" Shape file name cannot be read...")	}
	readShape(in);
	int counter=0;
	char liner[1000];
	Vector<std::string> tmm;

	int curp=in.tellg();
	while(!in.eof()){
		in.getline(liner,1000,'\n');
		if(std::string(liner)=="START MultiBiotCoilBfield")	break;
	}
	if(in.eof()){ in.seekg(curp); return;	}

	while(!in.eof()){
		in.getline(liner,1000,'\n');
		if(std::string(liner)=="END MultiBiotCoilBfield")	break;
		tmm=parse_param(std::string(liner));


		if(tmm.size()>0)
			if(tmm[0]=="size")
				Bfield.resize(std::atoi(tmm[1].c_str()),0.0);

		if(tmm.size()>=3){
			if(counter>Bfield.size()){
				std::cout<<std::endl<<"Error: MultiBiot<Grid>::read(std::ifstream)"<<std::endl;
				std::cout<<" Too many Magnetic Field points in the file...stopping at "<<Bfield.size()<<std::endl;
				return;
			}
			Bfield[counter][0]=std::atof(tmm[0].c_str());
			Bfield[counter][1]=std::atof(tmm[1].c_str());
			Bfield[counter][2]=std::atof(tmm[2].c_str());
			counter++;
		}
	}
	if(counter==0){
		std::cout<<std::endl<<"Warning: MultiBiot<Grid>::read(std::fstream)"<<std::endl;
		std::cout<<" nothing could be read from file "<<std::endl;
	}
}

END_BL_NAMESPACE

#endif
