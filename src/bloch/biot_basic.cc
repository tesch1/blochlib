

/* biot_basic.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-20-01
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
	biot_basic.cc--> some basic shape functins and
	a function registration class to be used inside
	BiotCoil...

*/

#ifndef _Biot_Basic_cc_
#define _Biot_Basic_cc_ 1

#include "utils/params.h"
#include "utils/matlab5.h"
#include "bloch/biot_basic.h"
#include <map>

BEGIN_BL_NAMESPACE

/************** Several 'basic' Shape functions **************/

//Our Master Biot Function List...
//map<std::string, void *>  BiotFunctions;
//BiotFuncMap<void *> BiotFunctions;

//void registerBiotShape( std::string name,  void (*genShape_t)(Parameters &pset, Vector<Vector<coord<> > > & Coil))
//{
//	BiotFunctions.insert(pair<std::string, void *>(name, (void (*)(void))genShape_t);
//}

//BiotFuncMap<genShape_t> BiotFunctions;

BiotFuncMap::BiotFuncMap()
{
	insert("circle", Biot_circle);
	insert("line", Biot_line);
	insert("helix", Biot_helix);
	insert("spiral", Biot_spiral);
	insert("helmholtz", Biot_helmholtz);
	insert("truehelmholtz", Biot_helmholtz_true);
}
//registerBiotShape("circle", Biot_circle);
/*registerBiotShape("line", Biot_line);
registerBiotShape("helix", Biot_helix);
registerBiotShape("spiral", Biot_spiral);
registerBiotShape("helmholtz", Biot_helmholtz);
registerBiotShape("truehelmholtz", Biot_helmholtz_true);
registerBiotShape("saddlecoil", Biot_saddle_coil);
*/

void BiotFuncMap::insert(std::string name, genShape_t in)
{	mydata.insert(std::pair<std::string, genShape_t>(name, in));	}

BiotFuncMap::Func_T BiotFuncMap::find(std::string type_)
{
	myMap_t::iterator myit;
	myit=mydata.find(type_);
	if(myit!=mydata.end()){
		return myit->second;
	}
	return 0;
}

void BiotFuncMap::print(std::ostream &oo)
{
	myMap_t::iterator myit;
 	for(myit = mydata.begin(); myit != mydata.end(); myit++)
    {
        oo<< "Register Biot Function: "<<myit->first<<std::endl;
    }
}

/*** DECLARE THE GLOBAL BIOT FUNCTIONS CONTAINER ***/
BiotFuncMap BiotFunctions;

void Biot_circle(Parameters &pset, Vector<Vector<coord<> > > &Coil)
{
	int i;
	std::string section="";
	double R=pset.getParamD("R", section);
	int numpts=pset.getParamI("numpts", section);
	char axis=pset.getParamC("axis", section, false, 'z');
	coord<> center=pset.getParamCoordD("center", section, ',',false);
	Coil.resize(1, Vector<coord<> >(numpts, 0.0));
	double angle=2.0*PI/double(numpts-1);

	double rotangle=pset.getParamC("rotangle", section, false, 0)*DEG2RAD;
	coord<> rotaxis(0,0,1);
	if(rotangle!=0){
		rotaxis=pset.getParamCoordD("rotaxis", section);
		rotaxis/=norm(rotaxis);
	}

	switch(axis){
		case 'y':
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=R*cos(i*angle);
				Coil[0][i][1]=0.0;
				Coil[0][i][2]=R*sin(i*angle);
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
		case 'x':
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=0.0;
				Coil[0][i][1]=R*cos(i*angle);
				Coil[0][i][2]=R*sin(i*angle);
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
		case 'z':
		default:
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=R*cos(i*angle);
				Coil[0][i][1]=R*sin(i*angle);
				Coil[0][i][2]=0.0;
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
	}

	std::cout<<"Circle: Radius "<<R<<" center: ["<<center<<"]"<<std::endl;
}

void Biot_line(Parameters &pset,Vector<Vector<coord<> > > &Coil)
{
	int i;
	std::string section="";
	coord<> begin=pset.getParamCoordD("begin", section);
	int numpts=pset.getParamI("numpts", section);
	coord<> end=pset.getParamCoordD("end", section);
	Coil.resize(1, Vector<coord<> >(numpts, 0.0));
	coord<> div=(end-begin)/double(numpts-1);
	coord<> ston=begin;
	for (i=0;i<numpts;i++){
		Coil[0][i][0]=ston.x();
		Coil[0][i][1]=ston.y();
		Coil[0][i][2]=ston.z();
		ston+=div;
	}
	std::cout<<"Line from ["<<begin<<"] to ["<<end<<"]"<<std::endl;
}

void Biot_helix(Parameters &pset, Vector<Vector<coord<> > > &Coil)
{
	int i;
	std::string section="";
	double R=pset.getParamD("R", section);
	double dist=pset.getParamD("Z", section);
	int turns=pset.getParamI("turns", section);
	int numpts=pset.getParamI("numpts", section);
	char axis=pset.getParamC("axis", section, false, 'z');
	coord<> center=pset.getParamCoordD("center", section, ',',false);
	double angle=2.0*PI*double(turns)/(numpts-1);
	Coil.resize(1, Vector<coord<> >(numpts,0.0));

	double rotangle=pset.getParamC("rotangle", section, false, 0)*DEG2RAD;
	coord<> rotaxis(0,0,1);
	if(rotangle!=0){
		rotaxis=pset.getParamCoordD("rotaxis", section);
		rotaxis/=norm(rotaxis);
	}

	switch(axis){
		case 'x':
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=i*double(dist)/double(numpts)-dist/2.0;
				Coil[0][i][1]=R*sin(i*angle);
				Coil[0][i][2]=R*cos(i*angle);
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
		case 'y':
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=R*cos(i*angle);
				Coil[0][i][1]=i*double(dist)/double(numpts)-dist/2.0;
				Coil[0][i][2]=R*sin(i*angle);
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
		case 'z':
		default:
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=R*cos(i*angle);
				Coil[0][i][1]=R*sin(i*angle);
				Coil[0][i][2]=i*double(dist)/double(numpts)-dist/2.0;
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
	}

	std::cout<<"Helix with length = "<<dist<<" cm and along "<<axis<<" axis with "<<turns<<" * 180"<<std::endl;
}

void Biot_spiral(Parameters &pset, Vector<Vector<coord<> > > &Coil)
{
	int i;
	std::string section="";
	double R1=pset.getParamD("R1", section);
	double R2=pset.getParamD("R2", section);
	int turns=pset.getParamI("turns", section);
	int numpts=pset.getParamI("numpts", section);
	char axis=pset.getParamC("axis", section, false, 'z');
	coord<> center=pset.getParamCoordD("center", section, ',',false);
	double angle=2.0*PI*double(turns)/double(numpts-1);
	Coil.resize(1, Vector<coord<> >(numpts, 0.0));

	double rotangle=pset.getParamC("rotangle", section, false, 0)*DEG2RAD;
	coord<> rotaxis(0,0,1);
	if(rotangle!=0){
		rotaxis=pset.getParamCoordD("rotaxis", section);
		rotaxis/=norm(rotaxis);
	}

	switch(axis){
		case 'x':
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=0.0;
				Coil[0][i][1]=((R2-R1)*double(i)/double(numpts)+R1)*cos(double(i)*angle);
				Coil[0][i][2]=((R2-R1)*double(i)/double(numpts)+R1)*sin(double(i)*angle);
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
		case 'y':
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=((R2-R1)*double(i)/double(numpts)+R1)*cos(double(i)*angle);
				Coil[0][i][1]=0.0;
				Coil[0][i][2]=((R2-R1)*double(i)/double(numpts)+R1)*sin(double(i)*angle);
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
		case 'z':
		default:
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=((R2-R1)*double(i)/double(numpts)+R1)*cos(double(i)*angle);
				Coil[0][i][1]=((R2-R1)*double(i)/double(numpts)+R1)*sin(double(i)*angle);
				Coil[0][i][2]=0.0;
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
			}
			break;
	}

	std::cout<<"Spiral Shape with "<<turns<<" spires"<<std::endl;
}


void Biot_helmholtz(Parameters &pset, Vector<Vector<coord<> > > &Coil)
{
	int i;
	std::string section="";
	double R=pset.getParamD("R", section);
	double distance=pset.getParamD("length", section);
	coord<> center=pset.getParamCoordD("center", section, ',', false);
	int numpts=pset.getParamI("numpts", section);
	char axis=pset.getParamC("axis", section, false, 'z');
	Coil.resize(2,Vector<coord<> >(numpts, 0.0));
	double div=360.0/double(numpts-1), st=0.0;

	double rotangle=pset.getParamC("rotangle", section, false, 0)*DEG2RAD;
	coord<> rotaxis(0,0,1);
	if(rotangle!=0){
		rotaxis=pset.getParamCoordD("rotaxis", section);
		rotaxis/=norm(rotaxis);
	}

	switch(axis){
		case 'x':
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=-distance/2.0;
				Coil[0][i][1]=R*cos(st*PI/180.0);
				Coil[0][i][2]=R*sin(st*PI/180.0);
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
				st+=div;
			}
			st=0.0;
			for (i=0;i<numpts;i++){
				Coil[1][i][0]=distance/2.0;
				Coil[1][i][1]=R*cos(st*PI/180.0);
				Coil[1][i][2]=R*sin(st*PI/180.0);
				Coil[1][i]+=center;
				if(rotangle!=0)	Coil[1][i].rotate(rotangle, rotaxis);
				st+=div;
			}
			break;
		case 'y':
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=R*cos(st*PI/180.0);
				Coil[0][i][1]=-distance/2.0;
				Coil[0][i][2]=R*sin(st*PI/180.0);
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
				st+=div;
			}
			st=0.0;
			for (i=0;i<numpts;i++){
				Coil[1][i][0]=R*cos(st*PI/180.0);
				Coil[1][i][1]=distance/2.0;
				Coil[1][i][2]=R*sin(st*PI/180.0);
				Coil[1][i]+=center;
				if(rotangle!=0)	Coil[1][i].rotate(rotangle, rotaxis);
				st+=div;
			}
			break;
		case 'z':
		default:
			for (i=0;i<numpts;i++){
				Coil[0][i][0]=R*cos(st*PI/180.0);
				Coil[0][i][1]=R*sin(st*PI/180.0);
				Coil[0][i][2]=-distance/2.0;
				Coil[0][i]+=center;
				if(rotangle!=0)	Coil[0][i].rotate(rotangle, rotaxis);
				st+=div;
			}
			st=0.0;
			for (i=0;i<numpts;i++){
				Coil[1][i][0]=R*cos(st*PI/180.0);
				Coil[1][i][1]=R*sin(st*PI/180.0);
				Coil[1][i][2]=distance/2.0;
				Coil[1][i]+=center;
				if(rotangle!=0)	Coil[1][i].rotate(rotangle, rotaxis);
				st+=div;
			}

			break;
	}
	std::cout<<"Helmholtz axis along the "<<axis<<" axis"<<std::endl;
  	std::cout<<"  Diameter of the coils = "<<2.0*R<<" cm"<<std::endl;
  	std::cout<<"  Distance between coils = "<<distance<<" cm"<<std::endl;
//	if(Symmetry!=BiotSym::none)	Symmetry=BiotSym::C1;//no symmetry (the turning make it so...)
}


void Biot_helmholtz_true(Parameters &pset, Vector<Vector<coord<> > > &Coil)
{
	int i=0;
	std::string section="";
	double R=pset.getParamD("R", section);
	double distance=pset.getParamD("length", section);
	double layerw=pset.getParamD("layerwidth", section);
	double layerh=pset.getParamD("layerheight", section, false, layerw);
	int numlayers=pset.getParamI("numlayers", section);
	int turns=pset.getParamI("turns", section);
	coord<> center=pset.getParamCoordD("center", section, ',', false);
	int numpts=pset.getParamI("numpts", section);
	int ptsperlayer=int(numpts/numlayers);
	char axis=pset.getParamC("axis", section, false, 'z');
	Coil.resize(2,Vector<coord<> >(numpts, 0.0));

	double rotangle=pset.getParamC("rotangle", section, false, 0)*DEG2RAD;
	coord<> rotaxis(0,0,1);
	if(rotangle!=0){
		rotaxis=pset.getParamCoordD("rotaxis", section);
		rotaxis/=norm(rotaxis);
	}

	double div=PI2*double(turns)/double(ptsperlayer);
	double lstep=layerh*double(turns)/double(ptsperlayer);

	int ct=0;
	switch(axis){
		case 'x':
			for(int j=0;j<numlayers;++j){
				double trueR=R+double(j)*layerw/2.0;
				double sig=((j+1)%2==0)?1.0:-1.0; //we start winding from the bottom then go down
				for (i=0;i<ptsperlayer;i++){
					Coil[0][ct][0]=-distance/2.0-sig*(-(layerh*double(turns))/2.0+sig*layerh/2.0+lstep*double(i));
					Coil[0][ct][1]=trueR*cos(double(i)*div);
					Coil[0][ct][2]=trueR*sin(double(i)*div);
					Coil[0][ct]+=center;
					if(rotangle!=0)	Coil[0][ct].rotate(rotangle, rotaxis);

					++ct;
				}
			}
			ct=0;
			for(int j=0;j<numlayers;++j){
				double trueR=R+double(j)*layerw/2.0;
				double sig=((j+1)%2==0)?1.0:-1.0; //we start winding from the bottom then go down
				for (i=0;i<ptsperlayer;i++){
					Coil[1][ct][0]=distance/2.0-sig*(-(layerh*double(turns))/2.0+sig*layerh/2.0+lstep*double(i));
					Coil[1][ct][1]=trueR*cos(double(i)*div);
					Coil[1][ct][2]=trueR*sin(double(i)*div);
					Coil[1][ct]+=center;
					if(rotangle!=0)	Coil[1][ct].rotate(rotangle, rotaxis);
					++ct;
				}
			}
			break;
		case 'y':
			for(int j=0;j<numlayers;++j){
				double trueR=R+double(j)*layerw/2.0;
				double sig=((j+1)%2==0)?1.0:-1.0; //we start winding from the bottom then go down
				for (i=0;i<ptsperlayer;i++){
					Coil[0][ct][0]=trueR*cos(double(i)*div);
					Coil[0][ct][1]=-distance/2.0-sig*(-(layerh*double(turns))/2.0+sig*layerh/2.0+lstep*double(i));
					Coil[0][ct][2]=trueR*sin(double(i)*div);
					Coil[0][ct]+=center;
					if(rotangle!=0)	Coil[0][ct].rotate(rotangle, rotaxis);
					++ct;
				}
			}
			ct=0;
			for(int j=0;j<numlayers;++j){
				double trueR=R+double(j)*layerw/2.0;
				double sig=((j+1)%2==0)?1.0:-1.0; //we start winding from the bottom then go down
				for (i=0;i<ptsperlayer;i++){
					Coil[1][ct][0]=trueR*cos(double(i)*div);
					Coil[1][ct][1]=distance/2.0-sig*(-(layerh*double(turns))/2.0+sig*layerh/2.0+lstep*double(i));
					Coil[1][ct][2]=trueR*sin(double(i)*div);
					Coil[1][ct]+=center;
					if(rotangle!=0)	Coil[1][ct].rotate(rotangle, rotaxis);
					++ct;
				}
			}
			break;
		case 'z':
		default:
			for(int j=0;j<numlayers;++j){
				double trueR=R+double(j)*layerw/2.0;
				double sig=((j+1)%2==0)?1.0:-1.0; //we start winding from the bottom then go down
				for (i=0;i<ptsperlayer;i++){
					Coil[0][ct][0]=trueR*cos(double(i)*div);
					Coil[0][ct][1]=trueR*sin(double(i)*div);
					Coil[0][ct][2]=-distance/2.0-sig*(-(layerh*double(turns))/2.0+sig*layerh/2.0+lstep*double(i));
					Coil[0][ct]+=center;
					if(rotangle!=0)	Coil[0][ct].rotate(rotangle, rotaxis);
					++ct;
				}
			}
			ct=0;
			for(int j=0;j<numlayers;++j){
				double trueR=R+double(j)*layerw/2.0;
				double sig=((j+1)%2==0)?1.0:-1.0; //we start winding from the bottom then go down
				for (i=0;i<ptsperlayer;i++){
					Coil[1][ct][0]=trueR*cos(double(i)*div);
					Coil[1][ct][1]=trueR*sin(double(i)*div);
					Coil[1][ct][2]=distance/2.0-sig*(-(layerh*double(turns))/2.0+sig*layerh/2.0+lstep*double(i));
					Coil[1][ct]+=center;
					if(rotangle!=0)	Coil[1][ct].rotate(rotangle, rotaxis);
					++ct;
				}
			}
			break;
	}
	//coord<> diff=(Max-Min)/2.0;
	//for(int i=0;i<Coil.size();++i){
	//	for(int j=0;j<Coil[i].size();++j){
	//		Coil[i][j]+=center-diff;
	//	}
	//}


	std::cout<<"'True' Helmholtz axis along the "<<axis<<" axis"<<std::endl;
  	std::cout<<"  Diameter of the coils = "<<2.0*R<<" cm"<<std::endl;
  	std::cout<<"  Distance between coils = "<<distance<<" cm"<<std::endl;
//	if(Symmetry!=BiotSym::none)	Symmetry=BiotSym::C1;//no symmetry (the turning make it so...)
}



void Biot_saddle_coil(Parameters &pset, Vector<Vector<coord<> > > &Coil)
{
	int i;
	std::string section="";
	double R1=pset.getParamD("R1", section);
	double R2=pset.getParamD("R2", section);
	double height=pset.getParamD("height", section);
	double dist=pset.getParamD("dist", section);
	double dev1=pset.getParamD("dev1",section);
	double dev2=pset.getParamD("dev2", section);
	int NP=pset.getParamI("npts1", section);
	int NPa=pset.getParamI("npts2", section);
	int numpts=4*NP+4*NPa+1;
	Coil.resize(8);

	std::cout<<"Double Saddle Coil...with "<<numpts<<" Points"<<std::endl;
	int counter=0;
	Coil[0].resize(NPa);
	for (i=0;i<NPa;i++){  /*lower arc - counter clockwise */
		Coil[0][i][0]=R1*cos((((180.0-2.0*dev1)*double(i))/double(NPa))*PI/180.0+(dev1*PI)/180.0);
		Coil[0][i][1]=dist/2+R1*sin((((180.0-2.0*dev1)*double(i))/double(NPa))*PI/180.0+(dev1*PI)/180.0);
		Coil[0][i][2]=-height/2.0;
		counter++;
	}

	std::cout<<"Lower arc stops at "<<counter<<std::endl;

	Coil[1].resize(NP);
	for (i=0;i<NP;i++){
		Coil[1][i][0]=R1*cos((180.0-dev1)*PI/180.0)+(R2*cos((180.0-dev2)*PI/180.0)-R1*cos((180.0-dev1)*PI/180.0))*double(i)/double(NP);
		Coil[1][i][1]=(R2*sin((180.0-dev2)*PI/180.0)-R1*sin((180.0-dev1)*PI/180.0))*double(i)/double(NP)+dist/2.0+R1*sin((180.0-dev1)*PI/180.0);
		Coil[1][i][2]=(height*double(i))/double(NP)-height/2.0;
		counter++;
	}

	Coil[2].resize(NPa);
	for (i=0;i<NPa;i++){  /*upper arc - clockwise */
		Coil[2][i][0]=R2*cos((((2*dev2-180)*double(i))/double(NPa))*PI/180.0+(180-dev2)*PI/180.0);
		Coil[2][i][1]=dist/2.0+R2*sin((((2*dev2-180)*double(i))/double(NPa))*PI/180.0+(180-dev2)*PI/180.0);
		Coil[2][i][2]=height/2.0;
		counter++;
	}

	std::cout<<"Upper arc stops at "<<counter<<std::endl;

	Coil[3].resize(NP);
	for (i=0;i<NP;i++){
		Coil[3][i][0]=R2*cos(dev2*PI/180.0)-(R2*cos(dev2*PI/180.0)-R1*cos(dev1*PI/180.0))*double(i)/double(NP);
		Coil[3][i][1]=dist/2+R2*sin(dev2*PI/180)-(R2*sin(dev2*PI/180.0)-R1*sin(dev1*PI/180.0))*double(i)/double(NP);
		Coil[3][i][2]=height/2.0-(height*double(i))/double(NP);
		counter++;
	}

	std::cout<<"Discontinuous segment at "<<counter<<std::endl;

	Coil[4].resize(NPa);
	for (i=0;i<NPa;i++){ /*lower arc - clockwise */
		Coil[4][i][0]=R1*cos((((-180+2*dev1)*double(i))/double(NPa))*PI/180.0-(dev1*PI)/180.0);
		Coil[4][i][1]=-dist/2.0+R1*sin((((-180+2*dev1)*double(i))/double(NPa))*PI/180.0-(dev1*PI)/180.0);
		Coil[4][i][2]=-height/2.0;
		counter++;
	}

	std::cout<<"Second lower arc stops at "<<counter<<std::endl;

	Coil[5].resize(NP);
	for (i=0;i<NP;i++){
		Coil[5][i][0]=R1*cos((180.0+dev1)*PI/180.0)+(R2*cos((180.0+dev2)*PI/180.0)-R1*cos((180.0+dev1)*PI/180.0))*double(i)/double(NP);
		Coil[5][i][1]=-dist/2.0+R1*sin((180.0+dev1)*PI/180.0)+(R2*sin((180.0+dev2)*PI/180.0)-R1*sin((180.0+dev1)*PI/180.0))*double(i)/double(NP);
		Coil[5][i][2]=(height*double(i))/double(NP)-height/2.0;
		counter++;
	}

	Coil[6].resize(NPa);
	for (i=0;i<NPa;i++){ /* upper arc- counter clockwise */
		Coil[6][i][0]=R2*cos((180-2*dev2)*i/NPa*PI/180.0+(180+dev2)*PI/180.0);
		Coil[6][i][1]=-dist/2.0+R2*sin((180-2*dev2)*i/NPa*PI/180.0+(180+dev2)*PI/180.0);
		Coil[6][i][2]=height/2;
		counter++;
	}
	std::cout<<"Second upper arc stops at "<<counter<<std::endl;

	Coil[7].resize(NP+1);
	for (i=0;i<NP+1;i++){
		Coil[7][i][0]=R2*cos(-dev2*PI/180.0)-(R2*cos(-dev2*PI/180.0)-R1*cos(-dev1*PI/180.0))*double(i)/double(NP);
		Coil[7][i][1]=-dist/2.0+R2*sin(-dev2*PI/180.0)-(R2*sin(-dev2*PI/180.0)-R1*sin(-dev1*PI/180.0))*double(i)/double(NP);
		Coil[7][i][2]=height/2.0-(height*double(i))/double(NP);
		counter++;
	}
	std::cout<<"The Saddle Coil"<<std::endl;
	std::cout<<"  R1 = "<<R1<<" cm , R2 = "<<R2<<" cm"<<std::endl;
	std::cout<<"  Height = "<<height<<" cm"<<std::endl;
	std::cout<<"  Distance between two saddles = "<<dist<<" cm along Y"<<std::endl;
	std::cout<<"  Deviations from 180 = "<<dev1<<" and "<<dev2<<" degrees"<<std::endl;
}

END_BL_NAMESPACE

#endif

