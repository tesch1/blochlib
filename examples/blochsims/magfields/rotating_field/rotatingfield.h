


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-23-01
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
several classes that use the BiotCoil and Biot classes to calculate
magnetic fields and make them time dependant to be used in the
'Offset' classes...
*/

#ifndef _bl_rotating_field_h_
#define _bl_rotating_field_h_ 1

#include "blochlib.h"



//This is the 'Time Dependant' Magnetic Field Coil Class Calculator
// for the assumption that we can change the currents in the
// field coils as fast as we wish at some oscilating frequency
//
// It attempts to model 3 perpendicular coils (liek a tri-helmholtz system)
//
// Need to inherit 'DyanamicField' for purposes of offset and dipole-dipole
// properly handleing the changed Bfield handleing in time

template<class TheGrid>
class WrBfield :
	public DynamicField
{
	public:
		TheGrid *grid;

		Vector<BiotCoil> Coils;

		Vector<coord<> > BxCoil; //field from X coil
		Vector<coord<> > ByCoil; //field from y coil
		Vector<coord<> > BzCoil; //field from Z coil

		//frequencies of rotation for each direction
		coord<> wr;
		coord<> angles; //

		WrBfield(){}
		WrBfield(TheGrid &ingrid, Parameters &pset, std::string sec="coil"):
			DynamicField()
		{
			grid=&ingrid;
			read(pset, sec, false);
		}

		~WrBfield(){	grid=NULL; 	}

		void read(Parameters &pset,std::string sec, bool readgrid)
		{
			pset.addSection(sec);
			Parameters subp(pset.section(sec));
			int numsec=pset.getParamI("numsec", sec);
			Coils.resize(numsec);
			std::string subbase=pset.getParamS("secbase", sec, false, "subcoil");
			for(int i=0;i<numsec;++i){
				Coils[i]=BiotCoil(subp, subbase+itost(i+1));
			}
		}

		void read(TheGrid &ingrid, Parameters &pset, std::string sec="coil")
		{
			grid=&ingrid;
			read(pset, sec, false);
		}


		int CoilSize() const
		{	int ll=0; for(int i=0;i<Coils.size();++i){	ll+=Coils[i].size();	}	return ll;	}

		void calculateField()
		{
			BxCoil.resize(grid->size(), 0.0);
			ByCoil.resize(grid->size(), 0.0);
			BzCoil.resize(grid->size(), 0.0);
			if(Coils.size()>=1){
				Biot<TheGrid> tmB(Coils[0], *grid);
				tmB.calculateField();
				BxCoil=tmB.Bfield;
			}
			if(Coils.size()>=2){
				Biot<TheGrid> tmB(Coils[1], *grid);
				tmB.calculateField();
				ByCoil=tmB.Bfield;
			}
			if(Coils.size()>=3){
				Biot<TheGrid> tmB(Coils[2], *grid);
				tmB.calculateField();
				BzCoil=tmB.Bfield;
			}
			if(Coils.size()>=4){
				std::cerr<<"Warning: WrBfield.calculateBfield()"<<std::endl;
				std::cerr<<" only 3 coils are allowed here, a Bx, By, Bz coils"<<std::endl;
				std::cerr<<" any other coils will not be concidered"<<std::endl;
			}
		}

		coord<> AverageField(coord<int> rotframe=OneType<coord<> >::one())
		{
			double Bz=0.0, By=Bz, Bx=Bz;

			for(int i=0;i<BxCoil.size();++i){
				if(rotframe.z()){
					Bz+=Bfield(i).z();
				}else{
					Bz=0.0;
				}
				if(rotframe.x()){
					Bx+=Bfield(i).x();
				}else{
					Bx=0.0;
				}
				if(rotframe.y()){
					By+=Bfield(i).y();
				}else{
					By=0.0;
				}
			}
			return coord<>(Bx/BxCoil.size(), By/BxCoil.size(), Bz/BxCoil.size());
		}

		coord<> MaxField(coord<int> rotframe){
			double sigz=1.0, sigx=1.0, sigy=1.0,
			Bz=-1e30, By=Bz, Bx=Bz,
			oldBz=Bz, oldBy=By, oldBx=Bx;

			for(int i=0;i<BxCoil.size();++i){
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


		Vector<coord<> > Bfield(){	return BxCoil+ByCoil+BzCoil;	}

		coord<> Bfield(double t, int indx)
		{
			coord<> tmp=BxCoil(indx)+ByCoil(indx)+BzCoil(indx);
			//cout<<"Pre: ["<<tmp<<"]"<<endl;
			tmp.Rotate3D(wr.y()*PI2*t, wr.x()*PI2*t, wr.z()*PI2*t);
			//cout<<"POST: ["<<tmp<<"]"<<endl;
			return tmp;
		}

		coord<> Bfield(int indx)
		{
			return BxCoil(indx)+ByCoil(indx)+BzCoil(indx);
		}


		void writeMatlab(std::string &matout)
		{
			matstream out(matout, ios::out | ios::binary);
			writeMatlab(out);
		}

		void writeMatlab(const char *matout)
		{
			matstream out(matout, ios::out | ios::binary);
			writeMatlab(out);
		}

		void writeMatlab(matstream &matout)
		{
			matout.put("B", Bfield());
			Vector<coord<> > Points(CoilSize(), 0.0);
			Vector<double> pdiv;
			int ct=0;
			for(int j=0;j<Coils.size();++j){
				for(int i=0;i<Coils[j].Coil.size();++i){
					pdiv.push_back(ct);
					for(int k=0;k<Coils[j].Coil[i].size();++k){
						Points[ct]=Coils[j].Coil[i][k];
						ct++;
					}
				}
			}
			pdiv.push_back(ct);
			matout.put("pdiv", pdiv);
			matout.put("P", Points);
			matout.put("grid", grid->data());
			matout.put("N", grid->dim());
			coord<> tmp=grid->Max()-grid->Min();
			matout.put("r", tmp);
			matout.put("res",grid->dr());
			matout.close();
		}


};


//This is the 'Time Dependant' Magnetic Field Coil Class Calculator
// this is a more realistic modle for the acctual aparatus...
// the field here can only be on or off (and + or -) there is no
// we cannot rotated the filed via a sin(wr t) due to the limitations
// of real hardware..
//
// we simply overwrite the 'Bfield(t, int)' function here. from WrField
// NOTE:  times are bounded like so t=[tbegin, tend)
//
// t=0---------t=t1---------t=t2------...
// |----dir0-----|---dir1----|----dir2...

template<class TheGrid>
class StepTimeBfield :
	public WrBfield<TheGrid>
{
	public:

		Vector<double> tsplit; //should be {0, 0.2, 0.4,0.6} meaningthe field s
		coord<> defaultdir; //starting magnetization direction
		 //for each time we have {(1,1,1), (-1,1,1) ...} meaning that
		 //the field switches directions at those times...
		Vector<coord<> > direction;

		StepTimeBfield(){}
		StepTimeBfield(TheGrid &ingrid, Parameters &pset, std::string sec="coil"):
			WrBfield<TheGrid>(ingrid, pset, sec)
		{
			defaultdir=coord<>(0.0,0.0,1.0);
		}

	//fills up the stepping vectors in 'MField' from a parameter file
		void readPulseParams(Parameters &pset, std::string sec="pulsefield")
		{
			pset.addSection(sec);
			readPulseParams(pset.section(sec));
		}

	//fills up the stepping vectors in 'MField' from a parameter file
		void readPulseParams(const Vector<std::string> &pulfie)
		{
			double maxtime=0.0;
			bool repeat=false;
			Vector<double> times;
			Vector<coord<> > dir;
			int on=0;
			for(int i=0;i<pulfie.size();++i)
			{
				Vector<std::string> parsT=parse_param(pulfie[i]);
				if(parsT.size()>=4 && parsT[0][0]!='#'){
					if(parsT[0]=="start"){
						addTimeDirection(0.0,coord<>(atof(parsT[1]),atof(parsT[2]),atof(parsT[3])));
					}else{
						times.push_back(atof(parsT[0]));
						dir.push_back(coord<>(atof(parsT[1]),atof(parsT[2]),atof(parsT[3])));
						addTimeDirection(times[on],  dir[on]);
						maxtime=max(maxtime, times[on]);
						on++;
					}
				}else if(parsT.size()>=1 && parsT[0][0]!='#'){
					if(parsT[0]=="repeat") repeat=true;
				}
			}
		}



		void addTimeDirection(double newt,const coord<> &dir)
		{

			//no elements yet, just push it back
			if(direction.size()==0 || tsplit.size()==0)
			{	tsplit.push_back(newt); direction.push_back(dir);	return;	}

			int tind=timeIndex(newt,1);

			//put element in the middle somewhere
			if(tind!= -1 && tind !=tsplit.size()-1){
				Vector<double> ttop(tsplit.size()+1);

				ttop(Range(Range::Start, tind))=tsplit(Range(Range::Start, tind));
				ttop(tind+1)=(newt);
				ttop(Range(tind+2, Range::End))=tsplit(Range(tind+1, Range::End));
				tsplit=ttop;


				Vector<coord<> > dirtop(direction.size()+1);
				dirtop(Range(Range::Start, tind))=direction(Range(Range::Start, tind));
				dirtop(tind+1)=(dir);
				dirtop(Range(tind+2, Range::End))=direction(Range(tind+1, Range::End));
				direction=dirtop;
				return;
			}

			//put element at begining
			if(tind==-1 && tsplit.size() !=0){
				Vector<double> ttop(tsplit.size()+1);
				ttop(0)=newt;
				ttop(Range(1, Range::End))=tsplit(Range::All);
				tsplit=ttop;


				Vector<coord<> > dirtop(direction.size()+1);
				dirtop(0)=(dir);
				dirtop(Range(1, Range::End))=direction(Range::All);
				direction=dirtop;
				return;
			}

			//elements are at the end...
			tsplit.push_back(newt);
			direction.push_back(dir);

		}

		//return the proper index for a given time
		//the adding flag tells me to skip the fmod if adding the time
		int timeIndex(double t, int adding=0)
		{
			if(tsplit.size()!=1) RunTimeAssert(tsplit.size()==direction.size());
			if(tsplit.size()==0)	return -1;
			double moo=((t>tsplit[tsplit.size()-1]) && !adding) ? fmod(t,tsplit[tsplit.size()-1]): t;
			//cout<<"fmod: "<<moo<<" T:"<<t<<" tspint: "<<tsplit[tsplit.size()-1]<<endl;
			if(moo<tsplit(0))	return -1;
			int i=0;
			for(i=0;i<tsplit.size()-1;++i)
			{
				if(moo>=tsplit(i) && moo<tsplit(i+1)){	return i;	}
			}
			if(moo>=tsplit(tsplit.size()-1)) return tsplit.size()-1;
			return i;
		}

		coord<> AveField(int tind)
		{
			if(tind==-1){
				return sum(defaultdir.x()*BxCoil+
					defaultdir.y()*ByCoil+
					defaultdir.z()*BzCoil)/BxCoil.size();
			}else{
				return sum(direction(tind).x()*BxCoil+
					direction(tind).y()*ByCoil+
					direction(tind).z()*BzCoil)/BxCoil.size();
			}
		}

		coord<> Bfield(double t, int indx)
		{
			int tind=timeIndex(t);
			if(tind!=-1){
				//cout<<"INDEX: "<<tind<<" TIME: "<<t<<" DIR: "<<direction(tind)<<endl;
				return (direction(tind).x()*BxCoil(indx)+
						direction(tind).y()*ByCoil(indx)+
						direction(tind).z()*BzCoil(indx));//-AveField(tind);
			}
			return (defaultdir.x()*BxCoil(indx)+
					defaultdir.y()*ByCoil(indx)+
					defaultdir.z()*BzCoil(indx));//-AveField(tind);
		}

		coord<> Bfield(int indx)
		{
			return (defaultdir.x()*BxCoil(indx)+
					defaultdir.y()*ByCoil(indx)+
					defaultdir.z()*BzCoil(indx));
		}

		Vector<coord<> > Bfield()
		{
			return defaultdir.x()*BxCoil+
					defaultdir.y()*ByCoil+
					defaultdir.z()*BzCoil;
		}

		double norm(double t){
			int tind=timeIndex(t);
			return std::norm(AveField(tind));
		}
};




#endif
