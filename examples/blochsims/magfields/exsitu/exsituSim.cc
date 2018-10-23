/* exsituSim.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06.8.02
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
	exsituSim.cc -> takes the grids, Bfields for a B0 and B1
	coil, and performs a spin simulation on them

	the names of the input saved Bfield/Coil files for each coil...
	Or the Vectors of Bfields

	it will run in parellel

*/

#include "exsituSim.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

//read in one file which is on grid 'which'
void exsituSim::read(std::string fname, FTypes type, int which){}

//read a chunk of fields
void exsituSim::read(Vector<std::string> fname, FTypes type){}

//calculates the FID
Vector<complex> exsituSim::FID(matrix &roeq, matrix &detect)
{


	if(B0fields.size()<1){
		std::cerr<<std::endl<<"Warning: exsituSim.FID()"<<std::endl;
		std::cerr<<" Empty B0Fields..cannot go on"<<std::endl;
		exit(1);
	}

/*	if(B1fields.size()<1){
		std::cerr<<std::endl<<"Warning: exsituSim.FID()"<<std::endl;
		std::cerr<<" Empty B1Fields..cannot go on"<<std::endl;
		exit(1);
	}
*/
	if(sys.size()<1){
		std::cerr<<std::endl<<"Warning: exsituSim.FID()"<<std::endl;
		std::cerr<<" Empty Spin system..."<<std::endl;
		sys[0]="1H";
	}
	sys.setSpinMats();
	Vector<double> BaseOffsets(sys.csa.size(), 0.0);
	for(int i=0;i<sys.csa.size();++i){
		BaseOffsets[i]=sys.csa[i].iso();
	}

//find the mean B0-Z field;
	int totSize=0;
	double B0meanZ=0.0;
	for(int i=0;i<B0fields.size();++i){
		totSize+=B0fields[i].size();
		B0meanZ+=sum(B0fields[i]).z();
	}
	B0meanZ/=double(totSize);
	if(MPIworld.master()) cout<<"B0 z-mean Field in Gauss: "<<B0meanZ<<endl;
	if(MPIworld.master()) cout<<"B0 z-mean Field in Hz: "<<B0meanZ*sys[0].gammaGauss()/PI2<<endl;

	Vector<complex> fid(npts, 0.0);
	//first make sure we are Parrellel
	if(MPIworld.parallel()){
	//Perform the "master/slave" loops
		if(MPIworld.master()){
		//First we see if we send an initial request to each
			int endI=B0fields.size(), ctJ=0, ctI=0,done=-1, go, endJ=B0fields[0].size();
			for(int rank=1;rank<MPIworld.size();++rank){
				MPIworld.put(ctI, rank);
				MPIworld.put(ctJ, rank);
				ctJ++;
				if(ctJ>=endJ){
					ctI++; ctJ=0;
					endJ=B0fields[ctI].size();
					if(ctI>=endI) break;
				}
			}

		//keep sending to any proc that is willing to get one
			while(1){
				int get=MPIworld.getAny(go); //get the calculated FID
				MPIworld.put(ctI,get); //send the proc a new point to do
				MPIworld.put(ctJ,get); //send the proc a new point to do
				ctJ++;
				if(ctJ>=endJ){
					ctI++; ctJ=0;
					if(ctI>=endI) break;
					endJ=B0fields[ctI].size();
				}
			}

		//get the last returns
			for(int qq=1;qq<MPIworld.size();++qq){ MPIworld.get(go,qq);	}

		//put the final kills
			for(int qq=1;qq<MPIworld.size();++qq)
				MPIworld.put(done, qq);


		}else{ //the slave procs...
			int curI, curJ;
			while(1)
			{
				MPIworld.get(curI,0);
				if(curI==-1) break;
				MPIworld.get(curJ,0);
				for(int csact=0;csact<BaseOffsets.size();csact++){
					double newIso=BaseOffsets[csact]+(B0fields[curI][curJ].z()-B0meanZ)*sys[sys.csa[csact].on()].gammaGauss()/PI2;
					sys.csa[csact].iso(newIso);
				}
				fid+=sys.FID(roeq, detect);
				MPIworld.put(curJ,0);
			}
		}

		MPIworld.reduce(fid, Reduce::Add);

	}else{ //Serial Mode
		for(int i=0;i<B0fields.size();++i){
			for(int j=0;j<B0fields[i].size();++j){
				for(int csact=0;csact<BaseOffsets.size();csact++){
					double newIso=BaseOffsets[csact]+(B0fields[i][j].z()-B0meanZ)*sys[sys.csa[csact].on()].gammaGauss()/PI2;
					sys.csa[csact].iso(newIso);
				}
				fid+=sys.FID(roeq, detect)/B0fields[i].size();
			}
		}
	}
	return fid;
}

