/* exsituSim.h ********/


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
	exsituSim.h -> takes the grids, Bfields for a B0 and B1
	coil, and performs a spin simulation on them

	the names of the input saved Bfield/Coil files for each coil...
	Or the Vectors of Bfields

	it will run in parellel

*/


#ifndef _exsituSim_h_
#define _exsituSim_h_ 1

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;


class exsituSim{

	public:

	//the number of points in the FID
		int npts;
	//the spectral sweep width
		double sw;

	//the spin system and FID calculator
		oneFID<SolidSys> sys;

	//The Vector of B0fields
		Vector<Vector<coord<> > > B0fields;

	//The Vector of B1fields
		Vector<Vector<coord<> > > B1fields;

	//the enumflag for readin purposes
		enum FTypes{
			B0=1,
			B1=2
		};

		exsituSim():
			npts(512), sw(20000.0)
		{}

		exsituSim(SolidSys &sys_i, int npts_i=512, double sw_i=20000.):
			npts(npts_i), sw(sw_i), sys(sys_i, npts, sw)
		{}

	//read in one file which is on grid 'which'
		void read(std::string fname, FTypes type, int which=0);

	//read a chunk of fields
		void read(Vector<std::string> fname, FTypes type);

	//calculates the FID
		Vector<complex> FID(matrix &ro, matrix &det);
};






#endif

