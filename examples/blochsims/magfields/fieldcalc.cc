
/* fieldcalc.cc ********/


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
	fieldcalc.cc--> calculates magnetic fields along grid point...

*/

#include "blochlib.h"


//the required 2 namespaces
using namespace BlochLib;
using namespace std;

//The main runner for calulating magnetic filds over grids a set number of grid points

int main(int argc, char **argv)
{
	MPIworld.start(argc, argv);
	std::string parse="";
	int q=1;
	if(MPIworld.master())
		query_parameter(argc, argv, q++, "\n\tEnter Parmeter set file name:", parse);

	MPIworld.scatter(parse);

	Parameters pset(parse);
	std::string choose=pset.getParamS("section");
	typedef XYZshape<XYZfull> TheGrid;
	pset.addSection("grid");
	Grid<UniformGrid> g1(pset.getParamCoordD("min", "grid"),
	pset.getParamCoordD("max", "grid"),
	pset.getParamCoordI("dim", "grid"));
	TheGrid g2(g1, XYZfull());
	MultiBiot<TheGrid> mycoil(g2,pset, choose);
	mycoil.Controller=MPIworld;
	mycoil.calculateField();
	pset.addSection("params");
	if(MPIworld.master()){
		mycoil.writeMatlab(pset.getParamS("matout", "params", false, "field.mat"));
		mycoil.write(pset.getParamS("textout", "params", false, "shape.boit"));
	}

	MPIworld.end();

}

