
/* singlefield.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 7.11.02
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
	singlefield.cc--> calculates magnetic fields along grid point...

*/

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;


//The main runner for calulating magnetic filds over grids a set number of grid points

int main(int argc, char **argv)
{

//start MPI
	MPIworld.start(argc, argv);
	std::string parse="";
	int q=1;
//get the parameter file
	if(MPIworld.master())
		query_parameter(argc, argv, q++, "\n\tEnter Parmeter set file name:", parse);

//distribute the file name to all the nodes
	MPIworld.scatter(parse);

//decalare our Parameter set
	Parameters pset(parse);

//get the desired coil to calculate the field over
	std::string choose=pset.getParamS("section");

//add this section to the parameter set
	pset.addSection(choose);

//a typdef to make our grid typing easier
	typedef XYZshape<XYZfull> TheGrid;

//get the number of grids we wish to calculate the field over
	int numGrids=pset.getParamI("numGrids");

//the grids sectiosn will look like 'grid1', 'grid2'
	std::string baseGname="grid";

	for(int i=1;i<=numGrids;++i)
	{
		std::string Gname=baseGname+itost(i);
	//add the grid section
		pset.addSection(Gname);

	//set up our base rectangular grids
		Grid<UniformGrid> g1(pset.getParamCoordD("min", Gname),
		pset.getParamCoordD("max", Gname),
		pset.getParamCoordI("dim", Gname));

	//set up thte master shape
		TheGrid g2(g1, XYZfull());

	//set up our Biot field calculator class
		MultiBiot<TheGrid> mycoil(g2,pset, choose);
		mycoil.Controller=MPIworld;

	//calculate the field
		mycoil.calculateField();

	//dump out the info to files
		if(MPIworld.master()){
			mycoil.writeMatlab(pset.getParamS("matout", "params", false, "field"+itost(i)+".mat"));
			mycoil.write(pset.getParamS("textout", "params", false, "shape"+itost(i)+".boit"));
		}
	}

//end out MPI session
	MPIworld.end();

}

