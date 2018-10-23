


/* biot.h ********/


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
	biot.h--> calculates magnetic fields along grid point...

	two classes here::
		pre) BiotSym --> several symmetry parameters to specify about the coil geometry
		1) the 'shape' generator--> generate the Coil shapes or reads them from a file
		2) the biot --> calculates the magnetic field along the grid points.....
*/

#ifndef _biot_h_
#define _biot_h_ 1

#include "container/containers.h"
#include "utils/params.h"
#include "utils/matlab5.h"
#include "bloch/biot_basic.h"
#include "mpi/mpi_controller.h"
#include <map>

BEGIN_BL_NAMESPACE


/********************BEGIN BIOTSYM ********/

//Biot Sym is a simple set of enums that allow the user to set various symetry parameters
// about the coil geometry...for instance a 'circle' has 'normal axial symmetry'..meaning
// if rotated around the direction normal to the circle's surface we have the magnetic field
//  this would be specified by  "BiotSym::Dinfh | Zdirection"
//  which means we only need to calculate a single slice above (or below) the coil
//  and rotated the result about Z and reflect it over


class BiotSym{
	public:
		enum Direction{
			X=0x00001,
			Y=0x00002,
			Z=0x00004
		};

		//not all the point groups, but some nice ones
		enum PointGroup{
			none=0x00008, //no sym
			C1=0x00010, //no sym
			C2=0x00020, //180 rotation/mirror plane (requires Axis of rotation)
			C3=0x00040, // 120 rotation (requires Axis of rotation)
			C4=0x00080, // 90 rotation (2 mirror planes) (requires Axis of rotation)
			Cinfh=0x00100, // aixal sym (infinite rotation) (requires Axis of rotation)
			Dinfh=0x00200, // aixal sym (infinite rotation) and middle mirror plane
						   // (requires Axis of rotation..mirror is assumed perpendicular to rotation)
			D2h=0x00400, //like the molecule C2H4 (requires Axis of rotation)
			D4h=0x00800 //like XeF4 (requires Axis of rotation) alot like C4
		};
};



/******************* END BIOTSYM ***********/


/********************BEGIN BIOTCOIL ********/
//
// this reads in the Parameter set "coil" from a file
// determins the shape, the number of points desired and
// caluclates the shape (or reads in a file)
//

/*  reads this type_ of filw...

coil{

	### GLOBAL OPTIONS

	#log file...
	#use "cout" for standard out, "cerr" for standard error
	# or specify a new file name...
	logfile cout

	#use 'type_ file' to indicate a read in of points from a file
	#use 'type_ circle' to indicate use of a built in shape
	type_ <circle, helix, file...>

	#needed only if "type_ file"
	filename <fname>

	#the number of wire loops
	loops 1000

	#the number of amps
	amps 2.5


	#These are options specific to a built in shape


	### CIRCLE

	#the number of points to chop up the coil into
	numpts 4000

	# the radius
	R 2

	### HELIX

	#the number of points to chop up the coil into
	numpts 4000

	#the raddius
	R 2

	#the Z height
	Z 4

	#number of turns in the helix
	turns 10

	### SPIRAL

	#the number of points to chop up the coil into
	numpts 4000

	#the top radius
	R1 2

	#the bottom radius
	R2 4

	#number of turns in the helix
	turns 10

	### OURCOIL --> A saddle coil of sorts..UNITS IN cm
	#the lower saddle point divisions
	npts1 800

	#the second saddle point divisions
	npts2 1600

	## The total number of points will be...4*npts1+4*npts2+1

	#height of the saddle
	height 6

	#inner radius
	R1 2

	#outer radius
	R2 4

	#distance b/w the two saddles
	dist 0

	#deviations from 180 degrees for the first saddle
	dev1 10

	#deviations from 180 degrees for the second saddle
	dev2 40

	### HELMHOLTZ
	# radius (cm)
	R 5

	#distance b/w the two coils (cm)
	dist 5
}

*/

class BiotCoil{

	//friend class Biot;
	private:
		static const int MAXPOINTS=32768;
		int numpts;
		std::string type_;
		std::string fname;

		//place to store the Bx, By,Bz constants fields if type_==constant
		coord<> constField_;

		Parameters pset;

		coord<> Min;
		coord<> Max;

		void FindMinMax();

		std::string section;

	public:

		double current; //the 'current' Here is acctually "Number of Loops*Amps"
		double amps;
		double nloops;


		BiotFuncMap::Func_T ShapeGenerate;

		int Symmetry; //the BiotSym::PointGroup
		int Direction; //the BiotSym::Direction

		Vector<Vector<coord<> > > Coil;

//		ofstream logfile; //a log file to dump things to...

		BiotCoil():
			numpts(0), type_(""), fname(""),constField_(0),
			 Min(0), Max(0), section("coil"),
			current(0.0), amps(0.0), nloops(0.0),
			ShapeGenerate(0)
		{}

		//reads a file
		BiotCoil(std::string infile, std::string sec="coil"):
			constField_(0),pset(infile),
			Min(0), Max(0), section(sec),
			current(0.0), amps(0.0), nloops(0.0),
			ShapeGenerate(0)
		{
			read(pset,sec);
		}

		//reads a Parameter chunk
		BiotCoil(Parameters &in,std::string sec="coil"):
			constField_(0),pset(in),
			Min(0), Max(0), section(sec),
			current(0.0), amps(0.0), nloops(0.0),
			ShapeGenerate(0)
		{
			read(pset,sec);
		}

		//reads a Parameter chunk
		BiotCoil(const Vector<std::string> &in,std::string sec="coil"):
			constField_(0),pset(in),
			Min(0), Max(0), section(sec),
			current(0.0), amps(0.0), nloops(0.0),
			ShapeGenerate(0)
		{
			read(pset,sec);
		}

		BiotCoil(std::string intype, int pts, double loops=1000, double inamps=2.5):
			numpts(pts),type_(intype),fname(""), constField_(0),
			Min(0), Max(0),
			current(loops*inamps), amps(inamps),nloops(loops),
			Symmetry(BiotSym::none),
			Direction(0)
		{
			calcShape();
		}

		BiotCoil(const Vector<coord<> > &pts, double loops=1000, double inamps=2.5):
			numpts(pts.size()),type_(""), fname(""),
			Min(0), Max(0),
			current(loops*inamps), amps(inamps),nloops(loops),
			Symmetry(BiotSym::none),
			Direction(0), Coil(1,pts)
		{}

		BiotCoil(const Vector<Vector<coord<> > > &pts, double loops=1000, double inamps=2.5):
			numpts(pts.size()),type_(""), fname(""),
			Min(0), Max(0),
			current(loops*inamps), amps(inamps),nloops(loops),
			Symmetry(BiotSym::none),
			Direction(0), Coil(pts)
		{}

		void operator=(const BiotCoil &rhs);

		void read(std::string in,std::string sec="coil");
		void read(Parameters &in,std::string sec="coil");
		void read(const Vector<std::string> &in,std::string sec="coil");
		void read();

		void calcShape();

		inline int size()const{
			int si=0;
			for(int i=0;i<Coil.size();++i)
			{	si+=Coil[i].size();	}
			return si;
		}

		void Translate(const coord<> &dir);
		void Translate(double x, double y, double z);

		coord<> min() const;
		coord<> max() const;

		void print(std::ostream &out);
		void print(){	print(std::cout);	}

	//returns type
		inline std::string type() const {	return type_;	}

	//sets the type
		inline std::string type(std::string &typ)
		{	return type_=typ;	}

	//returns the 'constant Field' value if the
	// type is 'constant'
		inline coord<> constField() const {	return constField_;	}

	//sets the constant Field
		inline coord<> constField(const coord<> &inF)
		{	return constField_=inF;	}

	//reads coil points from files...
		void readShape(std::string fname) ;
		void readShape(std::ifstream &fname) ;
	//writes coil points from files...
		void writeShape(std::string &out) ;
		void writeShape(std::fstream &out) ;
};

/********************END BIOTCOIL *************/


/********************BEGIN BIOTSYMCALC ***************/

//this class takes a Uniform Grid and a BiotCoil, ans produces
//a list of 'Points' that need to be calculated via the 'Biot.calculateField()'
// and a list of the 'symmetry' points that acquire the same Bfield once the
// calculated points have been calculated....


//this is a little piece that contains the vector of calculating points
// and the indexes to 'symmetry points'
class BiotSymBit{
	public:
		int index;
		Vector<int> symindex;
		Vector<coord<> > sign; //signs of the x,y,z componets in the symetry..
		BiotSymBit():
			index(0)
		{}
};

std::ostream &operator<<(std::ostream &out, const BiotSymBit &oo);

//the calculator
class BiotSymCalc {
	private:

	public:
		Vector<BiotSymBit> syms;
		BiotSymCalc(){}

		template<class Grid_t>
		BiotSymCalc(const Grid_t &ing,const BiotCoil &coil)
		{
			calcSym(ing, coil);
		}

		inline int size() const	{	return syms.size();	}

		BiotSymBit &operator[](int i);
		BiotSymBit &operator()(int i);

		BiotSymBit operator[](int i) const;
		BiotSymBit operator()(int i) const;

		template<class Grid_t>
		void calcSym(const Grid_t &grid, const BiotCoil &coil);

		template<class Grid_t>
		void calcDinfh(const Grid_t &grid, const BiotCoil &coil);
};


/********************BEGIN BIOT ***************/
//the calculator of the B field given a Grid and a Biot Coil
// of course this class is useless without the 'BiotCoil' class
// soo, it is rightly inherited.
//The Grid is set up in the parameter file as such

// *****NOTE::**** The grid dused here is simply a UNIFORM GRID!! ****
// no fancy shapes....the BiotSymCalc uses a single UniformGrid to make the
// symmetry piece possible...thus..only UniformGrid's here
/*

grid{
	# the bottom left corner
	min -1, -1, -1
	# the bottom right corner
	max 1,1,1

	#the number of steps on either side
	dims 10,10,10
}

*/


template<class Grid_t=XYZshape<XYZfull> >
class Biot : public BiotCoil
{
	public:
		Grid_t grid;
		BiotSymCalc symcalc;
		Vector<coord<> > Bfield;
		MPIcontroller Controller;
		Biot():
			BiotCoil(),
			grid(), Bfield(),Controller()
		{}

		Biot(std::string fname, std::string sec="coil"):
			BiotCoil(fname,sec),Controller()
		{
			read(BiotCoil::pset);
		}

		Biot(Parameters &in, std::string sec="coil"):
			BiotCoil(in, sec),Controller()
		{
			read(BiotCoil::pset);
		}

		Biot(Grid_t &ingrid, Parameters &in, std::string sec="coil"):
			BiotCoil(in, sec), grid(ingrid),Controller()
		{
			read(BiotCoil::pset, false);
		}

		Biot(const BiotCoil &in, const Grid_t &gg):
			BiotCoil(in),
			grid(gg),
			symcalc(grid, *this),
			Bfield(gg.size(), 0.0)
		{}
		void read(std::string &in,bool readgrid=true);
		void read(Parameters &in, bool readgrid=true);
		void read();

		bool calculateField();
		coord<> calculateField(int ct);

	//writes a matlab binary
		void writeMatlab(std::string fname);
		void writeMatlab(matstream &out);

	//writes a text file with both the
	//Shape and the Bfield
		void write(std::string fname);
		void write(std::fstream &out);

	//reads a text file with both the
	//Shape and the Bfield
		void read(std::ifstream &out);
};

/********************************************************************/
/*** Multi  BIOT CLASS ***/
/********************************************************************/

//A class the 'vectorizes' the input and allows use of multi
//biol coil types in creation of a larger coil
template<class Grid_t=XYZshape<XYZfull> >
class MultiBiot{
	private:
		Vector<BiotCoil> Coils;

	public:
		Grid_t grid;

		Vector<coord<> > Bfield;

		MPIcontroller Controller;

		MultiBiot(){}

		MultiBiot(Parameters &pset, std::string sec="coil");

		MultiBiot(Grid_t &ingrid, Parameters &pset, std::string sec="coil");

		void read(Parameters &pset,std::string sec="", bool readgrid=true);

		int CoilSize() const;

		void calculateField();

		coord<> MaxField(coord<int> rotframe=OneType<coord<int> >::one());
		coord<> AverageField(coord<int> rotframe=OneType<coord<int> >::one());

		coord<> maxField(coord<int> rotframe=OneType<coord<int> >::one())
		{	return MaxField(rotframe);	}

		coord<> averageField(coord<int> rotframe=OneType<coord<int> >::one())
		{	return AverageField(rotframe);	}

		void writeMatlab(std::string fname) ;
		void writeMatlab(matstream &matout) ;

	//reads coil points from files...
		void readShape(std::string fname) ;
		void readShape(std::ifstream &fname) ;
	//writes coil points from files...
		void writeShape(std::string out) ;
		void writeShape(std::fstream &out) ;

	//writes a text file with both the
	//Shape and the Bfield
		void write(std::string fname) ;
		void write(std::fstream &out) ;

	//reads a text file with both the
	//Shape and the Bfield
		void read(std::string fname) ;
		void read(std::ifstream &out) ;
};

END_BL_NAMESPACE

#include "bloch/biot_meth.h"

#endif
