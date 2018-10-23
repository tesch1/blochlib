

/* biot_basic.h ********/


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
	biot_basic.h--> some basic shape functins and
	a function registration class to be used inside
	BiotCoil...

*/

#ifndef _Biot_Basic_h_
#define _Biot_Basic_h_ 1

#include "utils/params.h"
#include "utils/matlab5.h"
#include <map>

BEGIN_BL_NAMESPACE



/**** BASIC SHAPES **/
void Biot_circle(Parameters &pset, Vector<Vector<coord<> > > &Coil);
void Biot_line(Parameters &pset,Vector<Vector<coord<> > > &Coil);
void Biot_helix(Parameters &pset, Vector<Vector<coord<> > > &Coil);
void Biot_spiral(Parameters &pset, Vector<Vector<coord<> > > &Coil);
void Biot_helmholtz(Parameters &pset, Vector<Vector<coord<> > > &Coil);
void Biot_helmholtz_true(Parameters &pset, Vector<Vector<coord<> > > &Coil);
void Biot_saddle_coil(Parameters &pset, Vector<Vector<coord<> > > &Coil);

//this 'registers' the generation function with a string
typedef void (*genShape_t)(Parameters &pset, Vector<Vector<coord<> > > & Coil) ;
//typedef map<std::string, void *> genShapeCont_t;

class BiotFuncMap
{
	public:
		BiotFuncMap();

		typedef std::map<std::string, genShape_t> myMap_t;
		typedef genShape_t Func_T;
		myMap_t mydata;


		void insert(std::string name, genShape_t in);
		//{	mydata.insert(std::pair<std::string, genShape_t>(name, in));	}

		Func_T find(std::string type_);
		/*{
			typename myMap_t::iterator myit;
			myit=mydata.find(type_);
			if(myit!=mydata.end()){
				return myit->second;
			}
			return 0;
		}*/

		void print(std::ostream &oo);
};

extern BiotFuncMap BiotFunctions;
//
////Our Master Generation function container
//
//extern genShapeCont_t BiotFunctions;
//void registerBiotShape( std::string name, void *);

END_BL_NAMESPACE

#endif

