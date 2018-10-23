

/* Dcoil.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 05.12.02
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
	Dcoil.h--> calculates magnetic fields along grid point...
	and uses the registration ability of the BiotFunctions to
	add a new shape..namely a 'D' shaped cylindrical coil

*/

#ifndef _D_Coil_Biot_h_
#define _D_Coil_Biot_h_ 1

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

//the D coil is basically a helix with a flattened surface
// from theta=theta1...theta2
void Biot_Dcoil(Parameters &pset, Vector<Vector<coord<> > > &Coil);

//the D circle (or 'D') is basically a cirlce with a flattened surface
// from theta=theta1...theta2
void Biot_Dcircle(Parameters &pset, Vector<Vector<coord<> > > &Coil);


#endif

