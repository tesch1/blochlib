
/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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
	Complex<Ctype_T>.cc -->the few Complex<Ctype_T> methods that can be compiled..
*/

//whats left of the functions not yet defined for the Complex<Ctype_T> stuff

#include "container/complex.h"
#include <string>
#include <fstream>
#include <iostream>

BEGIN_BL_NAMESPACE


const std::string cmx_form="%6.2f";		// Output format std::string

const Complex<double> complex0(0,0);		// z = 0 : (0,0)
const Complex<double> complex1(1,0);		// z = 1 : (1,0)
const Complex<double> complexi(0,1);		// z = i : (0,1)

const Complex<float> scomplex0(0,0);		// z = 0 : (0,0)
const Complex<float> scomplex1(1,0);		// z = 1 : (1,0)
const Complex<float> scomplexi(0,1);		// z = i : (0,1)



END_BL_NAMESPACE
