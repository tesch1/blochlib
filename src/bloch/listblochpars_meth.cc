/* listblochpars_meth.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-8-01
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
 	listblochpars_meth.cc-->a large list of BlochParams<> <>  with Grids attached...compilable parts
 	ONLY NULL GRID PARTS
 */


#ifndef _listlochparas_meth_cc_
#define _listlochparas_meth_cc_ 1


/* maintains a list of bloch parameters
	there are two specializations,
		one with no grid (a 'NullGrid')
		and one with a GradientGrid (a grid snapped with a graident vector)
	all others inlcude simple grids.

*/
#include "bloch/listblochpars.h"


BEGIN_BL_NAMESPACE




//this finds the largest 'Mo' in Normal/real units of all the spins in the
//list...this is used for setting a Dimensionless Mo where one could have
//different concnetrations and gamma factors for each spin...

END_BL_NAMESPACE



#endif




