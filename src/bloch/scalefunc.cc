/* scalefunc.cc ********/


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
 	scalefunc.cc-->scaling functions for the 'Dipolar Field'
 	the Compilable pieces.....
 */

/* --the scale function 'scales' the caclulated Bd by some ammount as determined by the user
 * this is quite usefull to take into account 'other averaging effects' that are hard
 * to calcuate, or are too time consuming.
 *
 * for instance diffusion can roughly be acounted for by scaling the Bd by some
 * factor proportional to the distance away from the target point
 * the closer it is too the target point the LESS effect Bd has as diffusion will
 * certainly average the interaction out, but points far away will be almost
 * entirely inclded.
 *
 * few examples classes are found in 'scalefunc.h'::THEY MUST HAVE A FUNCTION CALLED 'FUNCTION'!!!
 * and take in the coord I the coord J, THe DISTANCE bewteen the
 * (or some norm you feel is good) and the current magnitization coord
 *
 * function(coord<> &r, coord<> &rp, doubele dist,  coord<> &Bd[i])
 *
 * Not all the argurments need to be used, this is just
 * to provide some uniformity to the resulting scaling functions
 *
 * */

#ifndef  _scalfunc_cc_
#define  _scalfunc_cc_ 1

#include "bloch/scalefunc.h"

BEGIN_BL_NAMESPACE



std::ostream &operator<<(std::ostream &oo, const NoScaleFunc &out)
{
	oo<<"ScaleFunction:: None"<<std::endl;
	return oo;
}


std::ostream &operator<<(std::ostream &oo, const HardSphere &out)
{
	oo<<"ScaleFunction:: HardSphere--cutoff="<<out.cutoff<<std::endl;
	return oo;
}




std::ostream &operator<<(std::ostream &oo, const InvHardSphere &out)
{
	oo<<"ScaleFunction:: Inverse HardSphere--cutoff="<<out.cutoff<<std::endl;
	return oo;
}





std::ostream &operator<<(std::ostream &oo, const HardShell &out)
{
	oo<<"ScaleFunction:: HardShell-- r1="<<out.r1<<" r2="<<out.r2<<std::endl;
	return oo;
}



std::ostream &operator<<(std::ostream &oo, const InvHardShell &out)
{
	oo<<"ScaleFunction:: Inverse HardShell(r1="<<out.r1<<" r2="<<out.r2<<"]"<<std::endl;
	return oo;
}






std::ostream &operator<<(std::ostream &oo, const Power &out)
{
	oo<<"ScaleFunction::Power("<<out.factor<<"*r^"<<out.power<<")"<<std::endl;
	return oo;
}



std::ostream &operator<<(std::ostream &oo, const SubtractPower &out)
{
	oo<<"ScaleFunction:: Subtract Power("<<out.factor<<"*1/r^"<<out.power<<")"<<std::endl;
	return oo;
}

std::ostream &operator<<(std::ostream &oo, const BoltzmannScale &out)
{
	oo<<"ScaleFunction:: Blotzmann("<<out.factor<<"^2*exp(r*"<<out.factor<<")*r^2)"<<std::endl;
	return oo;
}



std::ostream &operator<<(std::ostream &oo, const TanhScale &out)
{
	oo<<"ScaleFunction:: Tanh("<<out.factor<<"r)"<<std::endl;
	return oo;
}


END_BL_NAMESPACE


#endif



