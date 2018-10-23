
/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-11-01
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
	gentrain.cc--> methods for the base class of all the Time train engines
*/

#ifndef _gentimetrain_cc_
#define _gentimetrain_cc_ 1

#include "container/Vector/Vector.h"
#include <iostream>
#include "timetrain/gentrain.h"

BEGIN_BL_NAMESPACE


/***************** General Time Engine Base Class *************/
const void GenTimeEngine::VecAssignErr()
{
	BLEXCEPTION(" elements in rhs vector and out of assending order")
}

const void GenTimeEngine::Additerr()
{
	BLEXCEPTION(" itterated too far...")
}

const void GenTimeEngine::Adderr()
{
	BLEXCEPTION(std::string(" time, 't', inputed is smaller then the previous")+
			"\n  we require an absolute time, not a step time")
}

const void GenTimeEngine::AddStepErr()
{
	BLEXCEPTION(" step size, MUST be 1 or bigger...")
}

const void GenTimeEngine::SetErr()
{
	BLEXCEPTION(std::string(" either 'where' is too small or too big") +
			"\n  OR t<T(where-1) OR t>T(where+1)")
}


const void GenTimeEngine::Subiterr()
{
	BLEXCEPTION(" already at begining of times...")
}

const void GenTimeEngine::NstepErr()
{
	BLEXCEPTION(" number of steps MUST be 1 or more")
}





GenTimeEngine &GenTimeEngine::operator =(const GenTimeEngine &rhs)
{
	if(this==&rhs) return *this;
	RecordSep=rhs.RecordSep;
	nele_=rhs.nele_;
	return *this;
}



//---------------------General Iterator---------------

GenTimeIter &GenTimeIter::operator =(const GenTimeIter &rhs)
{
	if(this==&rhs) return *this;
	mye_=rhs.mye_;
	curpos_=rhs.curpos_;
	nextpos_=rhs.nextpos_;
	notended_=rhs.notended_;
	loops_=rhs.loops_;
	curloop_=rhs.curloop_;
	return *this;
}

void GenTimeIter::reset()
{
	curpos_=0;
	nextpos_=1;
	notended_=true;
}

END_BL_NAMESPACE


#endif
