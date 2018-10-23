/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-25-01
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
/* singlespinop.cc ****************************************************-*-c++-*

Generates simple ONE spin operators (given the quantum number of the system)

this is just a few functions that serve as starting points for the spin system
i do not recomend using these functions in any program....
instead use the "spin_sys" class with one spin as it stores theses output matrices
so they do not have to calculated again..

I'd like to thank to the makers of Gamma

  S.A. Smith, T.O. Levante, B.H. Meier and R.R. Ernst
  Computer Simulations in Magnetic Resonance:
  An Object-Oriented Programming Approach
  J. Magn. Reson., 106a, 75-105, (1994S, Series A)

http://gamma.magnet.fsu.edu/

for this intial set type of set up...

*/


#ifndef _singlespinop_h
#define _singlespinop_h 1

#include "container/complex.h"
#include "utils/constants.h"
#include "utils/utils.h"
#include "container/matrix/matrix.h"
#include "QMspins/spinsyspars.h"
#include <map>
#include <list>

BEGIN_BL_NAMESPACE


//these functions generate matrices of size hsXhs...

//these static lists are maintatined so we only calculate the matrix ONCE
extern std::map<int, rimatrix> 	IeList;		//real identity matrix Ie lists
extern std::map<int, smatrix> 	IxList;		//real symmetric matrix Ix lists
extern std::map<int, hmatrix>	IyList;		//complex hermieitan matrix Iy lists
extern std::map<int, rdmatrix>	IzList;		//real diagonal matrix Iz list
extern std::map<int, rmatrix>		IpList;		//complex full matrix Iplus list
extern std::map<int, rmatrix>		ImiList;	//complex full matrix Iminus list




//IDENTITY MATRIX
//a real identity matrix
rimatrix Ie(int hs);
rimatrix sps_Ie(int hs);

//Ix MATRIX
//a symmetric (real) matrix
smatrix Ix(int hs);
smatrix sps_Ix(int hs);

//Ix MATRIX
//a hermitian (complex) matrix
hmatrix Iy(int hs);
hmatrix sps_Iy(int hs);

//Iz MATRIX
//a diagonal (real) matrix
rdmatrix Iz(int hs);
rdmatrix sps_Iz(int hs);

//Iplus MATRIX
//a full (real) matrix
rmatrix Ip(int hs);
rmatrix sps_Ip(int hs);

//Iminus MATRIX
//a full (complex) matrix
rmatrix Imi(int hs);
rmatrix sps_Imi(int hs);

//Generates all the entries for that 'hs' dimension
void GenSingleOps(int hs);

END_BL_NAMESPACE



#endif

