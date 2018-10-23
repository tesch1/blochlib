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

/* multispinop.cc **********************

This class acts as the generator for spin operators for a mutli spin system

like the 'singlespinop' file i would not use these directly, use 'spin_sys'
as it will maintain the lists of each operator type for each spin and make sure
they are only generated once
*/

#ifndef _multispinop_h_
#define _multispinop_h_ 1

#include "QMspins/singlespinop.h"
#include "utils/constants.h"
#include "utils/utils.h"
//#include "matrix.h"
#include "QMspins/spinsyspars.h"

BEGIN_BL_NAMESPACE



//makes a full space matrix from a single spin matrix (Ix,Iy,...),
//a 'spin_syspars' and an int to specify where in the tensor product
//series the operator is...it is templated to ensure the matrix type is conserved
//it simply expanded versus the Ie of each spin
template<class mxT>
mxT mps_expand(spin_sysPars &A, mxT &inop, int where);


//The rest of the function below follow from 'singlespinop' file
//and use the above expand to generate each spin op given the spin
//desired

//Identity Matrix (an easy one)

rimatrix mps_Ie(spin_sysPars &A, int spin);
rimatrix Ie(spin_sysPars &A, int spin);

//Ix matrix (symmetric/real) matrix
smatrix mps_Ix(spin_sysPars &A, int spin);
smatrix Ix(spin_sysPars &A, int spin);

//Iy matrix (hermitian/complex) matrix
hmatrix mps_Iy(spin_sysPars &A, int spin);
hmatrix Iy(spin_sysPars &A, int spin);

//Iz matrix (real digonal) matrix
rdmatrix mps_Iz(spin_sysPars &A, int spin);
rdmatrix Iz(spin_sysPars &A, int spin);

//Ip matrix (real full) matrix
rmatrix mps_Ip(spin_sysPars &A, int spin);
rmatrix Ip(spin_sysPars &A, int spin);

//Im matrix (real full) matrix
rmatrix mps_Imi(spin_sysPars &A, int spin);
rmatrix Imi(spin_sysPars &A, int spin);

END_BL_NAMESPACE



#endif

