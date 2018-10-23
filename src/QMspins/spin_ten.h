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
/*
	spin_ten.h--> spin spherical tensors up to rank 2
*/


#ifndef _spin_ten_h_
#define _spin_ten_h_

#include "QMspins/space_ten.h"
#include "container/matrix/matrix.h"
#include "QMspins/spinsys.h"

BEGIN_BL_NAMESPACE


rimatrix T0(SpinSys &A, int on);
rmatrix T1(SpinSys &A , int on1, int m);

matrix T2(SpinSys &A, int on1, int on2, int m);

//quadrupole verions... (on1=on2)
matrix T2(SpinSys &A, int on1, int m);

matrix T2_rot(SpinSys &A, int on1, int on2, int m, double alpha, double beta, double gamma);
matrix T2_rot(SpinSys &A, int on1, int m, double alpha, double beta, double gamma);

rmatrix T_CSA(SpinSys &A, int on, int l, int m);
matrix T_Qq(SpinSys &A, int on1, int order,int m);
matrix T_Quad(SpinSys &A, int on1, int m,int order=1);
//gives quarupole tensors (second order) in the total space
// of 2 ro 4 th rank
rdmatrix T_QuadTotal(SpinSys &A, int on_, int rank);

hmatrix T_Dip(SpinSys &A, int on1, int on2, int m);
hmatrix T_J(SpinSys &A, int on1, int on2, int m=0, int strong=0);

rmatrix T11(SpinSys &A, int on);
rdmatrix T10(SpinSys &A, int on);
rmatrix T1m1(SpinSys &A, int on);

rmatrix T22(SpinSys &A, int on1, int on2);
rmatrix T21(SpinSys &A, int on1, int on2);
matrix T20(SpinSys &A, int on1, int on2);
rmatrix T2m1(SpinSys &A, int on1, int on2);
rmatrix T2m2(SpinSys &A, int on1, int on2);

rmatrix T22(SpinSys &A, int on);
rmatrix T21(SpinSys &A, int on);
matrix T20(SpinSys &A, int on);
rmatrix T2m1(SpinSys &A, int on);
rmatrix T2m2(SpinSys &A, int on);


rdmatrix T00_CSA(SpinSys &A, int on);

rdmatrix T10_CSA(SpinSys &A, int on1);
rdmatrix T11_CSA(SpinSys &A, int on1);
rdmatrix T1m1_CSA(SpinSys &A, int on1);

rdmatrix T20_CSA(SpinSys &A, int on);
rmatrix T21_CSA(SpinSys &A, int on);
rdmatrix T22_CSA(SpinSys &A, int on);
rmatrix T2m1_CSA(SpinSys &A, int on);
rdmatrix T2m2_CSA(SpinSys &A, int on);


rdmatrix T00_Dip(SpinSys &A, int on1, int on2, int homo);

rdmatrix T10_Dip(SpinSys &A, int on1, int on2, int homo);
rdmatrix T11_Dip(SpinSys &A, int on1, int on2, int homo);
rdmatrix T1m1_Dip(SpinSys &A, int on1, int on2, int homo);

matrix T20_Dip(SpinSys &A, int on, int on2, int homo);
rdmatrix T21_Dip(SpinSys &A, int on, int on2, int homo);
rdmatrix T22_Dip(SpinSys &A, int on, int on2, int homo);
rdmatrix T2m1_Dip(SpinSys &A, int on, int on2, int homo);
rdmatrix T2m2_Dip(SpinSys &A, int on, int on2, int homo);

bool checksys(SpinSys &A, int on);
void printerr(std::string, int);


END_BL_NAMESPACE


#endif

