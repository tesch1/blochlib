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

/* This file holds the 'daig' and 'Mexp' specialization that can be compiled...*/


#ifndef Mat_EXP_CC_
#define Mat_EXP_CC_

#include "container/matrix/matexp.h"
#include "container/matrix/_matrix.h"
#include <stdlib.h>



BEGIN_BL_NAMESPACE



//exp(matrix) ...Symmetric input form
_matrix<double, FullMatrix> Mexp(const _matrix<double, SymmetricMatrix> &in)
{
	_matrix<double, DiagonalMatrix> dia;
	_matrix<double, FullMatrix> U;

	diag(in,dia, U, true);
	dia=exp(dia);
	return prop(U, dia);

	//return Mexp(_matrix<double, FullMatrix>(in));
}

_matrix<double, FullMatrix> Mexp(const _matrix<double, SymmetricMatrix> &in,  double mul)
{
	_matrix<double, DiagonalMatrix> dia;
	_matrix<double, FullMatrix> U;

	diag(in,dia, U, true);
	dia=exp(mul*dia);
	return prop(U, dia);

	//return Mexp(_matrix<double, FullMatrix>(in)*mul);
}

//exp(matrix) ...Symmetric input form
_matrix<double, FullMatrix> Mlog(const _matrix<double, SymmetricMatrix> &in)
{
	/*
	_matrix<double, DiagonalMatrix> dia;
	_matrix<double, FullMatrix> U;

	diag(in,dia, U, true);
	dia=log(dia);
	return prop(U, _matrix<double, SymmetricMatrix>(dia));
	*/
	return Mlog(in, 1.0);
}


END_BL_NAMESPACE



#endif


