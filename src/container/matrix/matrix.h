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
/***************************************
	this acts as a frount end to the '_matrix.h'
	class...although one can use the _matrix as a functional object,
	this calss will contain many of the 'other' BLAS functions
	that are used over and over...it will hopefully do
	some automatic matrix type conversion

****************************************/


#ifndef _matrix_h_
#define _matrix_h_ 1

#include "container/matrix/_matrix.h"

BEGIN_BL_NAMESPACE



typedef _matrix<complex, FullMatrix> matrix;
typedef _matrix<double, FullMatrix> rmatrix;

typedef _matrix<complex, IdentityMatrix> imatrix;
typedef _matrix<double, IdentityMatrix> rimatrix;

typedef _matrix<complex, DiagonalMatrix> dmatrix;
typedef _matrix<double, DiagonalMatrix> rdmatrix;

typedef _matrix<complex, TriDiagonalMatrix> trimatrix;
typedef _matrix<double, TriDiagonalMatrix> rtrimatrix;

typedef _matrix<complex, HermitianMatrix> hmatrix;
//typedef _matrix<complex, FullMatrix> hmatrix;

typedef _matrix<double, SymmetricMatrix> smatrix;


typedef _matrix<scomplex, FullMatrix> matrixs;
typedef _matrix<float, FullMatrix> rmatrixs;

typedef _matrix<scomplex, IdentityMatrix> imatrixs;
typedef _matrix<float, IdentityMatrix> rimatrixs;

typedef _matrix<scomplex, DiagonalMatrix> dmatrixs;
typedef _matrix<float, DiagonalMatrix> rdmatrixs;

typedef _matrix<scomplex, TriDiagonalMatrix> trimatrixs;
typedef _matrix<float, TriDiagonalMatrix> rtrimatrixs;

typedef _matrix<scomplex, HermitianMatrix> hmatrixs;
//typedef _matrix<complex, FullMatrix> hmatrix;

typedef _matrix<float, SymmetricMatrix> smatrixs;



/*

class _matrixType {
	private:
		GeneralMatrix gm_;
		FullMatrix fm_;
		IdentityMatrix im_;
		HermitianMatrix hm_;
		SymmetricMatrix sm_;
		DiagonalMatrix dm_;

		MatrixType thetype;

	public:

		typedef GeneralMatrix structure;

		_matrixType(): thetype(general){}
		_matrixType(MatrixType in): thetype(in){}
		_matrixType(int r, int c): thetype(full){
			typedef FullMatrix structure;
			fm_.resize(r,c);
		}
		_matrixType(int r, int c, MatrixType in): thetype(in){
			switch(thetype){
				case full:
					//typedef FullMatrix structure;
					fm_.resize(r,c);
					break;
				case hermitian:
					//typedef HermitianMatrix structure;
					hm_.resize(r,c);
					break;
				case symmetric:
					//typedef SymmetricMatrix structure;
					sm_.resize(r,c);
					break;
				case identity:
					//typedef IdentityMatrix structure;
					im_.resize(r,c);
					break;
				case diagonal:
					//typedef DiagonalMatrix structure;
					dm_.resize(r,c);
					break;
				default:
					//typedef GeneralMatrix structure;
					gm_.resize(r,c);
					thetype=general;
					break;
			}
		}


};

*/



END_BL_NAMESPACE




#endif
