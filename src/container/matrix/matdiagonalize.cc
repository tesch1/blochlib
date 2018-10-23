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


#ifndef _matdigaonalize_cc_
#define _matdigaonalize_cc_

#include "container/matrix/matdiagonalize.h"
#include "container/matrix/_matrix.h"
#include <stdlib.h>

BEGIN_BL_NAMESPACE


//special cases for the diag function given number type AND matrix structure
//full complex
/*
void diag(const _matrix<complex,FullMatrix> &inmat, _matrix<complex, DiagonalMatrix> &dia, _matrix<complex, FullMatrix> &U, const bool des){
	if(!inmat.issquare()){
		std::cerr<<"Error::matrix:diag"<<std::endl;
		std::cerr<<" matrix must be square to diagonalize..."<<std::endl<<std::endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}
	ComplexFullDiag(inmat, dia, U, des);
}


//complex hermitian
void diag(const _matrix<complex,HermitianMatrix> &inmat, _matrix<complex, DiagonalMatrix> &dia, _matrix<complex, FullMatrix> &U, const bool des){
	if(!inmat.issquare()){
		std::cerr<<"Error::matrix:diag"<<std::endl;
		std::cerr<<" matrix must be square to diagonalize..."<<std::endl<<std::endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}
	ComplexHermitianDiag(inmat, dia, U,des);
}

//complex Symmetric (Not possible....)
void diag(const _matrix<complex,SymmetricMatrix> &inmat, _matrix<complex, DiagonalMatrix> &dia, _matrix<complex, FullMatrix> &U, const bool des){
	if(!inmat.issquare()){
		std::cerr<<"Error::matrix:diag"<<std::endl;
		std::cerr<<" matrix must be square to diagonalize..."<<std::endl<<std::endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}
	std::cerr<<std::endl<<"Error:: diag()"<<std::endl;
	std::cerr<<" Symmetric matrix should not be of complex type...use Hermitian instead"<<std::endl;
	    BLEXCEPTION(__FILE__,__LINE__)
}
*/


END_BL_NAMESPACE



#endif


