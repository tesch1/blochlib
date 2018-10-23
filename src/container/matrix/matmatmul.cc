/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 03-25-02
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

/****************************************************************************
 * matmatmul.cc
 	Contains all the various Matrix*Matrix multiply functions..
 */


#ifndef __matrix_MUL_cc_
#define __matrix_MUL_cc_ 1

#include "container/matrix/_matrix.h"
//#include "container/matrix/dgemm.h"
//#include <stdlib.h>

BEGIN_BL_NAMESPACE

#ifdef AUX_BLAS_LIB

/**** FULL * FULL  COMPLEX ***/

_matrix<Complex<double>, FullMatrix>
 operator *(const _matrix<Complex<double>, FullMatrix> &m1,const _matrix<Complex<double>, FullMatrix> &m2 ){
	typedef MatOutType(FullMatrix,FullMatrix) structure;
	typedef Complex<double> numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();
	_matrix<numtype, structure> dd(m1.rows(), m2.cols());

	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());
	int M=m1.rows(), N=m2.cols(), K=m1.cols();
	static numtype alp(1.0,0.0), beta(0.0,0.0);
	static char *cc="n";
	ZGEMM(cc,cc, &M, &N, &K,
			&alp, aD, &M,
			bD, &N,
			&beta, ddD, &M);
	return dd;
}
_matrix<Complex<float>, FullMatrix>
 operator *(const _matrix<Complex<float>, FullMatrix> &m1,const _matrix<Complex<float>, FullMatrix> &m2 ){
	typedef MatOutType(FullMatrix,FullMatrix) structure;
	typedef Complex<float> numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();
	_matrix<numtype, structure> dd(m1.rows(), m2.cols());

		numtype *ddD=const_cast<numtype *>(dd.data());
		numtype *aD=const_cast<numtype *>(m1.data());
		numtype *bD=const_cast<numtype *>(m2.data());
		int M=m1.rows(), N=m2.cols(), K=m1.cols();
		static Complex<float> alp(1.0,0.0), beta(0.0,0.0);
		static char *cc="n";
		CGEMM(cc,cc, &M, &N, &K,
				&alp, aD, &M,
				bD, &N,
				&beta, ddD, &M);
	//}
	return dd;
}

_matrix<float, FullMatrix>
 operator *(const _matrix<float, FullMatrix> &m1,const _matrix<float, FullMatrix> &m2 ){
	typedef MatOutType(FullMatrix,FullMatrix) structure;
	typedef float numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();
	_matrix<numtype, structure> dd(m1.rows(), m2.cols());

	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());
	int M=m1.rows(), N=m2.cols(), K=m1.cols();
	static numtype alp(1.0), beta(0.0);
	static char *cc="n";
	SGEMM(cc,cc, &M, &N, &K,
			&alp, aD, &M,
			bD, &N,
			&beta, ddD, &M);

	return dd;
}

_matrix<double, FullMatrix>
 operator *(const _matrix<double, FullMatrix> &m1,const _matrix<double, FullMatrix> &m2 ){
	typedef MatOutType(FullMatrix,FullMatrix) structure;
	typedef double numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();
	_matrix<numtype, structure> dd(m1.rows(), m2.cols());

	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());
	int M=m1.rows(), N=m2.cols(), K=m1.cols();
	static numtype alp(1.0), beta(0.0);
	static char *cc="n";
	DGEMM(cc,cc, &M, &N, &K,
			&alp, aD, &M,
			bD, &N,
			&beta, ddD, &M);

	return dd;
}

//internla multiply for the same type..should be faster
//DOUBLE MUL
//This does This=rhs*this
//NOTE:: THE NMR PROPOGATOR CASE!!! NOT THE NORMAL ONE!!!
template<>
 _matrix<double, FullMatrix> &
 _matrix<double, FullMatrix>::
 operator*=(const _matrix<double, FullMatrix> &m2)
{
	typedef FullMatrix structure;
	typedef double numtype;
	if(m2.rows() != cols()) MulErr();
	if(empty()) return *this;

	numtype *bD; bD=new numtype[thetype.numelements()];

	std::memcpy(&bD[0], &data()[0], thetype.numelements()*sizeof(numtype));

	this->resize(m2.rows(), cols());
	numtype *ddD=data();
	numtype *aD=const_cast<numtype *>(m2.data());

	int M=m2.rows(), N=cols(), K=m2.cols();
	static numtype alp(1);
	static numtype beta(0);
	static char *cc="n";
	DGEMM(cc,cc, &M, &N, &K,
			&alp, aD, &M,
			bD, &N,
	&beta, ddD, &M);

	delete [] bD;
	return *this;
}

//internla multiply for the same type..should be faster
//COMPLEX MUL
//This does This=rhs*this
//NOTE:: THE NMR PROPOGATOR CASE!!! NOT THE NORMAL ONE!!!
template<>
 _matrix<Complex<double>, FullMatrix> &
 _matrix<Complex<double>, FullMatrix>::
 operator*=(const _matrix<Complex<double>, FullMatrix> &m2)
{
	typedef FullMatrix structure;
	typedef Complex<double> numtype;
	if(rows() != m2.cols()) MulErr();
	if(empty()) return *this;
	numtype *bD; bD=new numtype[thetype.numelements()];

	std::memcpy(&bD[0], &data()[0], thetype.numelements()*sizeof(numtype));

	this->resize(m2.rows(), cols());
	numtype *ddD=data();
	numtype *aD=const_cast<numtype *>(m2.data());

	int M=m2.rows(), N=cols(), K=m2.cols();
	static numtype alp(1);
	static numtype beta(0);
	static char *cc="n";
	ZGEMM(cc,cc, &M, &N, &K,
			&alp, aD, &M,
			bD, &N,
	&beta, ddD, &M);

	delete [] bD;
	return *this;

}

//internla multiply for the same type..should be faster
//COMPLEX MUL
//This does This=rhs*this
//NOTE:: THE NMR PROPOGATOR CASE!!! NOT THE NORMAL ONE!!!
template<>
 _matrix<Complex<float>, FullMatrix> &
 _matrix<Complex<float>, FullMatrix>::
 operator*=(const _matrix<Complex<float>, FullMatrix> &m2)
{
	typedef FullMatrix structure;
	typedef Complex<float> numtype;
	if(rows() != m2.cols()) MulErr();
	if(empty()) return *this;
	numtype *bD; bD=new numtype[thetype.numelements()];

	std::memcpy(&bD[0], &data()[0], thetype.numelements()*sizeof(numtype));

	this->resize(m2.rows(), cols());
	numtype *ddD=data();
	numtype *aD=const_cast<numtype *>(m2.data());

	int M=m2.rows(), N=cols(), K=m2.cols();
	static numtype alp(1);
	static numtype beta(0);
	static char *cc="n";
	CGEMM(cc,cc, &M, &N, &K,
			&alp, aD, &M,
			bD, &N,
	&beta, ddD, &M);

	delete [] bD;
	return *this;

}

//end using aux blas lib
#endif

END_BL_NAMESPACE

#endif
