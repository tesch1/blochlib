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

/**************************************************
 matdiagonalize.h

 1) contains all the diagonalization routines...like diag, Mexp and Mlog
 should be included from '_matrix.h'

 **************************************************/


#ifndef Mat_EXP_H_
#define Mat_EXP_H_ 1

//#include "math.h"
#include "container/Vector/Vector.h"
#include "container/rankType.h"
#include "container/matrix/_matrix.h"
#include "utils/utils.h"

BEGIN_BL_NAMESPACE


template<class expr>
_matrix<typename expr::numtype, FullMatrix> PadeSeries(const _matrixExpr<expr> &Aa)
{
	return PadeSeries(_matrix<typename expr::numtype, FullMatrix>(Aa));
}

template<class expr, class nt1>
_matrix<typename expr::numtype, FullMatrix> PadeSeries(const _matrixExpr<expr> &Aa, nt1 mul)
{
	return PadeSeries(_matrix<typename expr::numtype, FullMatrix>(Aa), mul);
}

template<class T, class INstructure>
_matrix<T, FullMatrix> PadeSeries(const _matrix<T, INstructure> &A)
{
	if(A.empty()) return A;

	if(!A.issquare())
	{
		BLEXCEPTION(" Matrix Must be Square....")
	}

	//sadly we must copy Aa....
	//_matrix<T, FullMatrix> A=Aa;
	double f;
	int e;
	double Anorm=Re(sum(abs(A.row(0))));
	for(int i=1;i<A.rows();++i) Anorm=max(Anorm, Re(sum(abs(A.row(i)))));

	f=std::frexp(Anorm, &e);
	int s = max(0,e+1);
	double pp=1.0/pow(2.0,double(s));
	// Scale A by power of 2 so that its norm is < 1/2 .
	//A /=pow(2.0,double(s));

	// Pade approximation for exp(A)
	_matrix<T, FullMatrix> X = A*pp;
	double c = 0.5;
	_matrix<T, FullMatrix> E = _matrix<double, IdentityMatrix>(A.rows(), A.cols()) + c*pp*A;
	_matrix<T, FullMatrix> D  = _matrix<double, IdentityMatrix>(A.rows(), A.cols())-  c*pp*A;
	_matrix<T, FullMatrix> cX(A.rows(), A.cols());
	int q = 6;
	bool p = 1;
	for(int k=2;k<=q;++k){
		c *= double(q-k+1) / double(k*(2.0*q-k+1));
		X = A*X*pp;
		cX = c*X;
		E+= cX;
		if(p){
			D+= cX;
		}else{
			D-= cX;
		}
		p = !p;
	}

	E = inv(D)*E;
		// Undo scaling by repeated squaring
	for (int k=1;k<=s;++k) E *= E;
	return E;
}

template<class T, class INstructure, class nt1>
_matrix<T, FullMatrix> PadeSeries(const _matrix<T, INstructure> &A, nt1 mul)
{
	if(A.empty()) return A;

	if(!A.issquare())
	{
		BLEXCEPTION(" Matrix Must be Square....")
	}

	//sadly we must copy Aa....
	//_matrix<T, FullMatrix> A=Aa;
	double f;
	int e;
	double Anorm=Re(sum(abs(A.row(0))));
	for(int i=1;i<A.rows();++i) Anorm=max(Anorm, Re(sum(abs(A.row(i)))));

	f=std::frexp(Anorm, &e);
	int s = max(0,e+1);
	nt1 pp=mul/pow(2.0,double(s));
	// Scale A by power of 2 so that its norm is < 1/2 .
	//A /=pow(2.0,double(s));

	// Pade approximation for exp(A)
	_matrix<T, FullMatrix> X = A*pp;
	double c = 0.5;
	_matrix<T, FullMatrix> E = _matrix<double, IdentityMatrix>(A.rows(), A.cols()) + c*pp*A;
	_matrix<T, FullMatrix> D  = _matrix<double, IdentityMatrix>(A.rows(), A.cols())-  c*pp*A;
	_matrix<T, FullMatrix> cX(A.rows(), A.cols());
	int q = 6;
	bool p = 1;
	for(int k=2;k<=q;++k){
		c *= double(q-k+1) / double(k*(2.0*q-k+1));
		X = A*X*pp;
		cX = c*X;
		E+= cX;
		if(p){
			D+= cX;
		}else{
			D-= cX;
		}
		p = !p;
	}

	E = inv(D)*E;
		// Undo scaling by repeated squaring
	for (int k=1;k<=s;++k) E *= E;
	return E;
}


template<class expr>
_matrix<typename expr::numtype, FullMatrix> Mexp(const _matrixExpr<expr> &Aa)
{
	return PadeSeries(_matrix<typename expr::numtype, FullMatrix>(Aa));
}

template<class expr, class nt2>
_matrix<typename expr::numtype, FullMatrix> Mexp( const _matrixExpr<expr> &in,  nt2 mul)
{
	//_matrix<typename expr::numtype,FullMatrix> tm(in);
	return PadeSeries(_matrix<typename expr::numtype,FullMatrix>(in)*mul);
}

//exp(matrix) ...General form
template<class T>
_matrix<T, FullMatrix> Mexp(const _matrix<T, FullMatrix> &in)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	diag(in,dia, U, true);
	dia=exp(mul*dia);
	return prop(U,dia);
	//return PadeSeries(in);
}

//exp(matrix) ...General form
template<class Ctype,  class nt2>
_matrix<Complex<Ctype>, FullMatrix> Mexp(const _matrix<Complex<Ctype>, FullMatrix> &in,  nt2 mul)
{
	_matrix<Complex<Ctype>, DiagonalMatrix> dia;
	_matrix<Complex<Ctype>, FullMatrix> U;

	diag(in,dia, U, true);
	dia=exp(mul*dia);
	return prop(U,dia);

	//return PadeSeries(in,mul);
}

//exp(matrix) ...General form
template<class T,  class nt2>
_matrix<T, FullMatrix> Mexp(const _matrix<T, FullMatrix> &in,  nt2 mul)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	diag(in,dia, U, true);
	dia=exp(mul*dia);
	return prop(U,dia);

	//return PadeSeries(in,mul);
}

//exp(matrix) ...Hermitian input form
template<class T>
_matrix<complex, FullMatrix> Mexp(const _matrix<T, HermitianMatrix> &in)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	//cout<<"POST"<<endl<<in*1e3<<endl;
	diag(in,dia, U, false);
	dia=chop(exp(dia));
	return prop(U, dia);
	//return PadeSeries(in);
}

//exp(matrix) ...Hermitian input form
template<class Ctype>
_matrix<Complex<Ctype>, FullMatrix> Mexp(const _matrix<Complex<Ctype>, HermitianMatrix> &in)
{
	_matrix<Complex<Ctype>, DiagonalMatrix> dia;
	_matrix<Complex<Ctype>, FullMatrix> U;

	//cout<<"POST"<<endl<<in*1e3<<endl;
	diag(in,dia, U, false);
	dia=chop(exp(dia));
	return prop(U, dia);
	//return PadeSeries(in);
}

//exp(matrix) ...Hermitian input form
template<class nt2, class Ctype>
_matrix<Complex<Ctype>, FullMatrix> Mexp(const _matrix<Complex<Ctype>, HermitianMatrix> &in, nt2 mul)
{
	_matrix<Complex<Ctype>, DiagonalMatrix> dia;
	_matrix<Complex<Ctype>, FullMatrix> U;
	//cout<<"POST"<<endl<<in<<endl;

	diag(in,dia, U, false);
	dia=chop(exp(mul*dia));
	return prop(U, dia);

	//return Mexp(_matrix<complex, FullMatrix>(in)*mul);
}
//exp(matrix) ...Symmetric input form
_matrix<double, FullMatrix> Mexp(const _matrix<double, SymmetricMatrix> &in);

//exp(matrix) ...Symmetric input form
_matrix<double, FullMatrix> Mexp(const _matrix<double, SymmetricMatrix> &in,  double mul);

//exp(matrix) ...copmplex Diagonal input form
template<class T>
_matrix<T, DiagonalMatrix> Mexp(const _matrix<T, DiagonalMatrix> &in)
{
	return exp(in);
}

//exp(matrix) ...copmplex Diagonal input form
template<class T, class nt2>
_matrix<OutType(T, nt2), DiagonalMatrix> Mexp(const _matrix<T, DiagonalMatrix> &in, nt2 mul)
{
	return exp(mul*in);
}
//exp(matrix) ...copmplex Diagonal input form
template<class T>
_matrix<T, DiagonalMatrix> Mexp(const _matrix<T, IdentityMatrix> &in)
{
	return exp(in);
}

//exp(matrix) ...copmplex Diagonal input form
template<class T, class nt2>
_matrix<OutType(T, nt2), DiagonalMatrix> Mexp(const _matrix<T, IdentityMatrix> &in, nt2 mul)
{
	return exp(mul*in);
}





//Mlog(matrix) ...General form
template<class T, class INstructure>
_matrix<complex, FullMatrix> Mlog(const _matrix<T, INstructure> &in)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(dia);
	return prop(U, _matrix<complex, INstructure>(dia));
}

//Mlog(matrix) ...General form
template<class Ctype, class INstructure>
_matrix<Complex<Ctype>, FullMatrix> Mlog(const _matrix<Complex<Ctype>, INstructure> &in)
{
	_matrix<Complex<Ctype>, DiagonalMatrix> dia;
	_matrix<Complex<Ctype>, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(dia);
	return prop(U, _matrix<Complex<Ctype>, INstructure>(dia));
}


//Mlog(matrix) ...General form
template<class T, class INstructure, class T2>
_matrix<complex, FullMatrix> Mlog(const _matrix<T, INstructure> &in, const T2 &mul)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(mul*dia);
	return prop(U, _matrix<complex, INstructure>(dia));
}

//Mlog(matrix) ...General form
template<class Ctype, class INstructure, class T2>
_matrix<Complex<Ctype>, FullMatrix> Mlog(const _matrix<Complex<Ctype>, INstructure> &in, const T2 &mul)
{
	_matrix<Complex<Ctype>, DiagonalMatrix> dia;
	_matrix<Complex<Ctype>, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(mul*dia);
	return prop(U, _matrix<complex, INstructure>(dia));
}


//Mlog(matrix) ...Hermitian input form
template<class T>
_matrix<complex, FullMatrix> Mlog(const _matrix<T, HermitianMatrix> &in)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(complexi*dia);
	return prop(U, _matrix<complex, FullMatrix>(dia));
}

//Mlog(matrix) ...Hermitian input form
template<class Ctype>
_matrix<Complex<Ctype>, FullMatrix> Mlog(const _matrix<Complex<Ctype>, HermitianMatrix> &in)
{
	_matrix<Complex<Ctype>, DiagonalMatrix> dia;
	_matrix<Complex<Ctype>, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(Complex<Ctype>(0,1)*dia);
	return prop(U, _matrix<Complex<Ctype>, FullMatrix>(dia));
}

//Mlog(matrix) ...Hermitian input form
template<class T, class nt2>
_matrix<complex, FullMatrix> Mlog(const _matrix<T, HermitianMatrix> &in, const nt2 &mul)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(mul*dia);
	return prop(U, _matrix<complex, FullMatrix>(dia));
}

//Mlog(matrix) ...Hermitian input form
template<class Ctype, class nt2>
_matrix<Complex<Ctype>, FullMatrix> Mlog(const _matrix<Complex<Ctype>, HermitianMatrix> &in, const nt2 &mul)
{
	_matrix<Complex<Ctype>, DiagonalMatrix> dia;
	_matrix<Complex<Ctype>, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(mul*dia);
	return prop(U, _matrix<Complex<Ctype>, FullMatrix>(dia));
}
//Mlog(matrix) ...Symmetric input form
_matrix<double, FullMatrix> Mlog(const _matrix<double, SymmetricMatrix> &in);

//exp(matrix) ...Symmetric input form
template<class nt2>
_matrix<double, FullMatrix> Mlog(const _matrix<double, SymmetricMatrix> &in, const nt2 &mul)
{
	_matrix<double, DiagonalMatrix> dia;
	_matrix<double, FullMatrix> U;

	diag(in,dia, U, false);
	dia=log(mul*dia);
	return prop(U, _matrix<double, FullMatrix>(dia));
}

//Mlog(matrix) ...copmplex Diagonal input form
template<class T>
_matrix<T, DiagonalMatrix> Mlog(const _matrix<T, DiagonalMatrix> &in)
{
	return log(in);
}

//Mlog(matrix) ...copmplex Diagonal input form
template<class T, class nt2>
_matrix<T, DiagonalMatrix> Mlog(const _matrix<T, DiagonalMatrix> &in, const nt2 &mul)
{
	return log(mul*in);
}

//Mlog(matrix) ...copmplex Diagonal input form
template<class T>
_matrix<T, DiagonalMatrix> Mlog(const _matrix<T, IdentityMatrix> &in)
{
	return log(in);
}

//Mlog(matrix) ...copmplex Diagonal input form
template<class T, class nt2>
_matrix<T, DiagonalMatrix> Mlog(const _matrix<T, IdentityMatrix> &in, const nt2 &mul)
{
	return log(mul*in);
}

template<class expr>
_matrix<typename expr::numtype, typename expr::structure> Mlog( const _matrixExpr<expr> &in)
{
	_matrix<typename expr::numtype,typename  expr::structure> tm(in);
	return Mlog(tm);
}

template<class expr, class nt2>
_matrix<typename expr::numtype, typename expr::structure> Mlog( const _matrixExpr<expr> &in, const nt2 &mul)
{
	_matrix<typename expr::numtype,typename  expr::structure> tm(in);
	return Mlog(tm, mul);
}

/*template<class T, class INstructure>
_matrix<complex, HermitianMatrix> _matrix<T, INstructure>::exp(double tt)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	diag(dia, U, true);
	return prop(U, _matrix<complex, HermitianMatrix>(exp(-OneType<complex>::one()*tt*dia)));
}

template<class T, class INstructure>
_matrix<complex, HermitianMatrix> _matrix<T, INstructure>::exp(complex tt)
{
	_matrix<complex, DiagonalMatrix> dia;
	_matrix<complex, FullMatrix> U;

	diag(dia, U, true);
	*this=prop(U, _matrix<complex, HermitianMatrix>(exp(tt*dia)));
}
*/

END_BL_NAMESPACE


#endif

