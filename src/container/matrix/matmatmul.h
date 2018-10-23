/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 03-20-02
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
 * matmul.h
 	Contains all the various Matrix*Matrix multiply functions..
 */


#ifndef __matrix_MUL_h_
#define __matrix_MUL_h_ 1


#include "blochconfig.h"

BEGIN_BL_NAMESPACE


#ifdef AUX_BLAS_LIB

/* SPECIFIC MATRIX MULTIPLY ***/
enum CBLAS_ORDER {
	CblasRowMajor=101, /* row-major arrays */
	CblasColMajor=102}; /* column-major arrays */
enum CBLAS_TRANSPOSE {
	CblasNoTrans=111, /* trans='N' */
	CblasTrans=112, /* trans='T' */
	CblasConjTrans=113}; /* trans='C' */
enum CBLAS_UPLO {
	CblasUpper=121, /* uplo ='U' */
	CblasLower=122}; /* uplo ='L' */
enum CBLAS_DIAG {
	CblasNonUnit=131, /* diag ='N' */
	CblasUnit=132}; /* diag ='U' */
enum CBLAS_SIDE {
	CblasLeft=141, /* side ='L' */
	CblasRight=142}; /* side ='R' */


#ifdef __cplusplus
extern "C"{
#endif

void cblas_zgemm(
	const CBLAS_ORDER Order,
	const CBLAS_TRANSPOSE TransA,
	const CBLAS_TRANSPOSE TransB,
	const int M, const int N, const int K, const void *alpha, const
	void *A, const int lda, const void *B, const int ldb, const
	void *beta, void *C, const int ldc);

void cblas_cgemm(
	const CBLAS_ORDER Order,
	const CBLAS_TRANSPOSE TransA,
	const CBLAS_TRANSPOSE TransB,
	const int M, const int N, const int K, const void *alpha, const
	void *A, const int lda, const void *B, const int ldb, const
	void *beta, void *C, const int ldc);

void cblas_sgemm(
	const CBLAS_ORDER Order,
	const CBLAS_TRANSPOSE TransA,
	const CBLAS_TRANSPOSE TransB,
	const int M, const int N, const int K, const void *alpha, const
	void *A, const int lda, const void *B, const int ldb, const
	void *beta, void *C, const int ldc);

void cblas_dgemm(
	const CBLAS_ORDER Order,
	const CBLAS_TRANSPOSE TransA,
	const CBLAS_TRANSPOSE TransB,
	const int M, const int N, const int K, const double alpha,
	const double *A, const int lda, const double *B, const int ldb,
	const double beta, double *C, const int ldc);


// Bo's wrappers for the 'generic' matrix multiply
#ifdef __cplusplus
}
#endif



template<class T1, class T2>
void
DGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, float * A, int *Astride,
	   float* B, int *Bstride,
  	   T2 *beta, float* C, int *Cstride)
{
	CBLAS_TRANSPOSE tA=(transA[0]=='t' || transA[0]=='T')?CblasTrans:
						((transA[0]=='c' || transA[0]=='C')?CblasConjTrans:CblasNoTrans);
	CBLAS_TRANSPOSE tB=(transB[0]=='t' || transB[0]=='T')?CblasTrans:
						((transB[0]=='c' || transB[0]=='C')?CblasConjTrans:CblasNoTrans);

	float a=*alpha, b=*beta;
	cblas_sgemm(CblasColMajor, tA, tB,
       *M, *N, *K,
	  &a,  A, *Astride,
	    B, *Bstride,
  	   &b, C, *Cstride);
}

template<class T1, class T2>
void
DGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, double * A, int *Astride,
	   double* B, int *Bstride,
  	   T2 *beta, double * C, int *Cstride)
{
	double a= *alpha, b= *beta;
	CBLAS_TRANSPOSE tA=(transA[0]=='t' || transA[0]=='T')?CblasTrans:
						((transA[0]=='c' || transA[0]=='C')?CblasConjTrans:CblasNoTrans);
	CBLAS_TRANSPOSE tB=(transB[0]=='t' || transB[0]=='T')?CblasTrans:
						((transB[0]=='c' || transB[0]=='C')?CblasConjTrans:CblasNoTrans);

//	static int ct=0;
//	static std::clock_t clock;
//	double t1=std::clock();
	cblas_dgemm(CblasColMajor, tA, tB,
       *M, *N, *K,
	 a,
	    A, *Astride,
	    B, *Bstride,
	 b,
  	    C, *Cstride);
//  	 cout<<"time: "<<std::clock()-t1<<" ct: "<<ct++<<endl;
}

template<class T1, class T2>
void
DGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, Complex<double>* A, int *Astride,
	   Complex<double>* B, int *Bstride,
  	   T2 *beta, Complex<double>* C, int *Cstride)
{

	Complex<double> a=*alpha, b=*beta;
	CBLAS_TRANSPOSE tA=(transA[0]=='t' || transA[0]=='T')?CblasTrans:
						((transA[0]=='c' || transA[0]=='C')?CblasConjTrans:CblasNoTrans);
	CBLAS_TRANSPOSE tB=(transB[0]=='t' || transB[0]=='T')?CblasTrans:
						((transB[0]=='c' || transB[0]=='C')?CblasConjTrans:CblasNoTrans);
	cblas_zgemm(CblasColMajor, tA, tB,
       *M, *N, *K,
	  &a,  A, *Astride,
	    B, *Bstride,
  	   &b, C, *Cstride);
}

template<class T1, class T2>
void
DGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, Complex<float>* A, int *Astride,
	   Complex<float>* B, int *Bstride,
  	   T2 *beta, Complex<float>* C, int *Cstride)
{

	Complex<float> a=*alpha, b=*beta;
	CBLAS_TRANSPOSE tA=(transA[0]=='t' || transA[0]=='T')?CblasTrans:
						((transA[0]=='c' || transA[0]=='C')?CblasConjTrans:CblasNoTrans);
	CBLAS_TRANSPOSE tB=(transB[0]=='t' || transB[0]=='T')?CblasTrans:
						((transB[0]=='c' || transB[0]=='C')?CblasConjTrans:CblasNoTrans);
	cblas_cgemm(CblasColMajor, tA, tB,
       *M, *N, *K,
	  &a,  A, *Astride,
	    B, *Bstride,
  	   &b, C, *Cstride);
}

template<class T1, class T2>
void
CGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, Complex<float>* A, int *Astride,
	   Complex<float>* B, int *Bstride,
  	   T2 *beta, Complex<float>* C, int *Cstride)
{

	Complex<float> a=*alpha, b=*beta;
	CBLAS_TRANSPOSE tA=(transA[0]=='t' || transA[0]=='T')?CblasTrans:
						((transA[0]=='c' || transA[0]=='C')?CblasConjTrans:CblasNoTrans);
	CBLAS_TRANSPOSE tB=(transB[0]=='t' || transB[0]=='T')?CblasTrans:
						((transB[0]=='c' || transB[0]=='C')?CblasConjTrans:CblasNoTrans);
	cblas_cgemm(CblasColMajor, tA, tB,
       *M, *N, *K,
	 &a,  A, *Astride,
	    B, *Bstride,
  	   &b, C, *Cstride);
}

template<class T1, class T2>
void
SGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, float* A, int *Astride,
	   float* B, int *Bstride,
  	   T2 *beta, float* C, int *Cstride)
{

	float a=*alpha, b=*beta;
	CBLAS_TRANSPOSE tA=(transA[0]=='t' || transA[0]=='T')?CblasTrans:
						((transA[0]=='c' || transA[0]=='C')?CblasConjTrans:CblasNoTrans);
	CBLAS_TRANSPOSE tB=(transB[0]=='t' || transB[0]=='T')?CblasTrans:
						((transB[0]=='c' || transB[0]=='C')?CblasConjTrans:CblasNoTrans);
	cblas_sgemm(CblasColMajor, tA, tB,
       *M, *N, *K,
	  &a,  A, *Astride,
	    B, *Bstride,
  	   &b, C, *Cstride);
}


template<class T1, class T2>
void
ZGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, Complex<double>* A, int *Astride,
	   Complex<double>* B, int *Bstride,
  	   T2 *beta, Complex<double>* C, int *Cstride)
{

	complex a=*alpha, b=*beta;
	CBLAS_TRANSPOSE tA=(transA[0]=='t' || transA[0]=='T')?CblasTrans:
						((transA[0]=='c' || transA[0]=='C')?CblasConjTrans:CblasNoTrans);
	CBLAS_TRANSPOSE tB=(transB[0]=='t' || transB[0]=='T')?CblasTrans:
						((transB[0]=='c' || transB[0]=='C')?CblasConjTrans:CblasNoTrans);
	cblas_zgemm(CblasColMajor, tA, tB,
       *M, *N, *K,
	  &a,  A, *Astride,
	    B, *Bstride,
  	   &b, C, *Cstride);
}


#else
#include "container/matrix/dgemm.h"

/*template<class T,class T1, class T2>
void
DGEMM( char *transA, char *transB, int *M, int *N, int * K,
	  const T1 *alpha,const T* A, int *Astride,
	  const T* B,const int *Bstride,
	 const T2 *beta, T* C, int *Cstride);*/
#endif


//****************************SPECIAL MULTIPLICATION CASE!!***************
//	handles the matrix * matrix types, becuase we have an extra sum in the
//	acctuall operation we have to create a temporary....
// 	_matrix * _matrix
// 	_matrix * _matrixExpr
//	_matrixExpr * _matrix
//	_matrixExpr * _matrixExpr

#define MUL_STRUCTURE(mat1, mat2) \
_matrix<typename _matrixExprBinOpMUL<mat1, mat2 >::numtype, \
        typename _matrixExprBinOpMUL<mat1,mat2 >::structure> \


#define MUL_STRUCTURE2(mat1frag1, mat1frag2, mat2) \
_matrix<typename _matrixExprBinOpMUL< mat1frag1, mat1frag2, mat2 >::numtype, \
        typename _matrixExprBinOpMUL< mat1frag1, mat1frag2, mat2 >::structure> \


#define MUL_STRUCTURE3(mat1frag1, mat1frag2, mat2frag1, mat2frag2) \
_matrix<typename _matrixExprBinOpMUL< mat1frag1, mat1frag2,  mat2frag1, mat2frag2 >::numtype, \
        typename _matrixExprBinOpMUL< mat1frag1, mat1frag2, mat2frag1, mat2frag2 >::structure> \

/************* Basic Multi Function.... ********/


/******** IDENTITY MATRIX * 'Any' Matrix ***********/
//IdentityMatrix * Matrix
template<class nt1, class nt2, class st2>
inline
MUL_STRUCTURE3(_matrix<nt1, IdentityMatrix>, _matrix<nt2, st2>)
 operator *(const _matrix<nt1, IdentityMatrix> &m1,const _matrix<nt2, st2> &m2 ){
	return m2;
}

//Matrix*IdentityMatrix
template<class nt1, class st1,class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, st1>, _matrix<nt2, IdentityMatrix>)
 operator *(const _matrix<nt1, st1> &m1,const _matrix<nt2, IdentityMatrix> &m2 ){
	return m1;
}

//IdentityMatrix * Diagonal
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, IdentityMatrix>, _matrix<nt2, DiagonalMatrix>)
 operator *(const _matrix<nt1, IdentityMatrix> &m1,const _matrix<nt2, DiagonalMatrix> &m2 ){
	return m2;
}

//Diagonal*IdentityMatrix
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, DiagonalMatrix>, _matrix<nt2, IdentityMatrix>)
 operator *(const _matrix<nt1, DiagonalMatrix> &m1,const _matrix<nt2, IdentityMatrix> &m2 ){
	return m1;
}

//IdentityMatrix * Symmetric
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, IdentityMatrix>, _matrix<nt2, SymmetricMatrix>)
 operator *(const _matrix<nt1, IdentityMatrix> &m1,const _matrix<nt2, SymmetricMatrix> &m2 ){
	return m2;
}

//Symmetric*IdentityMatrix
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, SymmetricMatrix>, _matrix<nt2, IdentityMatrix>)
 operator *(const _matrix<nt1, SymmetricMatrix> &m1,const _matrix<nt2, IdentityMatrix> &m2 ){
	return m1;
}

//IdentityMatrix * Hermitian
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, IdentityMatrix>, _matrix<nt2, HermitianMatrix>)
 operator *(const _matrix<nt1, IdentityMatrix> &m1,const _matrix<nt2, HermitianMatrix> &m2 ){
	return m2;
}

//HermitianMatrix*IdentityMatrix
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, HermitianMatrix>, _matrix<nt2, IdentityMatrix>)
 operator *(const _matrix<nt1, HermitianMatrix> &m1,const _matrix<nt2, IdentityMatrix> &m2 ){
	return m1;
}

//IdentityMatrix * Full
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, IdentityMatrix>, _matrix<nt2, FullMatrix>)
 operator *(const _matrix<nt1, IdentityMatrix> &m1,const _matrix<nt2, FullMatrix> &m2 ){
	return m2;
}

//Full*IdentityMatrix
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, FullMatrix>, _matrix<nt2, IdentityMatrix>)
 operator *(const _matrix<nt1, FullMatrix> &m1,const _matrix<nt2, IdentityMatrix> &m2 ){
	return m1;
}

//--------------Expresions ---------------

template<class nt1, class ex2>
inline
_matrixExpr<ex2>
operator *(const _matrix<nt1, IdentityMatrix> & m1,const _matrixExpr<ex2> &expr2 ){
	return expr2;
}

template<class ex1, class nt2>
inline
_matrixExpr<ex1>
operator *(const _matrixExpr<ex1> &expr1, const _matrix<nt2, IdentityMatrix> & m2 ){
	return expr1;
}


/******** Diagonal MATRIX * 'Any' Matrix ***********/

//diagonal*diagonal
template<class nt1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, DiagonalMatrix>, _matrix<nt2, DiagonalMatrix>)
 operator *(const _matrix<nt1, DiagonalMatrix> &m1,const _matrix<nt2, DiagonalMatrix> &m2 ){
	typedef typename MatOutType(DiagonalMatrix,DiagonalMatrix) structure;
	typedef SumType(OutType(nt1, nt2)) numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();
	_matrix<numtype, DiagonalMatrix> dd(m1.rows(), m2.cols());
	for(int i=0;i<m1.rows();++i){
		dd.put(i,i, m1(i,i)*m2(i,i));
	}
	return dd;
}

//diagonal*Any
template<class nt1, class nt2, class st2>
inline
MUL_STRUCTURE3(_matrix<nt1, DiagonalMatrix>, _matrix<nt2, st2>)
 operator *(const _matrix<nt1, DiagonalMatrix> &m1,const _matrix<nt2, st2> &m2 ){
	typedef typename MatOutType(DiagonalMatrix,st2) structure;
	typedef SumType(OutType(nt1, nt2)) numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();

	_matrix<numtype, structure> dd(m1.rows(), m2.cols());
	typename structure::iterator iter(dd.rows(), dd.cols());
	while(iter){
		// <i|out|j> = <i|D|i><i|M|j>
		dd.put(iter.row(),iter.col(), m1(iter.row(),iter.row())*m2(iter.row(),iter.col()));
		++iter;
	}
	return dd;
}

//diagonal*Any
template<class nt1, class st1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, st1>, _matrix<nt2, DiagonalMatrix>)
 operator *(const _matrix<nt1, st1> &m1,const _matrix<nt2, DiagonalMatrix> &m2 ){
	typedef typename MatOutType(DiagonalMatrix,st1) structure;
	typedef SumType(OutType(nt1, nt2)) numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();

	_matrix<numtype, structure> dd(m1.rows(), m2.cols());
	typename structure::iterator iter(dd.rows(), dd.cols());
	while(iter){
		// <i|out|j> = <i|M|j><j|D|j>
		dd.put(iter.row(),iter.col(), m1(iter.row(),iter.col())*m2(iter.col(),iter.col()));
		++iter;
	}
	return dd;
}




/************ Hermitian MATRIX ***********/

/*
No difference between these and the 'Any' Any' one

//HermitianMatrix*Any
template<class nt1,  class nt2, class st2>
inline
MUL_STRUCTURE3(_matrix<nt1, HermitianMatrix>, _matrix<nt2, st2>)
 operator *(const _matrix<nt1, HermitianMatrix> &m1,const _matrix<nt2, st2> &m2 ){
	typedef typename MatOutType(HermitianMatrix,st2) structure;
	typedef SumType(OutType(nt1, nt2)) numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();

	_matrix<numtype, structure> dd(m1.rows(), m2.cols());
	typename structure::iterator iter(dd.rows(), dd.cols());
	numtype tmp=ZeroType<numtype>::zero();
	while(iter){
		// <i|out|j> += <i|H|k><k|f|j>
		tmp=ZeroType<numtype>::zero();
		for(int k=0;k<dd.cols();++k)	 tmp+=m1(iter.row(),k) * m2(k,iter.col());
		dd.put(iter.row(), iter.col(), tmp);
		++iter;
	}
	return dd;
}

//Any*Hermitian
template<class nt1,  class st1, class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, st1>, _matrix<nt2, HermitianMatrix>)
 operator *(const _matrix<nt1, st1> &m1,const _matrix<nt2, HermitianMatrix> &m2 ){
	typedef typename MatOutType(HermitianMatrix,st1) structure;
	typedef SumType(OutType(nt1, nt2)) numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();

	_matrix<numtype, structure> dd(m1.rows(), m2.cols());
	typename structure::iterator iter(dd.rows(), dd.cols());
	numtype tmp=ZeroType<numtype>::zero();
	while(iter){
		// <i|out|j> += <i|H|k><k|f|j>
		tmp=ZeroType<numtype>::zero();
		for(int k=0;k<dd.cols();++k)	 tmp+=m1(iter.row(),k) * m2(k,iter.col());
		dd.put(iter.row(), iter.col(), tmp);
		++iter;
	}
	return dd;
}

*/
/**** FULL * FULL ***/
template<class nt1,  class nt2>
inline
MUL_STRUCTURE3(_matrix<nt1, FullMatrix>, _matrix<nt2, FullMatrix>)
 operator *(const _matrix<nt1, FullMatrix> &m1,const _matrix<nt2, FullMatrix> &m2 ){
	typedef typename MatOutType(FullMatrix,FullMatrix) structure;
	typedef SumType(OutType(nt1, nt2)) numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();
	_matrix<numtype, structure> dd(m1.rows(), m2.cols());

		numtype *ddD=const_cast<numtype *>(dd.data());
		numtype *aD=const_cast<numtype *>(m1.data());
		numtype *bD=const_cast<numtype *>(m2.data());
		int M=m1.rows(), N=m2.cols(), K=m1.cols();
		static numtype alp(1);
		static numtype beta(0);
		static char *cc="n";
		BL_NAMESPACE::DGEMM(cc,cc, &M, &N, &K,
				&alp, aD, &M,
				bD, &N,
			&beta, ddD, &M);

	return dd;
}

#ifdef AUX_BLAS_LIB

/**** FULL * FULL  COMPLEX( CAN BE COMPILED!) ***/

_matrix<Complex<double>, FullMatrix>
 operator *(const _matrix<Complex<double>, FullMatrix> &m1,
 			const _matrix<Complex<double>, FullMatrix> &m2 );

_matrix<Complex<float>, FullMatrix>
 operator *(const _matrix<Complex<float>, FullMatrix> &m1,
 			const _matrix<Complex<float>, FullMatrix> &m2 );

_matrix<double, FullMatrix>
 operator *(const _matrix<double, FullMatrix> &m1,
 			const _matrix<double, FullMatrix> &m2 );

_matrix<float, FullMatrix>
 operator *(const _matrix<float, FullMatrix> &m1,
 			const _matrix<float, FullMatrix> &m2 );

#endif

/**** ANY * ANY ***/
//Any*Any
template<class nt1,  class st1, class nt2, class st2>
inline
MUL_STRUCTURE3(_matrix<nt1, st1>, _matrix<nt2, st2>)
 operator *(const _matrix<nt1, st1> &m1,const _matrix<nt2, st2> &m2 ){
	typedef typename MatOutType(st2,st1) structure;
	typedef SumType(OutType(nt1, nt2)) numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();

	_matrix<numtype, structure> dd(m1.rows(), m2.cols());
	/*typename structure::iterator iter(dd.rows(), dd.cols());
	//numtype tmp=ZeroType<numtype>::zero();
	while(iter){
		// <i|out|j> += <i|M1|k><k|M2|j>
		dd(iter.row(), iter.col())=ZeroType<numtype>::zero();
		for(int k=0;k<dd.cols();++k)
			dd(iter.row(), iter.col())+=m1(iter.row(),k) * m2(k,iter.col());
		//dd.put(iter.row(), iter.col(), tmp);
		++iter;
	}*/
	static int i,j,k;
	for(i=0;i<dd.rows();++i){
		for(j=0;j<dd.cols();++j){
			dd(i,j)=m1(i,0) * m2(0,j);
			for(k=1; k<dd.cols();++k){
				dd(i,j)+=m1(i,k) * m2(k,j);
			}
		}
	}


	return dd;
}




//___________ Expresions __________________

template<class nt1, class st1, class ex2>
inline
MUL_STRUCTURE2(_matrix<nt1, st1>, _matrixExpr<ex2>)
 operator *(const _matrix<nt1, st1> & m1,const _matrixExpr<ex2> &expr2 )
{
	typedef _matrix<typename _matrixExpr<ex2>::numtype, typename _matrixExpr<ex2>::structure> Mat2;
	return m1*Mat2(expr2);
}

template<class ex1, class nt2, class st2>
inline
_matrix<typename _matrixExprBinOpMUL< _matrixExpr<ex1>, _matrix<nt2, st2> >::numtype,
        typename _matrixExprBinOpMUL< _matrixExpr<ex1>, _matrix<nt2, st2> >::structure >
operator *(const _matrixExpr<ex1> &expr1, const _matrix<nt2, st2> & m2 )
{
	typedef _matrix<typename _matrixExpr<ex1>::numtype, typename _matrixExpr<ex1>::structure> Mat1;
	return Mat1(expr1)*m2;
}



/******* EXPR * EXPR *****/
template<class ex1, class ex2>
inline
 MUL_STRUCTURE(_matrixExpr<ex1>, _matrixExpr<ex2>)
operator *(const _matrixExpr<ex1> &m1,const _matrixExpr<ex2> &m2)
{
	typedef _matrix<typename _matrixExpr<ex1>::numtype, typename _matrixExpr<ex1>::structure> Mat1;
	typedef _matrix<typename _matrixExpr<ex2>::numtype, typename _matrixExpr<ex2>::structure> Mat2;

	return Mat1(m1)*Mat2(m2);
}

END_BL_NAMESPACE


#endif
