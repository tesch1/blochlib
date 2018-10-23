/* A BLAS compatible interface for matrix-matrix generator
 *
 * This interface compatible with the BLAS DGEMM.
 *
 * Written by Dominic LAM <ctlam@po.eecs.berkeley.edu>
 *
 * $Id: dgemm.c,v 1.1 1995/10/10 05:51:03 ctlam Exp $
 */

#ifndef _BLPHPIP_GEMM_H_
#define _BLPHPIP_GEMM_H_ 1

//#include <stdio.h>
#include "container/matrix/mdmd_a1bc_md.h"
#include "container/matrix/mdmd_acbc_md.h"
#include "container/matrix/mdmdt_a1bc_md.h"
#include "container/matrix/mdmdt_acbc_md.h"
#include "container/matrix/mdtmd_a1bc_md.h"
#include "container/matrix/mdtmd_acbc_md.h"
#include "container/matrix/mdtmdt_a1bc_md.h"
#include "container/matrix/mdtmdt_acbc_md.h"
//#include "container/matrix/zgemm.h"

//complex mat *mat
//#include "container/matrix/mdmd_cm.h"


BEGIN_BL_NAMESPACE


#define maxFF(a,b)     (((a) > (b)) ? (a) : (b))




/* error handling routine */
template<class T>
void xerbla(char *dd, T *info)
{
	 printf("** On entry to %6s, parameter number %2i had an illegal value\n",
		dd, *info);
}

/*
template<class T, class T1>
void DGEMM(const char &transa,const  char &transb,
		const  int &M,const  int &N,const  int &K,
		const  T alpha,const complex *a,const  int lda,
		const complex *b,const  int &ldb,
		const  T1 &beta, complex *c, const int &ldc)
{
	ZGEMM(transa,transb,
		M,N,K,
		complex(alpha), a,lda,
		b,ldb,
		complex(beta), c, ldc);
}
*/
template<class T,class T1, class T2>
void DGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, T* A, int *Astride,
	   T* B, int *Bstride,
  	   T2 *beta, T* C, int *Cstride)
{
	dgemm_(transA, transB,
       M, N, K,
	   alpha,  A, Astride,
	    B, Bstride,
  	   beta, C, Cstride);
}

template<class T,class T1, class T2>
void ZGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, T* A, int *Astride,
	   T* B, int *Bstride,
  	   T2 *beta, T* C, int *Cstride)
{
	dgemm_(transA, transB,
       M, N, K,
	   alpha,  A, Astride,
	    B, Bstride,
  	   beta, C, Cstride);
}
/*
template<class T,class T1, class T2>
void SGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, T* A, int *Astride,
	   T* B, int *Bstride,
  	   T2 *beta, T* C, int *Cstride)
{
	dgemm_(transA, transB,
       M, N, K,
	   alpha,  A, Astride,
	    B, Bstride,
  	   beta, C, Cstride);
}

template<class T,class T1, class T2>
void CGEMM(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, T* A, int *Astride,
	   T* B, int *Bstride,
  	   T2 *beta, T* C, int *Cstride)
{
	dgemm_(transA, transB,
       M, N, K,
	   alpha,  A, Astride,
	    B, Bstride,
  	   beta, C, Cstride);
}
*/
template<class T,class T1, class T2>
void
zgemm_(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, T* A, int *Astride,
	   T* B, int *Bstride,
  	   T2 *beta, T* C, int *Cstride)
{
	dgemm_(transA, transB,
       M, N, K,
	   alpha,  A, Astride,
	    B, Bstride,
  	   beta, C, Cstride);
}

template<class T,class T1, class T2>
void
dgemm_(char *transA, char *transB,
       int *M, int *N, int *K,
	   T1 *alpha, T* A, int *Astride,
	   T* B, int *Bstride,
  	   T2 *beta, T* C, int *Cstride)
{

  int info = 0;

	T betaT=T(*beta);
	T alphaT=T(*alpha);
	if(!transA) transA="N";
	if(!transB) transB="N";

  /* error checking */

  if (*M < 0) {
    info = 3;
	xerbla ("DGEMM ", &info);
    return;
  }
  else if (*N < 0) {
    info = 4;
	xerbla ("DGEMM ", &info);
    return;
  }
  else if (*K < 0) {
    info = 5;
	xerbla ("DGEMM ", &info);
    return;
  }

  if (transA[0] == 'n' || transA[0] == 'N'){
    if (transB[0] == 'n' || transB[0] == 'N'){

      /* error checking */
      if (*Astride < maxFF(1,*M)){
	info = 8;
	xerbla ("DGEMM ", &info);
      } else if (*Bstride < maxFF(1,*K)) {
	info = 10;
	xerbla ("DGEMM ", &info);
      } else if (*Cstride < maxFF(1,*M)) {
	info = 13;
	xerbla ("DGEMM ", &info);
      }

      /* C = AB + beta*C */
      else if (*alpha == 1)
	mdmd_a1bc_md
	  (*N, *K, *M, B, A, C, *Bstride, *Astride, *Cstride, betaT);
      /* C = alpha*AB + beta*C */
      else
	mdmd_acbc_md
	  (*N, *K, *M, B, A, C, *Bstride, *Astride, *Cstride, alphaT, betaT);
    } else {

      /* error checking */
      if (transB[0] != 'c' && transB[0] != 'C' &&
	  transB[0] != 't' && transB[0] != 'T'){
	info = 2;
	xerbla ("DGEMM ", &info);
      } else if (*Astride < maxFF(1,*M)) {
	info = 8;
	xerbla ("DGEMM ", &info);
      } else if (*Bstride < maxFF(1,*N)) {
	info = 10;
	xerbla ("DGEMM ", &info);
      } else if (*Cstride < maxFF(1,*M)) {
	info = 13;
	xerbla ("DGEMM ", &info);
      }
      /* C = AB' + beta*C */
      else if (*alpha == 1)
	mdtmd_a1bc_md
	  (*N, *K, *M, B, A, C, *Bstride, *Astride, *Cstride, betaT);
      /* C = alpha*AB' + beta*C */
      else
	mdtmd_acbc_md
	  (*N, *K, *M, B, A, C, *Bstride, *Astride, *Cstride, alphaT, betaT);
    }
  } else {
    if (transB[0] == 'n' || transB[0]== 'N'){

      /* error checking */
      if (transA[0] != 'c' && transA[0] != 'C' &&
	  transA[0] != 't' && transA[0] != 'T') {
		info = 1;
		xerbla ("DGEMM ", &info);
      } else if (*Astride < maxFF(1,*K)) {
		info = 8;
		xerbla ("DGEMM ", &info);
      } else if (*Bstride < maxFF(1,*K)) {
		info = 10;
		xerbla ("DGEMM ", &info);
      } else if (*Cstride < maxFF(1,*M)) {
		info = 13;
		xerbla ("DGEMM ", &info);
      }

      /* C = A'B + beta*C */
      else if (*alpha == 1)
	mdmdt_a1bc_md
	  (*N, *K, *M, B, A, C, *Bstride, *Astride, *Cstride, betaT);

      /* C = alpha*A'B + beta*C */
      else
	mdmdt_acbc_md
	  (*N, *K, *M, B, A, C, *Bstride, *Astride, *Cstride, alphaT, betaT);

    } else {

      /* error checking */
      if (transA[0] != 'c' && transA[0] != 'C' &&
	  transA[0] != 't' && transA[0] != 'T') {
		info = 1;
		xerbla ("DGEMM ", &info);
      } else if (transB[0]!= 'c' && transB[0] != 'C' &&
	       transB[0] != 't' && transB[0] != 'T') {
			info = 2;
			xerbla ("DGEMM ", &info);
      } else if (*Astride < maxFF(1,*K)) {
			info = 8;
			xerbla ("DGEMM ", &info);
      } else if (*Bstride < maxFF(1,*N)) {
			info = 10;
			xerbla ("DGEMM ", &info);
      } else if (*Cstride < maxFF(1,*M)) {
			info = 13;
			xerbla ("DGEMM ", &info);
      }

      /* C = A'B' * beta*C */
      else if (*alpha == 1)
	mdtmdt_a1bc_md
	  (*N, *K, *M, B, A, C, *Bstride, *Astride, *Cstride, betaT);
     /* C = alpha*A'B' * beta*C */
      else
	mdtmdt_acbc_md
	  (*N, *K, *M, B, A, C, *Bstride, *Astride, *Cstride, alphaT, betaT);
   }
  }
}
/*

template<class T1, class T2>
void
DGEMM(const char &transA,const char &transB,const int &M,const int &N,const int& K,
	  const T1 &alpha,const complex* A,const int &Astride,
	  const complex* B,const int &Bstride,
	 const T2 &beta, complex* C,const int &Cstride)
{

  int info = 0;

  //error checking

  if (M < 0) {
    info = 3;
	xerbla ("DGEMM ", &info);
    return;
  }
  else if (N < 0) {
    info = 4;
	xerbla ("DGEMM ", &info);
    return;
  }
  else if (K < 0) {
    info = 5;
	xerbla ("DGEMM ", &info);
    return;
  }

  if (transA == 'n' || transA == 'N'){
    if (transB == 'n' || transB == 'N'){

      // error checking
      if (Astride < maxFF(1,M)){
	info = 8;
	xerbla ("DGEMM ", &info);
      } else if (Bstride < maxFF(1,K)) {
	info = 10;
	xerbla ("DGEMM ", &info);
      } else if (Cstride < maxFF(1,M)) {
	info = 13;
	xerbla ("DGEMM ", &info);
      }

      // C = AB + beta*C
      else if (alpha == 1)
	mdmd_a1bc_md_cm
	  (N, K, M, B, A, C, Bstride, Astride, Cstride, complex(beta));
      // C = alpha*AB + beta*C
      else
	mdmd_acbc_md_cm
	  (N, K, M, B, A, C, Bstride, Astride, Cstride, complex(alpha), complex(beta));
    } else {

      // error checking
      if (transB != 'c' && transB != 'C' &&
	  transB != 't' && transB != 'T'){
	info = 2;
	xerbla ("DGEMM ", &info);
      } else if (Astride < maxFF(1,M)) {
	info = 8;
	xerbla ("DGEMM ", &info);
      } else if (Bstride < maxFF(1,N)) {
	info = 10;
	xerbla ("DGEMM ", &info);
      } else if (Cstride < maxFF(1,M)) {
	info = 13;
	xerbla ("DGEMM ", &info);
      }
      // C = AB' + beta*C
      else if (alpha == 1)
	mdtmd_a1bc_md
	  (N, K, M, B, A, C, Bstride, Astride, Cstride, complex(beta));
      // C = alpha*AB' + beta*C
      else
	mdtmd_acbc_md
	  (N, K, M, B, A, C, Bstride, Astride, Cstride, complex(alpha), complex(beta));
    }
  } else {
    if (transB == 'n' || transB == 'N'){

      // error checking
      if (transA != 'c' && transA != 'C' &&
	  transA != 't' && transA != 'T') {
		info = 1;
		xerbla ("DGEMM ", &info);
      } else if (Astride < maxFF(1,K)) {
		info = 8;
		xerbla ("DGEMM ", &info);
      } else if (Bstride < maxFF(1,K)) {
		info = 10;
		xerbla ("DGEMM ", &info);
      } else if (Cstride < maxFF(1,M)) {
		info = 13;
		xerbla ("DGEMM ", &info);
      }

      // C = A'B + beta*C
      else if (alpha == 1)
	mdmdt_a1bc_md
	  (N, K, M, B, A, C, Bstride, Astride, Cstride, complex(beta));

      // C = alpha*A'B + beta*C
      else
	mdmdt_acbc_md
	  (N, K, M, B, A, C, Bstride, Astride, Cstride, complex(alpha), complex(beta));

    } else {

      // error checking
      if (transA != 'c' && transA != 'C' &&
	  transA != 't' && transA != 'T') {
	info = 1;
	xerbla ("DGEMM ", &info);
      } else if (transB != 'c' && transB != 'C' &&
	       transB != 't' && transB != 'T') {
	info = 2;
	xerbla ("DGEMM ", &info);
      } else if (Astride < maxFF(1,K)) {
	info = 8;
	xerbla ("DGEMM ", &info);
      } else if (Bstride < maxFF(1,N)) {
	info = 10;
	xerbla ("DGEMM ", &info);
      } else if (Cstride < maxFF(1,M)) {
	info = 13;
	xerbla ("DGEMM ", &info);
      }

      // C = A'B' * beta*C
      else if (alpha == 1)
	mdtmdt_a1bc_md
	  (N, K, M, B, A, C, Bstride, Astride, Cstride, complex(beta));
      // C = alpha*A'B' * beta*C
      else
	mdtmdt_acbc_md
	  (N, K, M, B, A, C, Bstride, Astride, Cstride, complex(alpha), complex(beta));
    }
  }
}

*/
END_BL_NAMESPACE



#endif
