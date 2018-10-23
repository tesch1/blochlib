
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifndef _MAT_MUL_BL_F2C_
#define _MAT_MUL_BL_F2C_ 1

#include<stdio.h>

int _blError(char *srname, int *info);
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))

static int c__1 = 1;

template<class T, class T1, class T2, class T3, class T4>
int dgemm_(char transa, char transb,
			int M, int N, int K,
			T alpha, T1 *a, int lda, T2 *b, int ldb,
			T3 beta,T4 *c__, int ldc)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
	    i__3, i__4, i__5, i__6, i__7, i__8;

    /* Local variables */
    static int ilen, jlen, info;
    static bool nota, notb;
    static int i__, j, l, ncola;

  	static int ispan, jspan, lspan, nrowa, nrowb;
    static T4 ch[4096]	/* was [64][64] */;
    static int ii, jj;
    static T4 t11, t12, t21, t22;
    static int ll;
    static int idepth, jdepth;
    static T4 ch1[64], ch2[64];

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */

/*     Purpose */
/*     ======= */

/*     DGEMM  performs one of the matrix-matrix operations */

/*     C := alpha*op( A )*op( B ) + beta*C, */

/*     where  op( X ) is one of */

/*     op( X ) = X   or   op( X ) = X', */

/*     alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/*     an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */

/*     Parameters */
/*     ========== */

/*     TRANSA - CHARACTER*1. */
/*     On entry, TRANSA specifies the form of op( A ) to be used in */
/*     the matrix multiplication as follows: */

/*     TRANSA = 'N' or 'n',  op( A ) = A. */

/*     TRANSA = 'T' or 't',  op( A ) = A'. */

/*     TRANSA = 'C' or 'c',  op( A ) = A'. */

/*     Unchanged on exit. */

/*     TRANSB - CHARACTER*1. */
/*     On entry, TRANSB specifies the form of op( B ) to be used in */
/*     the matrix multiplication as follows: */

/*     TRANSB = 'N' or 'n',  op( B ) = B. */

/*     TRANSB = 'T' or 't',  op( B ) = B'. */

/*     TRANSB = 'C' or 'c',  op( B ) = B'. */

/*     Unchanged on exit. */

/*     M      - INTEGER. */
/*     On entry,  M  specifies  the number  of rows  of the  matrix */
/*     op( A )  and of the  matrix  C.  M  must  be at least  zero. */
/*     Unchanged on exit. */

/*     N      - INTEGER. */
/*     On entry,  N  specifies the number  of columns of the matrix */
/*     op( B ) and the number of columns of the matrix C. N must be */
/*     at least zero. */
/*     Unchanged on exit. */

/*     K      - INTEGER. */
/*     On entry,  K  specifies  the number of columns of the matrix */
/*     op( A ) and the number of rows of the matrix op( B ). K must */
/*     be at least  zero. */
/*     Unchanged on exit. */

/*     ALPHA  - DOUBLE PRECISION. */
/*     On entry, ALPHA specifies the scalar alpha. */
/*     Unchanged on exit. */

/*     A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is */
/*     k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/*     Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
/*     part of the array  A  must contain the matrix  A,  otherwise */
/*     the leading  k by m  part of the array  A  must contain  the */
/*     matrix A. */
/*     Unchanged on exit. */

/*     LDA    - INTEGER. */
/*     On entry, LDA specifies the first dimension of A as declared */
/*     in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
/*     LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/*     least  max( 1, k ). */
/*     Unchanged on exit. */

/*     B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is */
/*     n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/*     Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
/*     part of the array  B  must contain the matrix  B,  otherwise */
/*     the leading  n by k  part of the array  B  must contain  the */
/*     matrix B. */
/*     Unchanged on exit. */

/*     LDB    - INTEGER. */
/*     On entry, LDB specifies the first dimension of B as declared */
/*     in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
/*     LDB must be at least  max( 1, k ), otherwise  LDB must be at */
/*     least  max( 1, n ). */
/*     Unchanged on exit. */

/*     BETA   - DOUBLE PRECISION. */
/*     On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/*     supplied as zero then C need not be set on input. */
/*     Unchanged on exit. */

/*     C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
/*     Before entry, the leading  m by n  part of the array  C must */
/*     contain the matrix  C,  except when  beta  is zero, in which */
/*     case C need not be set on entry. */
/*     On exit, the array  C  is overwritten by the  m by n  matrix */
/*     ( alpha*op( A )*op( B ) + beta*C ). */

/*     LDC    - INTEGER. */
/*     On entry, LDC specifies the first dimension of C as declared */
/*     in  the  calling  (sub)  program.   LDC  must  be  at  least */
/*     max( 1, m ). */
/*     Unchanged on exit. */


/*     Level 3 Blas routine. */

/*     -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */

/*     This code comes from a report entitled: */
/*     The IBM RISC System/6000 and Linear Algebra Operations, by */
/*     Jack J. Dongarra, Peter Mayes, and Giuseppe Radicati di Brozolo, */
/*     University of Tennessee Computer Science Tech Report: CS - 90 - 122. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. Executable Statements .. */

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
/*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows */
/*     and  columns of  A  and the  number of  rows  of  B  respectively. */

    /* Parameter adjustments */
    a_dim1 = lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;

    /* Function Body */
    nota = transa=='n';
    notb = transb=='n';
    if (nota) {
	nrowa = M;
	ncola = K;
    } else {
	nrowa = K;
	ncola = M;
    }
    if (notb) {
	nrowb = K;
    } else {
	nrowb = N;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && transa!='C' &&  transa!='T') {
	info = 1;
    } else if (! notb && transb!='C' &&  transb!='T') {
	info = 2;
    } else if (M < 0) {
	info = 3;
    } else if (N < 0) {
	info = 4;
    } else if (K < 0) {
	info = 5;
    } else if (lda < max(1,nrowa)) {
	info = 8;
    } else if (ldb < max(1,nrowb)) {
	info = 10;
    } else if (ldc < max(1,M)) {
	info = 13;
    }
    if (info != 0) {
	_blError("DGEMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (M == 0 || N == 0 || (alpha == 0. || K == 0) && beta == 1.) {
	return 0;
    }
    if (beta == 0.) {
	i__1 = N;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = M;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		c__[i__ + j * c_dim1] = 0.;
/* L20: */
	    }
/* L40: */
	}
    } else {
	i__1 = N;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = M;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		c__[i__ + j * c_dim1] = beta * c__[i__ + j * c_dim1];
/* L60: */
	    }
/* L80: */
	}
    }

/*     And if  alpha.eq.zero. */

    if (alpha == 0.) {
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := C + alpha*A*B. */

	    i__1 = K;
	    for (l = 1; l <= i__1; l += 64) {
/* Computing MIN */
		i__2 = 64, i__3 = K - l + 1;
		lspan = min(i__2,i__3);
		i__2 = M;
		for (i__ = 1; i__ <= i__2; i__ += 64) {
		    idepth = 2;
/* Computing MIN */
		    i__3 = 64, i__4 = M - i__ + 1;
		    ispan = min(i__3,i__4);
		    ilen = idepth * (ispan / idepth);
		    i__3 = i__ + ispan - 1;
		    for (ii = i__; ii <= i__3; ++ii) {
			i__4 = l + lspan - 1;
			for (ll = l; ll <= i__4; ++ll) {
			    ch[ll - l + 1 + (ii - i__ + 1 << 6) - 65] =
				    alpha * a[ii + ll * a_dim1];
/* L100: */
			}
/* L120: */
		    }
		    i__3 = N;
		    for (j = 1; j <= i__3; j += 64) {
			jdepth = 2;
/* Computing MIN */
			i__4 = 64, i__5 = N - j + 1;
			jspan = min(i__4,i__5);
			jlen = jdepth * (jspan / jdepth);
			i__4 = j + jlen - 1;
			i__5 = jdepth;
			for (jj = j; i__5 < 0 ? jj >= i__4 : jj <= i__4; jj +=
				 i__5) {
			    i__6 = i__ + ilen - 1;
			    i__7 = idepth;
			    for (ii = i__; i__7 < 0 ? ii >= i__6 : ii <= i__6;
				     ii += i__7) {
				t11 = 0.;
				t21 = 0.;
				t12 = 0.;
				t22 = 0.;
				i__8 = l + lspan - 1;
				for (ll = l; ll <= i__8; ++ll) {
				    t11 += ch[ll - l + 1 + (ii - i__ + 1 << 6)
					     - 65] * b[ll + jj * b_dim1];
				    t21 += ch[ll - l + 1 + (ii - i__ + 2 << 6)
					     - 65] * b[ll + jj * b_dim1];
				    t12 += ch[ll - l + 1 + (ii - i__ + 1 << 6)
					     - 65] * b[ll + (jj + 1) * b_dim1]
					    ;
				    t22 += ch[ll - l + 1 + (ii - i__ + 2 << 6)
					     - 65] * b[ll + (jj + 1) * b_dim1]
					    ;
/* L140: */
				}
				c__[ii + jj * c_dim1] += t11;
				c__[ii + 1 + jj * c_dim1] += t21;
				c__[ii + (jj + 1) * c_dim1] += t12;
				c__[ii + 1 + (jj + 1) * c_dim1] += t22;
/* L160: */
			    }
			    if (ilen < ispan) {
				i__7 = i__ + ispan - 1;
				for (ii = i__ + ilen; ii <= i__7; ++ii) {
				    t11 = 0.;
				    t12 = 0.;
				    i__6 = l + lspan - 1;
				    for (ll = l; ll <= i__6; ++ll) {
					t11 += ch[ll - l + 1 + (ii - i__ + 1
						<< 6) - 65] * b[ll + jj *
						b_dim1];
					t12 += ch[ll - l + 1 + (ii - i__ + 1
						<< 6) - 65] * b[ll + (jj + 1)
						* b_dim1];
/* L180: */
				    }
				    c__[ii + jj * c_dim1] += t11;
				    c__[ii + (jj + 1) * c_dim1] += t12;
/* L200: */
				}
			    }
/* L220: */
			}
			if (jlen < jspan) {
			    i__5 = j + jspan - 1;
			    for (jj = j + jlen; jj <= i__5; ++jj) {
				i__4 = i__ + ilen - 1;
				i__7 = idepth;
				for (ii = i__; i__7 < 0 ? ii >= i__4 : ii <=
					i__4; ii += i__7) {
				    t11 = 0.;
				    t21 = 0.;
				    i__6 = l + lspan - 1;
				    for (ll = l; ll <= i__6; ++ll) {
					t11 += ch[ll - l + 1 + (ii - i__ + 1
						<< 6) - 65] * b[ll + jj *
						b_dim1];
					t21 += ch[ll - l + 1 + (ii - i__ + 2
						<< 6) - 65] * b[ll + jj *
						b_dim1];
/* L240: */
				    }
				    c__[ii + jj * c_dim1] += t11;
				    c__[ii + 1 + jj * c_dim1] += t21;
/* L260: */
				}
				if (ilen < ispan) {
				    i__7 = i__ + ispan - 1;
				    for (ii = i__ + ilen; ii <= i__7; ++ii) {
					t11 = 0.;
					i__4 = l + lspan - 1;
					for (ll = l; ll <= i__4; ++ll) {
					    t11 += ch[ll - l + 1 + (ii - i__
						    + 1 << 6) - 65] * b[ll +
						    jj * b_dim1];
/* L280: */
					}
					c__[ii + jj * c_dim1] += t11;
/* L300: */
				    }
				}
/* L320: */
			    }
			}
/* L340: */
		    }
/* L360: */
		}
/* L380: */
	    }
	} else {

/*           Form  C := C + alpha*A'*B */

	    i__1 = M;
	    for (i__ = 1; i__ <= i__1; i__ += 64) {
		idepth = 2;
/* Computing MIN */
		i__2 = 64, i__3 = M - i__ + 1;
		ispan = min(i__2,i__3);
		ilen = idepth * (ispan / idepth);
		i__2 = K;
		for (l = 1; l <= i__2; l += 64) {
/* Computing MIN */
		    i__3 = 64, i__5 = K - l + 1;
		    lspan = min(i__3,i__5);
		    i__3 = i__ + ispan - 1;
		    for (ii = i__; ii <= i__3; ++ii) {
			i__5 = l + lspan - 1;
			for (ll = l; ll <= i__5; ++ll) {
			    ch[ll - l + 1 + (ii - i__ + 1 << 6) - 65] =
				    alpha * a[ll + ii * a_dim1];
/* L400: */
			}
/* L420: */
		    }
		    i__3 = N;
		    for (j = 1; j <= i__3; j += 64) {
			jdepth = 2;
/* Computing MIN */
			i__5 = 64, i__7 = N - j + 1;
			jspan = min(i__5,i__7);
			jlen = jdepth * (jspan / jdepth);
			i__5 = j + jlen - 1;
			i__7 = jdepth;
			for (jj = j; i__7 < 0 ? jj >= i__5 : jj <= i__5; jj +=
				 i__7) {
			    i__4 = i__ + ilen - 1;
			    i__6 = idepth;
			    for (ii = i__; i__6 < 0 ? ii >= i__4 : ii <= i__4;
				     ii += i__6) {
				t11 = 0.;
				t21 = 0.;
				t12 = 0.;
				t22 = 0.;
				i__8 = l + lspan - 1;
				for (ll = l; ll <= i__8; ++ll) {
				    t11 += ch[ll - l + 1 + (ii - i__ + 1 << 6)
					     - 65] * b[ll + jj * b_dim1];
				    t21 += ch[ll - l + 1 + (ii - i__ + 2 << 6)
					     - 65] * b[ll + jj * b_dim1];
				    t12 += ch[ll - l + 1 + (ii - i__ + 1 << 6)
					     - 65] * b[ll + (jj + 1) * b_dim1]
					    ;
				    t22 += ch[ll - l + 1 + (ii - i__ + 2 << 6)
					     - 65] * b[ll + (jj + 1) * b_dim1]
					    ;
/* L440: */
				}
				c__[ii + jj * c_dim1] += t11;
				c__[ii + 1 + jj * c_dim1] += t21;
				c__[ii + (jj + 1) * c_dim1] += t12;
				c__[ii + 1 + (jj + 1) * c_dim1] += t22;
/* L460: */
			    }
			    if (ilen < ispan) {
				i__6 = i__ + ispan - 1;
				for (ii = i__ + ilen; ii <= i__6; ++ii) {
				    t11 = 0.;
				    t12 = 0.;
				    i__4 = l + lspan - 1;
				    for (ll = l; ll <= i__4; ++ll) {
					t11 += ch[ll - l + 1 + (ii - i__ + 1
						<< 6) - 65] * b[ll + jj *
						b_dim1];
					t12 += ch[ll - l + 1 + (ii - i__ + 1
						<< 6) - 65] * b[ll + (jj + 1)
						* b_dim1];
/* L480: */
				    }
				    c__[ii + jj * c_dim1] += t11;
				    c__[ii + (jj + 1) * c_dim1] += t12;
/* L500: */
				}
			    }
/* L520: */
			}
			if (jlen < jspan) {
			    i__7 = j + jspan - 1;
			    for (jj = j + jlen; jj <= i__7; ++jj) {
				i__5 = i__ + ilen - 1;
				i__6 = idepth;
				for (ii = i__; i__6 < 0 ? ii >= i__5 : ii <=
					i__5; ii += i__6) {
				    t11 = 0.;
				    t21 = 0.;
				    i__4 = l + lspan - 1;
				    for (ll = l; ll <= i__4; ++ll) {
					t11 += ch[ll - l + 1 + (ii - i__ + 1
						<< 6) - 65] * b[ll + jj *
						b_dim1];
					t21 += ch[ll - l + 1 + (ii - i__ + 2
						<< 6) - 65] * b[ll + jj *
						b_dim1];
/* L540: */
				    }
				    c__[ii + jj * c_dim1] += t11;
				    c__[ii + 1 + jj * c_dim1] += t21;
/* L560: */
				}
				if (ilen < ispan) {
				    i__6 = i__ + ispan - 1;
				    for (ii = i__ + ilen; ii <= i__6; ++ii) {
					t11 = 0.;
					i__5 = l + lspan - 1;
					for (ll = l; ll <= i__5; ++ll) {
					    t11 += ch[ll - l + 1 + (ii - i__
						    + 1 << 6) - 65] * b[ll +
						    jj * b_dim1];
/* L580: */
					}
					c__[ii + jj * c_dim1] += t11;
/* L600: */
				    }
				}
/* L620: */
			    }
			}
/* L640: */
		    }
/* L660: */
		}
/* L680: */
	    }
	}
    } else {
	if (nota) {

/*           Form  C := C + alpha*A*B' */

	    i__1 = N;
	    for (j = 1; j <= i__1; j += 64) {
		jdepth = 2;
/* Computing MIN */
		i__2 = 64, i__3 = N - j + 1;
		jspan = min(i__2,i__3);
		jlen = jdepth * (jspan / jdepth);
		i__2 = K;
		for (l = 1; l <= i__2; l += 64) {
/* Computing MIN */
		    i__3 = 64, i__7 = K - l + 1;
		    lspan = min(i__3,i__7);
		    i__3 = j + jspan - 1;
		    for (jj = j; jj <= i__3; ++jj) {
			i__7 = l + lspan - 1;
			for (ll = l; ll <= i__7; ++ll) {
			    ch[ll - l + 1 + (jj - j + 1 << 6) - 65] = alpha *
				     b[jj + ll * b_dim1];
/* L700: */
			}
/* L720: */
		    }
		    i__3 = M;
		    for (i__ = 1; i__ <= i__3; i__ += 64) {
			idepth = 2;
/* Computing MIN */
			i__7 = 64, i__6 = M - i__ + 1;
			ispan = min(i__7,i__6);
			ilen = idepth * (ispan / idepth);
			i__7 = i__ + ilen - 1;
			i__6 = idepth;
			for (ii = i__; i__6 < 0 ? ii >= i__7 : ii <= i__7; ii
				+= i__6) {
			    i__5 = l + lspan - 1;
			    for (ll = l; ll <= i__5; ++ll) {
				ch1[ll - l] = a[ii + ll * a_dim1];
				ch2[ll - l] = a[ii + 1 + ll * a_dim1];
/* L740: */
			    }
			    i__5 = j + jlen - 1;
			    i__4 = jdepth;
			    for (jj = j; i__4 < 0 ? jj >= i__5 : jj <= i__5;
				    jj += i__4) {
				t11 = 0.;
				t21 = 0.;
				t12 = 0.;
				t22 = 0.;
				i__8 = l + lspan - 1;
				for (ll = l; ll <= i__8; ++ll) {
				    t11 += ch1[ll - l] * ch[ll - l + 1 + (jj
					    - j + 1 << 6) - 65];
				    t21 += ch2[ll - l] * ch[ll - l + 1 + (jj
					    - j + 1 << 6) - 65];
				    t12 += ch1[ll - l] * ch[ll - l + 1 + (jj
					    - j + 2 << 6) - 65];
				    t22 += ch2[ll - l] * ch[ll - l + 1 + (jj
					    - j + 2 << 6) - 65];
/* L760: */
				}
				c__[ii + jj * c_dim1] += t11;
				c__[ii + 1 + jj * c_dim1] += t21;
				c__[ii + (jj + 1) * c_dim1] += t12;
				c__[ii + 1 + (jj + 1) * c_dim1] += t22;
/* L780: */
			    }
			    if (jlen < jspan) {
				i__4 = j + jspan - 1;
				for (jj = j + jlen; jj <= i__4; ++jj) {
				    t11 = 0.;
				    t21 = 0.;
				    i__5 = l + lspan - 1;
				    for (ll = l; ll <= i__5; ++ll) {
					t11 += a[ii + ll * a_dim1] * ch[ll -
						l + 1 + (jj - j + 1 << 6) -
						65];
					t21 += a[ii + 1 + ll * a_dim1] * ch[
						ll - l + 1 + (jj - j + 1 << 6)
						 - 65];
/* L800: */
				    }
				    c__[ii + jj * c_dim1] += t11;
				    c__[ii + 1 + jj * c_dim1] += t21;
/* L820: */
				}
			    }
/* L840: */
			}
			if (ilen < ispan) {
			    i__6 = i__ + ispan - 1;
			    for (ii = i__ + ilen; ii <= i__6; ++ii) {
				i__7 = j + jlen - 1;
				i__4 = jdepth;
				for (jj = j; i__4 < 0 ? jj >= i__7 : jj <=
					i__7; jj += i__4) {
				    t11 = 0.;
				    t12 = 0.;
				    i__5 = l + lspan - 1;
				    for (ll = l; ll <= i__5; ++ll) {
					t11 += a[ii + ll * a_dim1] * ch[ll -
						l + 1 + (jj - j + 1 << 6) -
						65];
					t12 += a[ii + ll * a_dim1] * ch[ll -
						l + 1 + (jj - j + 2 << 6) -
						65];
/* L860: */
				    }
				    c__[ii + jj * c_dim1] += t11;
				    c__[ii + (jj + 1) * c_dim1] += t12;
/* L880: */
				}
				if (jlen < jspan) {
				    i__4 = j + jspan - 1;
				    for (jj = j + jlen; jj <= i__4; ++jj) {
					t11 = 0.;
					i__7 = l + lspan - 1;
					for (ll = l; ll <= i__7; ++ll) {
					    t11 += a[ii + ll * a_dim1] * ch[
						    ll - l + 1 + (jj - j + 1
						    << 6) - 65];
/* L900: */
					}
					c__[ii + jj * c_dim1] += t11;
/* L920: */
				    }
				}
/* L940: */
			    }
			}
/* L960: */
		    }
/* L980: */
		}
/* L1000: */
	    }
	} else {

/*           Form  C := C + alpha*A'*B' */

	    i__1 = N;
	    for (j = 1; j <= i__1; j += 64) {
		jdepth = 2;
/* Computing MIN */
		i__2 = 64, i__3 = N - j + 1;
		jspan = min(i__2,i__3);
		jlen = jdepth * (jspan / jdepth);
		i__2 = K;
		for (l = 1; l <= i__2; l += 64) {
/* Computing MIN */
		    i__3 = 64, i__6 = K - l + 1;
		    lspan = min(i__3,i__6);
		    i__3 = j + jspan - 1;
		    for (jj = j; jj <= i__3; ++jj) {
			i__6 = l + lspan - 1;
			for (ll = l; ll <= i__6; ++ll) {
			    ch[ll - l + 1 + (jj - j + 1 << 6) - 65] = alpha *
				     b[jj + ll * b_dim1];
/* L1020: */
			}
/* L1040: */
		    }
		    i__3 = M;
		    for (i__ = 1; i__ <= i__3; i__ += 64) {
			idepth = 2;
/* Computing MIN */
			i__6 = 64, i__4 = M - i__ + 1;
			ispan = min(i__6,i__4);
			ilen = idepth * (ispan / idepth);
			i__6 = i__ + ilen - 1;
			i__4 = idepth;
			for (ii = i__; i__4 < 0 ? ii >= i__6 : ii <= i__6; ii
				+= i__4) {
			    i__7 = j + jlen - 1;
			    i__5 = jdepth;
			    for (jj = j; i__5 < 0 ? jj >= i__7 : jj <= i__7;
				    jj += i__5) {
				t11 = 0.;
				t21 = 0.;
				t12 = 0.;
				t22 = 0.;
				i__8 = l + lspan - 1;
				for (ll = l; ll <= i__8; ++ll) {
				    t11 += a[ll + ii * a_dim1] * ch[ll - l +
					    1 + (jj - j + 1 << 6) - 65];
				    t21 += a[ll + (ii + 1) * a_dim1] * ch[ll
					    - l + 1 + (jj - j + 1 << 6) - 65];
				    t12 += a[ll + ii * a_dim1] * ch[ll - l +
					    1 + (jj - j + 2 << 6) - 65];
				    t22 += a[ll + (ii + 1) * a_dim1] * ch[ll
					    - l + 1 + (jj - j + 2 << 6) - 65];
/* L1060: */
				}
				c__[ii + jj * c_dim1] += t11;
				c__[ii + 1 + jj * c_dim1] += t21;
				c__[ii + (jj + 1) * c_dim1] += t12;
				c__[ii + 1 + (jj + 1) * c_dim1] += t22;
/* L1080: */
			    }
			    if (jlen < jspan) {
				i__5 = j + jspan - 1;
				for (jj = j + jlen; jj <= i__5; ++jj) {
				    t11 = 0.;
				    t21 = 0.;
				    i__7 = l + lspan - 1;
				    for (ll = l; ll <= i__7; ++ll) {
					t11 += a[ll + ii * a_dim1] * ch[ll -
						l + 1 + (jj - j + 1 << 6) -
						65];
					t21 += a[ll + (ii + 1) * a_dim1] * ch[
						ll - l + 1 + (jj - j + 1 << 6)
						 - 65];
/* L1100: */
				    }
				    c__[ii + jj * c_dim1] += t11;
				    c__[ii + 1 + jj * c_dim1] += t21;
/* L1120: */
				}
			    }
/* L1140: */
			}
			if (ilen < ispan) {
			    i__4 = i__ + ispan - 1;
			    for (ii = i__ + ilen; ii <= i__4; ++ii) {
				i__6 = j + jlen - 1;
				i__5 = jdepth;
				for (jj = j; i__5 < 0 ? jj >= i__6 : jj <=
					i__6; jj += i__5) {
				    t11 = 0.;
				    t12 = 0.;
				    i__7 = l + lspan - 1;
				    for (ll = l; ll <= i__7; ++ll) {
					t11 += a[ll + ii * a_dim1] * ch[ll -
						l + 1 + (jj - j + 1 << 6) -
						65];
					t12 += a[ll + ii * a_dim1] * ch[ll -
						l + 1 + (jj - j + 2 << 6) -
						65];
/* L1160: */
				    }
				    c__[ii + jj * c_dim1] += t11;
				    c__[ii + (jj + 1) * c_dim1] += t12;
/* L1180: */
				}
				if (jlen < jspan) {
				    i__5 = j + jspan - 1;
				    for (jj = j + jlen; jj <= i__5; ++jj) {
					t11 = 0.;
					i__6 = l + lspan - 1;
					for (ll = l; ll <= i__6; ++ll) {
					    t11 += a[ll + ii * a_dim1] * ch[
						    ll - l + 1 + (jj - j + 1
						    << 6) - 65];
/* L1200: */
					}
					c__[ii + jj * c_dim1] += t11;
/* L1220: */
				    }
				}
/* L1240: */
			    }
			}
/* L1260: */
		    }
/* L1280: */
		}
/* L1300: */
	    }
	}
    }

    return 0;

/*     End of DGEMM . */

} /* dgemm_ */




//Error Message
int _blError(char *srname, int *info)
{
    printf("** On entry to %6s, parameter number %2i had an illegal value\n",
		srname, *info);
    return 0;
}

#endif
