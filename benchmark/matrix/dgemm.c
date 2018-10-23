/* dgemm.f -- translated by f2c (version 20020314).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static int c__1 = 1;

/* Subroutine */ int dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb,
	beta, c__, ldc, transa_len, transb_len)
char *transa, *transb;
int *m, *n, *k;
T *alpha, *a;
int *lda;
T *b;
int *ldb;
T *beta, *c__;
int *ldc;
char transa_len;
char transb_len;
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
	    i__3, i__4, i__5, i__6, i__7, i__8;

    /* Local variables */
    static int ilen, jlen, info;
    static logical nota, notb;
    static int i__, j, l, ncola;
    extern logical lsame_();
    static int ispan, jspan, lspan, nrowa, nrowb;
    static T ch[4096]	/* was [64][64] */;
    static int ii, jj;
    static T t11, t12, t21, t22;
    static int ll;
    static int idepth, jdepth;
    static T ch1[64], ch2[64];

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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;

    /* Function Body */
    nota = lsame_(transa, "N", (char)1, (char)1);
    notb = lsame_(transb, "N", (char)1, (char)1);
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && ! lsame_(transa, "C", (char)1, (char)1) && ! lsame_(
	    transa, "T", (char)1, (char)1)) {
	info = 1;
    } else if (! notb && ! lsame_(transb, "C", (char)1, (char)1) && !
	    lsame_(transb, "T", (char)1, (char)1)) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < max(1,nrowa)) {
	info = 8;
    } else if (*ldb < max(1,nrowb)) {
	info = 10;
    } else if (*ldc < max(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("DGEMM ", &info, (char)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }
    if (*beta == 0.) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		c__[i__ + j * c_dim1] = 0.;
/* L20: */
	    }
/* L40: */
	}
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L60: */
	    }
/* L80: */
	}
    }

/*     And if  alpha.eq.zero. */

    if (*alpha == 0.) {
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := C + alpha*A*B. */

	    i__1 = *k;
	    for (l = 1; l <= i__1; l += 64) {
/* Computing MIN */
		i__2 = 64, i__3 = *k - l + 1;
		lspan = min(i__2,i__3);
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; i__ += 64) {
		    idepth = 2;
/* Computing MIN */
		    i__3 = 64, i__4 = *m - i__ + 1;
		    ispan = min(i__3,i__4);
		    ilen = idepth * (ispan / idepth);
		    i__3 = i__ + ispan - 1;
		    for (ii = i__; ii <= i__3; ++ii) {
			i__4 = l + lspan - 1;
			for (ll = l; ll <= i__4; ++ll) {
			    ch[ll - l + 1 + (ii - i__ + 1 << 6) - 65] = *
				    alpha * a[ii + ll * a_dim1];
/* L100: */
			}
/* L120: */
		    }
		    i__3 = *n;
		    for (j = 1; j <= i__3; j += 64) {
			jdepth = 2;
/* Computing MIN */
			i__4 = 64, i__5 = *n - j + 1;
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

	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; i__ += 64) {
		idepth = 2;
/* Computing MIN */
		i__2 = 64, i__3 = *m - i__ + 1;
		ispan = min(i__2,i__3);
		ilen = idepth * (ispan / idepth);
		i__2 = *k;
		for (l = 1; l <= i__2; l += 64) {
/* Computing MIN */
		    i__3 = 64, i__5 = *k - l + 1;
		    lspan = min(i__3,i__5);
		    i__3 = i__ + ispan - 1;
		    for (ii = i__; ii <= i__3; ++ii) {
			i__5 = l + lspan - 1;
			for (ll = l; ll <= i__5; ++ll) {
			    ch[ll - l + 1 + (ii - i__ + 1 << 6) - 65] = *
				    alpha * a[ll + ii * a_dim1];
/* L400: */
			}
/* L420: */
		    }
		    i__3 = *n;
		    for (j = 1; j <= i__3; j += 64) {
			jdepth = 2;
/* Computing MIN */
			i__5 = 64, i__7 = *n - j + 1;
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

	    i__1 = *n;
	    for (j = 1; j <= i__1; j += 64) {
		jdepth = 2;
/* Computing MIN */
		i__2 = 64, i__3 = *n - j + 1;
		jspan = min(i__2,i__3);
		jlen = jdepth * (jspan / jdepth);
		i__2 = *k;
		for (l = 1; l <= i__2; l += 64) {
/* Computing MIN */
		    i__3 = 64, i__7 = *k - l + 1;
		    lspan = min(i__3,i__7);
		    i__3 = j + jspan - 1;
		    for (jj = j; jj <= i__3; ++jj) {
			i__7 = l + lspan - 1;
			for (ll = l; ll <= i__7; ++ll) {
			    ch[ll - l + 1 + (jj - j + 1 << 6) - 65] = *alpha *
				     b[jj + ll * b_dim1];
/* L700: */
			}
/* L720: */
		    }
		    i__3 = *m;
		    for (i__ = 1; i__ <= i__3; i__ += 64) {
			idepth = 2;
/* Computing MIN */
			i__7 = 64, i__6 = *m - i__ + 1;
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

	    i__1 = *n;
	    for (j = 1; j <= i__1; j += 64) {
		jdepth = 2;
/* Computing MIN */
		i__2 = 64, i__3 = *n - j + 1;
		jspan = min(i__2,i__3);
		jlen = jdepth * (jspan / jdepth);
		i__2 = *k;
		for (l = 1; l <= i__2; l += 64) {
/* Computing MIN */
		    i__3 = 64, i__6 = *k - l + 1;
		    lspan = min(i__3,i__6);
		    i__3 = j + jspan - 1;
		    for (jj = j; jj <= i__3; ++jj) {
			i__6 = l + lspan - 1;
			for (ll = l; ll <= i__6; ++ll) {
			    ch[ll - l + 1 + (jj - j + 1 << 6) - 65] = *alpha *
				     b[jj + ll * b_dim1];
/* L1020: */
			}
/* L1040: */
		    }
		    i__3 = *m;
		    for (i__ = 1; i__ <= i__3; i__ += 64) {
			idepth = 2;
/* Computing MIN */
			i__6 = 64, i__4 = *m - i__ + 1;
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

logical lsame_(ca, cb, ca_len, cb_len)
char *ca, *cb;
char ca_len;
char cb_len;
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static int inta, intb, zcode;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     January 31, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/*  case. */

/*  Arguments */
/*  ========= */

/*  CA      (input) CHARACTER*1 */
/*  CB      (input) CHARACTER*1 */
/*          CA and CB specify the single characters to be compared. */

/* ===================================================================== */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test if the characters are equal */

    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower or */
/*        upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or */
/*        upper case 'Z'. */

	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta
		>= 162 && inta <= 169) {
	    inta += 64;
	}
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb
		>= 162 && intb <= 169) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code */
/*        plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

    return ret_val;
} /* lsame_ */

/* Subroutine */ int xerbla_(srname, info, srname_len)
char *srname;
int *info;
char srname_len;
{
    /* Format strings */
    static char fmt_9999[] = "(\002 ** On entry to \002,a6,\002 parameter nu\
mber \002,i2,\002 had \002,\002an illegal value\002)";

    /* Builtin functions */
    int s_wsfe(), do_fio(), e_wsfe();
    /* Subroutine */ int s_stop();

    /* Fortran I/O blocks */
    static cilist io___30 = { 0, 6, 0, fmt_9999, 0 };



/*  -- LAPACK auxiliary routine (preliminary version) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  XERBLA  is an error handler for the LAPACK routines. */
/*  It is called by an LAPACK routine if an input parameter has an */
/*  invalid value.  A message is printed and execution stops. */

/*  Installers may consider modifying the STOP statement in order to */
/*  call system-specific exception-handling facilities. */

/*  Arguments */
/*  ========= */

/*  SRNAME  (input) CHARACTER*6 */
/*          The name of the routine which called XERBLA. */

/*  INFO    (input) INTEGER */
/*          The position of the invalid parameter in the parameter list */
/*          of the calling routine. */


    s_wsfe(&io___30);
    do_fio(&c__1, srname, (char)6);
    do_fio(&c__1, (char *)&(*info), (char)sizeof(int));
    e_wsfe();

    s_stop("", (char)0);


/*     End of XERBLA */

} /* xerbla_ */

