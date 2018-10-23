/*
**
** PHiPAC Matrix-Matrix Code for the operation:
**    C = transpose(A)*B + beta*C
**
** Automatically Generated by mm_gen ($Revision: 1.32 $) using the command:
**    ./../mm_gen/mm_gen -cb 16 1 2 -cb 1 32 32 -prec complex -file mdtmd_a1bc_md.c -routine_name mdtmd_a1bc_md -beta c -opA T
**
** Run './../mm_gen/mm_gen -help' for help.
**
** Generated on: Thursday March 21 2002, 19:48:10 PST
** Created by: Jeff Bilmes <bilmes@cs.berkeley.edu>
**             http://www.icsi.berkeley.edu/~bilmes/phipac
**
**
** Usage:
**    mdtmd_a1bc_md(const int M, const int K, const int N, const complex *const A, const complex *const B, complex *const C, const int Astride, const int Bstride, const int Cstride, const complex beta)
** where
**  transpose(A) is an MxK matrix
**  B is an KxN matrix
**  C is an MxN matrix
**  Astride is the number of entries between the start of each row of A
**  Bstride is the number of entries between the start of each row of B
**  Cstride is the number of entries between the start of each row of C
**
**
** "Copyright (c) 1995 The Regents of the University of California.  All
** rights reserved."  Permission to use, copy, modify, and distribute
** this software and its documentation for any purpose, without fee, and
** without written agreement is hereby granted, provided that the above
** copyright notice and the following two paragraphs appear in all copies
** of this software.
**
** IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
** DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
** OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
** CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
** THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
** INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
** AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
** ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO
** PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
**
*/

#include "container/complex.h"
#include "container/matrix/mdmd_cm.h"


BEGIN_BL_NAMESPACE


#define LOAD1x1(c00,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; \
}

#define STORE1x1(c00,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; \
}

/* Fixed M,K,N = 1,1,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd1x1md1x1_md1x1(c00,A,Astride,B,Bstride) \
{ \
    complex _b0; \
    complex _a0; \
   \
   \
   _b0 = B[0]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; \
   A += Astride; \
}


#define LOAD1x2(c00,c01,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; \
}

#define STORE1x2(c00,c01,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; \
}

/* Fixed M,K,N = 1,1,2 fully-unrolled matrix matrix multiply. */
#define mul_tmd1x1md1x2_md1x2(c00,c01,A,Astride,B,Bstride) \
{ \
    complex _b0,_b1; \
    complex _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; \
   A += Astride; \
}


#define LOAD2x1(c00,c10,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; \
   _cp += Cstride; \
   c10 = _cp[0]; \
}

#define STORE2x1(c00,c10,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; \
   _cp += Cstride; \
   _cp[0] = c10; \
}

/* Fixed M,K,N = 2,1,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd2x1md1x1_md2x1(c00,c10,A,Astride,B,Bstride) \
{ \
    complex _b0; \
    complex _a0,_a1; \
   \
   \
   _b0 = B[0]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; \
   _a1 = A[1]; \
   c10 += _a1*_b0; \
   A += Astride; \
}


#define LOAD2x2(c00,c01,c10,c11,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; \
   _cp += Cstride; \
   c10 = _cp[0]; c11 = _cp[1]; \
}

#define STORE2x2(c00,c01,c10,c11,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; \
   _cp += Cstride; \
   _cp[0] = c10; _cp[1] = c11; \
}

/* Fixed M,K,N = 2,1,2 fully-unrolled matrix matrix multiply. */
#define mul_tmd2x1md1x2_md2x2(c00,c01,c10,c11,A,Astride,B,Bstride) \
{ \
    complex _b0,_b1; \
    complex _a0,_a1; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; \
   _a1 = A[1]; \
   c10 += _a1*_b0; c11 += _a1*_b1; \
   A += Astride; \
}


#define LOAD4x1(c00,c10,c20,c30,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; \
   _cp += Cstride; \
   c10 = _cp[0]; \
   _cp += Cstride; \
   c20 = _cp[0]; \
   _cp += Cstride; \
   c30 = _cp[0]; \
}

#define STORE4x1(c00,c10,c20,c30,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; \
   _cp += Cstride; \
   _cp[0] = c10; \
   _cp += Cstride; \
   _cp[0] = c20; \
   _cp += Cstride; \
   _cp[0] = c30; \
}

/* Fixed M,K,N = 4,1,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd4x1md1x1_md4x1(c00,c10,c20,c30,A,Astride,B,Bstride) \
{ \
    complex _b0; \
    complex _a0,_a1,_a2,_a3; \
   \
   \
   _b0 = B[0]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; \
   _a1 = A[1]; \
   c10 += _a1*_b0; \
   _a2 = A[2]; \
   c20 += _a2*_b0; \
   _a3 = A[3]; \
   c30 += _a3*_b0; \
   A += Astride; \
}


#define LOAD4x2(c00,c01,c10,c11,c20,c21,c30,c31,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; \
   _cp += Cstride; \
   c10 = _cp[0]; c11 = _cp[1]; \
   _cp += Cstride; \
   c20 = _cp[0]; c21 = _cp[1]; \
   _cp += Cstride; \
   c30 = _cp[0]; c31 = _cp[1]; \
}

#define STORE4x2(c00,c01,c10,c11,c20,c21,c30,c31,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; \
   _cp += Cstride; \
   _cp[0] = c10; _cp[1] = c11; \
   _cp += Cstride; \
   _cp[0] = c20; _cp[1] = c21; \
   _cp += Cstride; \
   _cp[0] = c30; _cp[1] = c31; \
}

/* Fixed M,K,N = 4,1,2 fully-unrolled matrix matrix multiply. */
#define mul_tmd4x1md1x2_md4x2(c00,c01,c10,c11,c20,c21,c30,c31,A,Astride,B,Bstride) \
{ \
    complex _b0,_b1; \
    complex _a0,_a1,_a2,_a3; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; \
   _a1 = A[1]; \
   c10 += _a1*_b0; c11 += _a1*_b1; \
   _a2 = A[2]; \
   c20 += _a2*_b0; c21 += _a2*_b1; \
   _a3 = A[3]; \
   c30 += _a3*_b0; c31 += _a3*_b1; \
   A += Astride; \
}


#define LOAD8x1(c00,c10,c20,c30,c40,c50,c60,c70,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; \
   _cp += Cstride; \
   c10 = _cp[0]; \
   _cp += Cstride; \
   c20 = _cp[0]; \
   _cp += Cstride; \
   c30 = _cp[0]; \
   _cp += Cstride; \
   c40 = _cp[0]; \
   _cp += Cstride; \
   c50 = _cp[0]; \
   _cp += Cstride; \
   c60 = _cp[0]; \
   _cp += Cstride; \
   c70 = _cp[0]; \
}

#define STORE8x1(c00,c10,c20,c30,c40,c50,c60,c70,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; \
   _cp += Cstride; \
   _cp[0] = c10; \
   _cp += Cstride; \
   _cp[0] = c20; \
   _cp += Cstride; \
   _cp[0] = c30; \
   _cp += Cstride; \
   _cp[0] = c40; \
   _cp += Cstride; \
   _cp[0] = c50; \
   _cp += Cstride; \
   _cp[0] = c60; \
   _cp += Cstride; \
   _cp[0] = c70; \
}

/* Fixed M,K,N = 8,1,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd8x1md1x1_md8x1(c00,c10,c20,c30,c40,c50,c60,c70,A,Astride,B,Bstride) \
{ \
    complex _b0; \
    complex _a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7; \
   \
   \
   _b0 = B[0]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; \
   _a1 = A[1]; \
   c10 += _a1*_b0; \
   _a2 = A[2]; \
   c20 += _a2*_b0; \
   _a3 = A[3]; \
   c30 += _a3*_b0; \
   _a4 = A[4]; \
   c40 += _a4*_b0; \
   _a5 = A[5]; \
   c50 += _a5*_b0; \
   _a6 = A[6]; \
   c60 += _a6*_b0; \
   _a7 = A[7]; \
   c70 += _a7*_b0; \
   A += Astride; \
}


#define LOAD8x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; \
   _cp += Cstride; \
   c10 = _cp[0]; c11 = _cp[1]; \
   _cp += Cstride; \
   c20 = _cp[0]; c21 = _cp[1]; \
   _cp += Cstride; \
   c30 = _cp[0]; c31 = _cp[1]; \
   _cp += Cstride; \
   c40 = _cp[0]; c41 = _cp[1]; \
   _cp += Cstride; \
   c50 = _cp[0]; c51 = _cp[1]; \
   _cp += Cstride; \
   c60 = _cp[0]; c61 = _cp[1]; \
   _cp += Cstride; \
   c70 = _cp[0]; c71 = _cp[1]; \
}

#define STORE8x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; \
   _cp += Cstride; \
   _cp[0] = c10; _cp[1] = c11; \
   _cp += Cstride; \
   _cp[0] = c20; _cp[1] = c21; \
   _cp += Cstride; \
   _cp[0] = c30; _cp[1] = c31; \
   _cp += Cstride; \
   _cp[0] = c40; _cp[1] = c41; \
   _cp += Cstride; \
   _cp[0] = c50; _cp[1] = c51; \
   _cp += Cstride; \
   _cp[0] = c60; _cp[1] = c61; \
   _cp += Cstride; \
   _cp[0] = c70; _cp[1] = c71; \
}

/* Fixed M,K,N = 8,1,2 fully-unrolled matrix matrix multiply. */
#define mul_tmd8x1md1x2_md8x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,A,Astride,B,Bstride) \
{ \
    complex _b0,_b1; \
    complex _a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; \
   _a1 = A[1]; \
   c10 += _a1*_b0; c11 += _a1*_b1; \
   _a2 = A[2]; \
   c20 += _a2*_b0; c21 += _a2*_b1; \
   _a3 = A[3]; \
   c30 += _a3*_b0; c31 += _a3*_b1; \
   _a4 = A[4]; \
   c40 += _a4*_b0; c41 += _a4*_b1; \
   _a5 = A[5]; \
   c50 += _a5*_b0; c51 += _a5*_b1; \
   _a6 = A[6]; \
   c60 += _a6*_b0; c61 += _a6*_b1; \
   _a7 = A[7]; \
   c70 += _a7*_b0; c71 += _a7*_b1; \
   A += Astride; \
}


#define LOAD16x1(c00,c10,c20,c30,c40,c50,c60,c70,c80,c90,c100,c110,c120,c130,c140,c150,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; \
   _cp += Cstride; \
   c10 = _cp[0]; \
   _cp += Cstride; \
   c20 = _cp[0]; \
   _cp += Cstride; \
   c30 = _cp[0]; \
   _cp += Cstride; \
   c40 = _cp[0]; \
   _cp += Cstride; \
   c50 = _cp[0]; \
   _cp += Cstride; \
   c60 = _cp[0]; \
   _cp += Cstride; \
   c70 = _cp[0]; \
   _cp += Cstride; \
   c80 = _cp[0]; \
   _cp += Cstride; \
   c90 = _cp[0]; \
   _cp += Cstride; \
   c100 = _cp[0]; \
   _cp += Cstride; \
   c110 = _cp[0]; \
   _cp += Cstride; \
   c120 = _cp[0]; \
   _cp += Cstride; \
   c130 = _cp[0]; \
   _cp += Cstride; \
   c140 = _cp[0]; \
   _cp += Cstride; \
   c150 = _cp[0]; \
}

#define STORE16x1(c00,c10,c20,c30,c40,c50,c60,c70,c80,c90,c100,c110,c120,c130,c140,c150,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; \
   _cp += Cstride; \
   _cp[0] = c10; \
   _cp += Cstride; \
   _cp[0] = c20; \
   _cp += Cstride; \
   _cp[0] = c30; \
   _cp += Cstride; \
   _cp[0] = c40; \
   _cp += Cstride; \
   _cp[0] = c50; \
   _cp += Cstride; \
   _cp[0] = c60; \
   _cp += Cstride; \
   _cp[0] = c70; \
   _cp += Cstride; \
   _cp[0] = c80; \
   _cp += Cstride; \
   _cp[0] = c90; \
   _cp += Cstride; \
   _cp[0] = c100; \
   _cp += Cstride; \
   _cp[0] = c110; \
   _cp += Cstride; \
   _cp[0] = c120; \
   _cp += Cstride; \
   _cp[0] = c130; \
   _cp += Cstride; \
   _cp[0] = c140; \
   _cp += Cstride; \
   _cp[0] = c150; \
}

/* Fixed M,K,N = 16,1,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd16x1md1x1_md16x1(c00,c10,c20,c30,c40,c50,c60,c70,c80,c90,c100,c110,c120,c130,c140,c150,A,Astride,B,Bstride) \
{ \
    complex _b0; \
    complex _a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7,_a8,_a9,_a10,_a11,_a12,_a13,_a14,_a15; \
   \
   \
   _b0 = B[0]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; \
   _a1 = A[1]; \
   c10 += _a1*_b0; \
   _a2 = A[2]; \
   c20 += _a2*_b0; \
   _a3 = A[3]; \
   c30 += _a3*_b0; \
   _a4 = A[4]; \
   c40 += _a4*_b0; \
   _a5 = A[5]; \
   c50 += _a5*_b0; \
   _a6 = A[6]; \
   c60 += _a6*_b0; \
   _a7 = A[7]; \
   c70 += _a7*_b0; \
   _a8 = A[8]; \
   c80 += _a8*_b0; \
   _a9 = A[9]; \
   c90 += _a9*_b0; \
   _a10 = A[10]; \
   c100 += _a10*_b0; \
   _a11 = A[11]; \
   c110 += _a11*_b0; \
   _a12 = A[12]; \
   c120 += _a12*_b0; \
   _a13 = A[13]; \
   c130 += _a13*_b0; \
   _a14 = A[14]; \
   c140 += _a14*_b0; \
   _a15 = A[15]; \
   c150 += _a15*_b0; \
   A += Astride; \
}


#define LOAD16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,C,Cstride) \
{\
   complex * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; \
   _cp += Cstride; \
   c10 = _cp[0]; c11 = _cp[1]; \
   _cp += Cstride; \
   c20 = _cp[0]; c21 = _cp[1]; \
   _cp += Cstride; \
   c30 = _cp[0]; c31 = _cp[1]; \
   _cp += Cstride; \
   c40 = _cp[0]; c41 = _cp[1]; \
   _cp += Cstride; \
   c50 = _cp[0]; c51 = _cp[1]; \
   _cp += Cstride; \
   c60 = _cp[0]; c61 = _cp[1]; \
   _cp += Cstride; \
   c70 = _cp[0]; c71 = _cp[1]; \
   _cp += Cstride; \
   c80 = _cp[0]; c81 = _cp[1]; \
   _cp += Cstride; \
   c90 = _cp[0]; c91 = _cp[1]; \
   _cp += Cstride; \
   c100 = _cp[0]; c101 = _cp[1]; \
   _cp += Cstride; \
   c110 = _cp[0]; c111 = _cp[1]; \
   _cp += Cstride; \
   c120 = _cp[0]; c121 = _cp[1]; \
   _cp += Cstride; \
   c130 = _cp[0]; c131 = _cp[1]; \
   _cp += Cstride; \
   c140 = _cp[0]; c141 = _cp[1]; \
   _cp += Cstride; \
   c150 = _cp[0]; c151 = _cp[1]; \
}

#define STORE16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,C,Cstride) \
{\
   complex *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; \
   _cp += Cstride; \
   _cp[0] = c10; _cp[1] = c11; \
   _cp += Cstride; \
   _cp[0] = c20; _cp[1] = c21; \
   _cp += Cstride; \
   _cp[0] = c30; _cp[1] = c31; \
   _cp += Cstride; \
   _cp[0] = c40; _cp[1] = c41; \
   _cp += Cstride; \
   _cp[0] = c50; _cp[1] = c51; \
   _cp += Cstride; \
   _cp[0] = c60; _cp[1] = c61; \
   _cp += Cstride; \
   _cp[0] = c70; _cp[1] = c71; \
   _cp += Cstride; \
   _cp[0] = c80; _cp[1] = c81; \
   _cp += Cstride; \
   _cp[0] = c90; _cp[1] = c91; \
   _cp += Cstride; \
   _cp[0] = c100; _cp[1] = c101; \
   _cp += Cstride; \
   _cp[0] = c110; _cp[1] = c111; \
   _cp += Cstride; \
   _cp[0] = c120; _cp[1] = c121; \
   _cp += Cstride; \
   _cp[0] = c130; _cp[1] = c131; \
   _cp += Cstride; \
   _cp[0] = c140; _cp[1] = c141; \
   _cp += Cstride; \
   _cp[0] = c150; _cp[1] = c151; \
}

/* Fixed M,K,N = 16,1,2 fully-unrolled matrix matrix multiply. */
#define mul_tmd16x1md1x2_md16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,A,Astride,B,Bstride) \
{ \
    complex _b0,_b1; \
    complex _a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7,_a8,_a9,_a10,_a11,_a12,_a13,_a14,_a15; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; \
   B += Bstride; \
   _a0 = A[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; \
   _a1 = A[1]; \
   c10 += _a1*_b0; c11 += _a1*_b1; \
   _a2 = A[2]; \
   c20 += _a2*_b0; c21 += _a2*_b1; \
   _a3 = A[3]; \
   c30 += _a3*_b0; c31 += _a3*_b1; \
   _a4 = A[4]; \
   c40 += _a4*_b0; c41 += _a4*_b1; \
   _a5 = A[5]; \
   c50 += _a5*_b0; c51 += _a5*_b1; \
   _a6 = A[6]; \
   c60 += _a6*_b0; c61 += _a6*_b1; \
   _a7 = A[7]; \
   c70 += _a7*_b0; c71 += _a7*_b1; \
   _a8 = A[8]; \
   c80 += _a8*_b0; c81 += _a8*_b1; \
   _a9 = A[9]; \
   c90 += _a9*_b0; c91 += _a9*_b1; \
   _a10 = A[10]; \
   c100 += _a10*_b0; c101 += _a10*_b1; \
   _a11 = A[11]; \
   c110 += _a11*_b0; c111 += _a11*_b1; \
   _a12 = A[12]; \
   c120 += _a12*_b0; c121 += _a12*_b1; \
   _a13 = A[13]; \
   c130 += _a13*_b0; c131 += _a13*_b1; \
   _a14 = A[14]; \
   c140 += _a14*_b0; c141 += _a14*_b1; \
   _a15 = A[15]; \
   c150 += _a15*_b0; c151 += _a15*_b1; \
   A += Astride; \
}


/* Fixed M,N = 16,64, Arbitrary K L0-blocked matrix matrix multiply. */
static void
mdtmd_a1bc_md_l1_arb_k(int K, const complex *const A, const complex *const B, complex *const C, const int Astride, const int Bstride, const int Cstride)
{
   const complex *a0,*b0;
   complex *c0;
   const complex *ap0;
   const complex *bp0;
   complex *cp0;
   const int C_sbs_stride = Cstride*16;
   const int k_marg_el = K & 0;
   const int k_norm = (K - k_marg_el)*Astride;
   complex *const c0_endp = C+16*Cstride;
    complex c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=16) {
      const complex* const ap0_endp = a0 + k_norm;
      complex* const cp0_endp = c0 + 64;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=2,cp0+=2) {
         ap0=a0;
         bp0=b0;
         LOAD16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd16x1md1x2_md16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,ap0,Astride,bp0,Bstride);
         }
         STORE16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,cp0,Cstride);
      }
   }
}

/* Arbitrary M,K,N L0-blocked matrix matrix multiply. */
static void
mdtmd_a1bc_md_l1_arb_all(const int M, const int K, const int N, const complex *const A, const complex *const B, complex *const C, const int Astride, const int Bstride, const int Cstride)
{
   const complex *a0,*b0;
   complex *c0;
   const complex *ap0;
   const complex *bp0;
   complex *cp0;
   const int C_sbs_stride = Cstride*16;
   const int k_marg_el = K & 0;
   const int k_norm = (K - k_marg_el)*Astride;
   const int m_marg_el = M & 15;
   const int m_norm = M - m_marg_el;
   const int n_marg_el = N & 1;
   const int n_norm = N - n_marg_el;
   complex *const c0_endp = C+m_norm*Cstride;
    complex c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=16) {
      const complex* const ap0_endp = a0 + k_norm;
      complex* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=2,cp0+=2) {
         ap0=a0;
         bp0=b0;
         LOAD16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd16x1md1x2_md16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,ap0,Astride,bp0,Bstride);
         }
         STORE16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,cp0,Cstride);
      }
   }
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=16) {
      const complex* const ap0_endp = a0 + k_norm;
      b0 = B+n_norm;
      cp0 = c0+n_norm;
      if (n_marg_el & 0x1) {
         ap0=a0;
         bp0=b0;
         LOAD16x1(c00,c10,c20,c30,c40,c50,c60,c70,c80,c90,c100,c110,c120,c130,c140,c150,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd16x1md1x1_md16x1(c00,c10,c20,c30,c40,c50,c60,c70,c80,c90,c100,c110,c120,c130,c140,c150,ap0,Astride,bp0,Bstride);
         }
         STORE16x1(c00,c10,c20,c30,c40,c50,c60,c70,c80,c90,c100,c110,c120,c130,c140,c150,cp0,Cstride);
      }
   }
   if (m_marg_el & 0x8) {
      const complex* const ap0_endp = a0 + k_norm;
      complex* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=2,cp0+=2) {
         ap0=a0;
         bp0=b0;
         LOAD8x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd8x1md1x2_md8x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,ap0,Astride,bp0,Bstride);
         }
         STORE8x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,cp0,Cstride);
      }
      if (n_marg_el & 0x1) {
         ap0=a0;
         bp0=b0;
         LOAD8x1(c00,c10,c20,c30,c40,c50,c60,c70,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd8x1md1x1_md8x1(c00,c10,c20,c30,c40,c50,c60,c70,ap0,Astride,bp0,Bstride);
         }
         STORE8x1(c00,c10,c20,c30,c40,c50,c60,c70,cp0,Cstride);
      }
      c0 += Cstride*8;
      a0 += 8;
   }
   if (m_marg_el & 0x4) {
      const complex* const ap0_endp = a0 + k_norm;
      complex* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=2,cp0+=2) {
         ap0=a0;
         bp0=b0;
         LOAD4x2(c00,c01,c10,c11,c20,c21,c30,c31,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd4x1md1x2_md4x2(c00,c01,c10,c11,c20,c21,c30,c31,ap0,Astride,bp0,Bstride);
         }
         STORE4x2(c00,c01,c10,c11,c20,c21,c30,c31,cp0,Cstride);
      }
      if (n_marg_el & 0x1) {
         ap0=a0;
         bp0=b0;
         LOAD4x1(c00,c10,c20,c30,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd4x1md1x1_md4x1(c00,c10,c20,c30,ap0,Astride,bp0,Bstride);
         }
         STORE4x1(c00,c10,c20,c30,cp0,Cstride);
      }
      c0 += Cstride*4;
      a0 += 4;
   }
   if (m_marg_el & 0x2) {
      const complex* const ap0_endp = a0 + k_norm;
      complex* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=2,cp0+=2) {
         ap0=a0;
         bp0=b0;
         LOAD2x2(c00,c01,c10,c11,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd2x1md1x2_md2x2(c00,c01,c10,c11,ap0,Astride,bp0,Bstride);
         }
         STORE2x2(c00,c01,c10,c11,cp0,Cstride);
      }
      if (n_marg_el & 0x1) {
         ap0=a0;
         bp0=b0;
         LOAD2x1(c00,c10,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd2x1md1x1_md2x1(c00,c10,ap0,Astride,bp0,Bstride);
         }
         STORE2x1(c00,c10,cp0,Cstride);
      }
      c0 += Cstride*2;
      a0 += 2;
   }
   if (m_marg_el & 0x1) {
      const complex* const ap0_endp = a0 + k_norm;
      complex* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=2,cp0+=2) {
         ap0=a0;
         bp0=b0;
         LOAD1x2(c00,c01,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd1x1md1x2_md1x2(c00,c01,ap0,Astride,bp0,Bstride);
         }
         STORE1x2(c00,c01,cp0,Cstride);
      }
      if (n_marg_el & 0x1) {
         ap0=a0;
         bp0=b0;
         LOAD1x1(c00,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd1x1md1x1_md1x1(c00,ap0,Astride,bp0,Bstride);
         }
         STORE1x1(c00,cp0,Cstride);
      }
   }
}

/* Fixed M,K,N = 16,32,64 L0-blocked matrix matrix multiply. */
static void
mdtmd_a1bc_md_l1(const complex *const A, const complex *const B, complex *const C, const int Astride, const int Bstride, const int Cstride)
{
   const complex *a0,*b0;
   complex *c0;
   const complex *ap0;
   const complex *bp0;
   complex *cp0;
   const int C_sbs_stride = Cstride*16;
   const int k_norm = 32*Astride;
   complex *const c0_endp = C+16*Cstride;
    complex c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=16) {
      const complex* const ap0_endp = a0 + k_norm;
      complex* const cp0_endp = c0 + 64;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=2,cp0+=2) {
         ap0=a0;
         bp0=b0;
         LOAD16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,cp0,Cstride);
         for (; ap0!=ap0_endp; ) {
            mul_tmd16x1md1x2_md16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,ap0,Astride,bp0,Bstride);
         }
         STORE16x2(c00,c01,c10,c11,c20,c21,c30,c31,c40,c41,c50,c51,c60,c61,c70,c71,c80,c81,c90,c91,c100,c101,c110,c111,c120,c121,c130,c131,c140,c141,c150,c151,cp0,Cstride);
      }
   }
}

void
mdtmd_a1bc_md_cm(const int M, const int K, const int N, const complex *const A, const complex *const B, complex *const C, const int Astride, const int Bstride, const int Cstride, const complex beta)
{
   /* Code for L1-blocked routine. */
   int m2,k2,n2;
   const complex *a2,*b2;
   complex *c2;
   const complex *ap2,*bp2;
   complex *cp2;
   {
      complex *cprb,*cpre,*cp,*cpe;
      cpre = C + M*Cstride;
      for (cprb = C; cprb != cpre; cprb += Cstride) {
         cpe = cprb + N;
         for (cp = cprb; cp != cpe; cp++) {
            *cp *= (beta);
         }
      }
   }
   if (M < 17 && K < 33 && N < 65) {
      mdtmd_a1bc_md_l1_arb_all(M,K,N,A,B,C,Astride,Bstride,Cstride);
      return;
   }
   for (m2=0; m2<=M-16; m2+=16) {
      c2 = C + m2*Cstride;
      a2 = A + m2;
      for (n2=0,b2=B,cp2=c2; n2<=N-64; n2+=64,b2+=64,cp2+=64) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-32; k2+=32,bp2+=32*Bstride,ap2+=32*Astride) {
            mdtmd_a1bc_md_l1(ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmd_a1bc_md_l1_arb_k(K-k2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
      if (n2 < N) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-32; k2+=32,bp2+=32*Bstride,ap2+=32*Astride) {
            mdtmd_a1bc_md_l1_arb_all(16,32,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmd_a1bc_md_l1_arb_all(16,K-k2,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
   }
   if (m2 < M) {
      c2 = C + m2*Cstride;
      a2 = A + m2;
      for (n2=0,b2=B,cp2=c2; n2<=N-64; n2+=64,b2+=64,cp2+=64) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-32; k2+=32,bp2+=32*Bstride,ap2+=32*Astride) {
            mdtmd_a1bc_md_l1_arb_all(M-m2,32,64,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmd_a1bc_md_l1_arb_all(M-m2,K-k2,64,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
      if (n2 < N) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-32; k2+=32,bp2+=32*Bstride,ap2+=32*Astride) {
            mdtmd_a1bc_md_l1_arb_all(M-m2,32,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmd_a1bc_md_l1_arb_all(M-m2,K-k2,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
   }
}

END_BL_NAMESPACE
