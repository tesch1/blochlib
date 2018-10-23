/*
**
** PHiPAC Matrix-Matrix Code for the operation:
**    C = transpose(A)*transpose(B) + beta*C
**
** Automatically Generated by mm_gen ($Revision: 1.32 $) using the command:
**    ./../mm_gen/mm_gen -cb 4 3 1 -cb 2 26 40 -cb 73 295 147 -prec T -file mdtmdt_a1bc_md.c -routine_name mdtmdt_a1bc_md -beta c -opA T -opB T
**
** Run './../mm_gen/mm_gen -help' for help.
**
** Generated on: Thursday March 21 2002, 09:15:49 PST
** Created by: Jeff Bilmes <bilmes@cs.berkeley.edu>
**             http://www.icsi.berkeley.edu/~bilmes/phipac
**
**
** Usage:
**    mdtmdt_a1bc_md(const int M, const int K, const int N, const T *const A, const T *const B, T *const C, const int Astride, const int Bstride, const int Cstride, const T beta)
** where
**  transpose(A) is an MxK matrix
**  transpose(B) is an KxN matrix
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


#ifndef _BLPHIP_mdtmdt_a1bc_md_h_
#define _BLPHIP_mdtmdt_a1bc_md_h_ 1


BEGIN_BL_NAMESPACE



#define LOAD1x1(c00,C,Cstride) \
{\
   T * _cp = C; \
   c00 = _cp[0]; \
}

#define STORE1x1(c00,C,Cstride) \
{\
   T *_cp = C; \
   _cp[0] = c00; \
}

/* Fixed M,K,N = 1,1,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd1x1tmd1x1_md1x1(c00,A,Astride,B0) \
{ \
    T _a0; \
    T _b0; \
   \
   \
   _a0 = A[0]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; \
}


/* Fixed M,K,N = 1,2,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd1x2tmd2x1_md1x1(c00,A,Astride,B0) \
{ \
    T _a0; \
    T _b0; \
   \
   \
   _a0 = A[0]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; \
   \
   _a0 = A[0]; \
   A += Astride; \
   _b0 = B0[1]; \
   c00 += _b0*_a0; \
}


/* Fixed M,K,N = 1,3,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd1x3tmd3x1_md1x1(c00,A,Astride,B0) \
{ \
    T _a0; \
    T _b0; \
   \
   \
   _a0 = A[0]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; \
   \
   _a0 = A[0]; \
   A += Astride; \
   _b0 = B0[1]; \
   c00 += _b0*_a0; \
   \
   _a0 = A[0]; \
   A += Astride; \
   _b0 = B0[2]; \
   c00 += _b0*_a0; \
}


#define LOAD2x1(c00,c10,C,Cstride) \
{\
   T * _cp = C; \
   c00 = _cp[0]; \
   _cp += Cstride; \
   c10 = _cp[0]; \
}

#define STORE2x1(c00,c10,C,Cstride) \
{\
   T *_cp = C; \
   _cp[0] = c00; \
   _cp += Cstride; \
   _cp[0] = c10; \
}

/* Fixed M,K,N = 2,1,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd2x1tmd1x1_md2x1(c00,c10,A,Astride,B0) \
{ \
    T _a0,_a1; \
    T _b0; \
   \
   \
   _a0 = A[0]; _a1 = A[1]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; c10 += _b0*_a1; \
}


/* Fixed M,K,N = 2,2,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd2x2tmd2x1_md2x1(c00,c10,A,Astride,B0) \
{ \
    T _a0,_a1; \
    T _b0; \
   \
   \
   _a0 = A[0]; _a1 = A[1]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; c10 += _b0*_a1; \
   \
   _a0 = A[0]; _a1 = A[1]; \
   A += Astride; \
   _b0 = B0[1]; \
   c00 += _b0*_a0; c10 += _b0*_a1; \
}


/* Fixed M,K,N = 2,3,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd2x3tmd3x1_md2x1(c00,c10,A,Astride,B0) \
{ \
    T _a0,_a1; \
    T _b0; \
   \
   \
   _a0 = A[0]; _a1 = A[1]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; c10 += _b0*_a1; \
   \
   _a0 = A[0]; _a1 = A[1]; \
   A += Astride; \
   _b0 = B0[1]; \
   c00 += _b0*_a0; c10 += _b0*_a1; \
   \
   _a0 = A[0]; _a1 = A[1]; \
   A += Astride; \
   _b0 = B0[2]; \
   c00 += _b0*_a0; c10 += _b0*_a1; \
}


#define LOAD4x1(c00,c10,c20,c30,C,Cstride) \
{\
   T * _cp = C; \
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
   T *_cp = C; \
   _cp[0] = c00; \
   _cp += Cstride; \
   _cp[0] = c10; \
   _cp += Cstride; \
   _cp[0] = c20; \
   _cp += Cstride; \
   _cp[0] = c30; \
}

/* Fixed M,K,N = 4,1,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd4x1tmd1x1_md4x1(c00,c10,c20,c30,A,Astride,B0) \
{ \
    T _a0,_a1,_a2,_a3; \
    T _b0; \
   \
   \
   _a0 = A[0]; _a1 = A[1]; _a2 = A[2]; _a3 = A[3]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; c10 += _b0*_a1; c20 += _b0*_a2; c30 += _b0*_a3; \
}


/* Fixed M,K,N = 4,2,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd4x2tmd2x1_md4x1(c00,c10,c20,c30,A,Astride,B0) \
{ \
    T _a0,_a1,_a2,_a3; \
    T _b0; \
   \
   \
   _a0 = A[0]; _a1 = A[1]; _a2 = A[2]; _a3 = A[3]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; c10 += _b0*_a1; c20 += _b0*_a2; c30 += _b0*_a3; \
   \
   _a0 = A[0]; _a1 = A[1]; _a2 = A[2]; _a3 = A[3]; \
   A += Astride; \
   _b0 = B0[1]; \
   c00 += _b0*_a0; c10 += _b0*_a1; c20 += _b0*_a2; c30 += _b0*_a3; \
}


/* Fixed M,K,N = 4,3,1 fully-unrolled matrix matrix multiply. */
#define mul_tmd4x3tmd3x1_md4x1(c00,c10,c20,c30,A,Astride,B0) \
{ \
    T _a0,_a1,_a2,_a3; \
    T _b0; \
   \
   \
   _a0 = A[0]; _a1 = A[1]; _a2 = A[2]; _a3 = A[3]; \
   A += Astride; \
   _b0 = B0[0]; \
   c00 += _b0*_a0; c10 += _b0*_a1; c20 += _b0*_a2; c30 += _b0*_a3; \
   \
   _a0 = A[0]; _a1 = A[1]; _a2 = A[2]; _a3 = A[3]; \
   A += Astride; \
   _b0 = B0[1]; \
   c00 += _b0*_a0; c10 += _b0*_a1; c20 += _b0*_a2; c30 += _b0*_a3; \
   \
   _a0 = A[0]; _a1 = A[1]; _a2 = A[2]; _a3 = A[3]; \
   A += Astride; \
   _b0 = B0[2]; \
   c00 += _b0*_a0; c10 += _b0*_a1; c20 += _b0*_a2; c30 += _b0*_a3; \
}


/* Fixed M,N = 8,40, Arbitrary K L0-blocked matrix matrix multiply. */
template<class T>
static void
mdtmdt_a1bc_md_l1_arb_k(int K, const T *const A, const T *const B, T *const C, const int Astride, const int Bstride, const int Cstride)
{
   const T *a0,*b0;
   T *c0;
   const T *ap0;
   const T *bp0_0;
   T *cp0;
   const int B_sbs_stride = Bstride*1;
   const int C_sbs_stride = Cstride*4;
   const int k_marg_el = K % 3;
   const int k_norm = (K - k_marg_el)*Astride;
   T *const c0_endp = C+8*Cstride;
    T c00,c10,c20,c30;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=4) {
      const T* const ap0_endp = a0 + k_norm;
      T* const cp0_endp = c0 + 40;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=B_sbs_stride,cp0+=1) {
         ap0=a0;
         bp0_0 = b0;
         LOAD4x1(c00,c10,c20,c30,cp0,Cstride);
         for (; ap0!=ap0_endp; bp0_0+=3) {
            mul_tmd4x3tmd3x1_md4x1(c00,c10,c20,c30,ap0,Astride,bp0_0);
         }
         if (k_marg_el & 0x2) {
            mul_tmd4x2tmd2x1_md4x1(c00,c10,c20,c30,ap0,Astride,bp0_0);
            bp0_0+=2;
         }
         if (k_marg_el & 0x1) {
            mul_tmd4x1tmd1x1_md4x1(c00,c10,c20,c30,ap0,Astride,bp0_0);
         }
         STORE4x1(c00,c10,c20,c30,cp0,Cstride);
      }
   }
}

/* Arbitrary M,K,N L0-blocked matrix matrix multiply. */
template<class T>
static void
mdtmdt_a1bc_md_l1_arb_all(const int M, const int K, const int N, const T *const A, const T *const B, T *const C, const int Astride, const int Bstride, const int Cstride)
{
   const T *a0,*b0;
   T *c0;
   const T *ap0;
   const T *bp0_0;
   T *cp0;
   const int B_sbs_stride = Bstride*1;
   const int C_sbs_stride = Cstride*4;
   const int k_marg_el = K % 3;
   const int k_norm = (K - k_marg_el)*Astride;
   const int m_marg_el = M & 3;
   const int m_norm = M - m_marg_el;
   const int n_marg_el = N & 0;
   const int n_norm = N - n_marg_el;
   T *const c0_endp = C+m_norm*Cstride;
    T c00,c10,c20,c30;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=4) {
      const T* const ap0_endp = a0 + k_norm;
      T* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=B_sbs_stride,cp0+=1) {
         ap0=a0;
         bp0_0 = b0;
         LOAD4x1(c00,c10,c20,c30,cp0,Cstride);
         for (; ap0!=ap0_endp; bp0_0+=3) {
            mul_tmd4x3tmd3x1_md4x1(c00,c10,c20,c30,ap0,Astride,bp0_0);
         }
         if (k_marg_el & 0x2) {
            mul_tmd4x2tmd2x1_md4x1(c00,c10,c20,c30,ap0,Astride,bp0_0);
            bp0_0+=2;
         }
         if (k_marg_el & 0x1) {
            mul_tmd4x1tmd1x1_md4x1(c00,c10,c20,c30,ap0,Astride,bp0_0);
         }
         STORE4x1(c00,c10,c20,c30,cp0,Cstride);
      }
   }
   if (m_marg_el & 0x2) {
      const T* const ap0_endp = a0 + k_norm;
      T* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=B_sbs_stride,cp0+=1) {
         ap0=a0;
         bp0_0 = b0;
         LOAD2x1(c00,c10,cp0,Cstride);
         for (; ap0!=ap0_endp; bp0_0+=3) {
            mul_tmd2x3tmd3x1_md2x1(c00,c10,ap0,Astride,bp0_0);
         }
         if (k_marg_el & 0x2) {
            mul_tmd2x2tmd2x1_md2x1(c00,c10,ap0,Astride,bp0_0);
            bp0_0+=2;
         }
         if (k_marg_el & 0x1) {
            mul_tmd2x1tmd1x1_md2x1(c00,c10,ap0,Astride,bp0_0);
         }
         STORE2x1(c00,c10,cp0,Cstride);
      }
      c0 += Cstride*2;
      a0 += 2;
   }
   if (m_marg_el & 0x1) {
      const T* const ap0_endp = a0 + k_norm;
      T* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=B_sbs_stride,cp0+=1) {
         ap0=a0;
         bp0_0 = b0;
         LOAD1x1(c00,cp0,Cstride);
         for (; ap0!=ap0_endp; bp0_0+=3) {
            mul_tmd1x3tmd3x1_md1x1(c00,ap0,Astride,bp0_0);
         }
         if (k_marg_el & 0x2) {
            mul_tmd1x2tmd2x1_md1x1(c00,ap0,Astride,bp0_0);
            bp0_0+=2;
         }
         if (k_marg_el & 0x1) {
            mul_tmd1x1tmd1x1_md1x1(c00,ap0,Astride,bp0_0);
         }
         STORE1x1(c00,cp0,Cstride);
      }
   }
}

/* Fixed M,K,N = 8,78,40 L0-blocked matrix matrix multiply. */
template<class T>
static void
mdtmdt_a1bc_md_l1(const T *const A, const T *const B, T *const C, const int Astride, const int Bstride, const int Cstride)
{
   const T *a0,*b0;
   T *c0;
   const T *ap0;
   const T *bp0_0;
   T *cp0;
   const int B_sbs_stride = Bstride*1;
   const int C_sbs_stride = Cstride*4;
   const int k_norm = 78*Astride;
   T *const c0_endp = C+8*Cstride;
    T c00,c10,c20,c30;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=4) {
      const T* const ap0_endp = a0 + k_norm;
      T* const cp0_endp = c0 + 40;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=B_sbs_stride,cp0+=1) {
         ap0=a0;
         bp0_0 = b0;
         LOAD4x1(c00,c10,c20,c30,cp0,Cstride);
         for (; ap0!=ap0_endp; bp0_0+=3) {
            mul_tmd4x3tmd3x1_md4x1(c00,c10,c20,c30,ap0,Astride,bp0_0);
         }
         STORE4x1(c00,c10,c20,c30,cp0,Cstride);
      }
   }
}

/* Fixed M,N = 584,5880, Arbitrary K L1-blocked matrix matrix multiply. */
template<class T>
static void
mdtmdt_a1bc_md_l2_arb_k(int K, const T *const A, const T *const B, T *const C, const int Astride, const int Bstride, const int Cstride)
{
   /* Code for L1-blocked routine. */
   int m2,k2,n2;
   const T *a2,*b2;
   T *c2;
   const T *ap2,*bp2;
   T *cp2;
   for (m2=0; m2<584; m2+=8) {
      c2 = C + m2*Cstride;
      a2 = A + m2;
      for (n2=0,b2=B,cp2=c2; n2<5880; n2+=40,b2+=40*Bstride,cp2+=40) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-78; k2+=78,bp2+=78,ap2+=78*Astride) {
            mdtmdt_a1bc_md_l1(ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmdt_a1bc_md_l1_arb_k(K-k2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
   }
}

/* Arbitrary M,K,N L1-blocked matrix matrix multiply. */
template<class T>
static void
mdtmdt_a1bc_md_l2_arb_all(const int M, const int K, const int N, const T *const A, const T *const B, T *const C, const int Astride, const int Bstride, const int Cstride)
{
   /* Code for L1-blocked routine. */
   int m2,k2,n2;
   const T *a2,*b2;
   T *c2;
   const T *ap2,*bp2;
   T *cp2;
   if (M < 9 && K < 79 && N < 41) {
      mdtmdt_a1bc_md_l1_arb_all(M,K,N,A,B,C,Astride,Bstride,Cstride);
      return;
   }
   for (m2=0; m2<=M-8; m2+=8) {
      c2 = C + m2*Cstride;
      a2 = A + m2;
      for (n2=0,b2=B,cp2=c2; n2<=N-40; n2+=40,b2+=40*Bstride,cp2+=40) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-78; k2+=78,bp2+=78,ap2+=78*Astride) {
            mdtmdt_a1bc_md_l1(ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmdt_a1bc_md_l1_arb_k(K-k2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
      if (n2 < N) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-78; k2+=78,bp2+=78,ap2+=78*Astride) {
            mdtmdt_a1bc_md_l1_arb_all(8,78,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmdt_a1bc_md_l1_arb_all(8,K-k2,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
   }
   if (m2 < M) {
      c2 = C + m2*Cstride;
      a2 = A + m2;
      for (n2=0,b2=B,cp2=c2; n2<=N-40; n2+=40,b2+=40*Bstride,cp2+=40) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-78; k2+=78,bp2+=78,ap2+=78*Astride) {
            mdtmdt_a1bc_md_l1_arb_all(M-m2,78,40,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmdt_a1bc_md_l1_arb_all(M-m2,K-k2,40,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
      if (n2 < N) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-78; k2+=78,bp2+=78,ap2+=78*Astride) {
            mdtmdt_a1bc_md_l1_arb_all(M-m2,78,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mdtmdt_a1bc_md_l1_arb_all(M-m2,K-k2,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
   }
}

/* Fixed M,K,N = 584,23010,5880 L1-blocked matrix matrix multiply. */
template<class T>
static void
mdtmdt_a1bc_md_l2(const T *const A, const T *const B, T *const C, const int Astride, const int Bstride, const int Cstride)
{
   /* Code for L1-blocked routine. */
   int m2,k2,n2;
   const T *a2,*b2;
   T *c2;
   const T *ap2,*bp2;
   T *cp2;
   for (m2=0; m2<584; m2+=8) {
      c2 = C + m2*Cstride;
      a2 = A + m2;
      for (n2=0,b2=B,cp2=c2; n2<5880; n2+=40,b2+=40*Bstride,cp2+=40) {
         for (k2=0,bp2=b2,ap2=a2; k2<23010; k2+=78,bp2+=78,ap2+=78*Astride) {
            mdtmdt_a1bc_md_l1(ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
   }
}

template<class T>
void
mdtmdt_a1bc_md(const int M, const int K, const int N, const T *const A, const T *const B, T *const C, const int Astride, const int Bstride, const int Cstride, const T beta)
{
   /* Code for L2-blocked routine. */
   int m3,k3,n3;
   const T *a3,*b3;
   T *c3;
   const T *ap3,*bp3;
   T *cp3;
   {
      T *cprb,*cpre,*cp,*cpe;
      cpre = C + M*Cstride;
      for (cprb = C; cprb != cpre; cprb += Cstride) {
         cpe = cprb + N;
         for (cp = cprb; cp != cpe; cp++) {
            *cp *= (beta);
         }
      }
   }
   if (M < 585 && K < 23011 && N < 5881) {
      mdtmdt_a1bc_md_l2_arb_all(M,K,N,A,B,C,Astride,Bstride,Cstride);
      return;
   }
   for (m3=0; m3<=M-584; m3+=584) {
      c3 = C + m3*Cstride;
      a3 = A + m3;
      for (n3=0,b3=B,cp3=c3; n3<=N-5880; n3+=5880,b3+=5880*Bstride,cp3+=5880) {
         for (k3=0,bp3=b3,ap3=a3; k3<=K-23010; k3+=23010,bp3+=23010,ap3+=23010*Astride) {
            mdtmdt_a1bc_md_l2(ap3,bp3,cp3,Astride,Bstride,Cstride);
         }
         if (k3 < K) {
            mdtmdt_a1bc_md_l2_arb_k(K-k3,ap3,bp3,cp3,Astride,Bstride,Cstride);
         }
      }
      if (n3 < N) {
         for (k3=0,bp3=b3,ap3=a3; k3<=K-23010; k3+=23010,bp3+=23010,ap3+=23010*Astride) {
            mdtmdt_a1bc_md_l2_arb_all(584,23010,N-n3,ap3,bp3,cp3,Astride,Bstride,Cstride);
         }
         if (k3 < K) {
            mdtmdt_a1bc_md_l2_arb_all(584,K-k3,N-n3,ap3,bp3,cp3,Astride,Bstride,Cstride);
         }
      }
   }
   if (m3 < M) {
      c3 = C + m3*Cstride;
      a3 = A + m3;
      for (n3=0,b3=B,cp3=c3; n3<=N-5880; n3+=5880,b3+=5880*Bstride,cp3+=5880) {
         for (k3=0,bp3=b3,ap3=a3; k3<=K-23010; k3+=23010,bp3+=23010,ap3+=23010*Astride) {
            mdtmdt_a1bc_md_l2_arb_all(M-m3,23010,5880,ap3,bp3,cp3,Astride,Bstride,Cstride);
         }
         if (k3 < K) {
            mdtmdt_a1bc_md_l2_arb_all(M-m3,K-k3,5880,ap3,bp3,cp3,Astride,Bstride,Cstride);
         }
      }
      if (n3 < N) {
         for (k3=0,bp3=b3,ap3=a3; k3<=K-23010; k3+=23010,bp3+=23010,ap3+=23010*Astride) {
            mdtmdt_a1bc_md_l2_arb_all(M-m3,23010,N-n3,ap3,bp3,cp3,Astride,Bstride,Cstride);
         }
         if (k3 < K) {
            mdtmdt_a1bc_md_l2_arb_all(M-m3,K-k3,N-n3,ap3,bp3,cp3,Astride,Bstride,Cstride);
         }
      }
   }
}

END_BL_NAMESPACE



#endif


