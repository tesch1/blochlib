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

I'd like to thank to the makers of Gamma

  S.A. Smith, T.O. Levante, B.H. Meier and R.R. Ernst
  Computer Simulations in Magnetic Resonance:
  An Object-Oriented Programming Approach
  J. Magn. Reson., 106a, 75-105, (1994S, Series A)

http://gamma.magnet.fsu.edu/

for many of the algorithms used here
 **************************************************/


#ifndef _matdiagonalize_h_
#define _matdiagonalize_h_ 1

//#include "math.h"
#include "container/Vector/Vector.h"
#include "container/rankType.h"
#include "container/matrix/_matrix.h"
#include <iostream>
BEGIN_BL_NAMESPACE


//*************************************************************************
// reduces a REAL Symmetric	matrix to a Re	**
// tridiagonal form using the Housholder algorithm.		**

template<class nt1>
void tred( _matrix<nt1, FullMatrix> &in,
		 _matrix<nt1, DiagonalMatrix> &d,
		  Vector<nt1> &e)
{
	const int n=in.rows();
	FloatType(nt1) scale, h, f,g;
	int i,j,k;
	for (i=n-1;i>=1;i--) {
		const int l=i-1;
		h=scale=ZeroType<FloatType(nt1)>::zero();

		if (l>0) {
			for (k=l;k--;) scale+=abs(in(i,k));

			if (scale==ZeroType<FloatType(nt1)>::zero()){
				e(i)=in(i,l);
			}else {
				scale+=abs(in(i,l));
				for (k=0;k<i;k++) {
					in(i,k) /= scale;
					h+=in(i,k)*in(i,k);
				}
				f=in(i,l);
				g=-sign(std::sqrt(h),f);
				e(i)=scale*g;
				h-=f*g;
				in(i,l)=f-g;
				f=ZeroType<FloatType(nt1)>::zero();
				for (j=0;j<i;j++) {
					in(j,i)=in(i,j)/h;
					g=ZeroType<FloatType(nt1)>::zero();
					for (k=0;k<=j;k++)	g+=in(j,k)*in(i,k);
					if (l>j) {
						for (k=j+1;k<i;k++)	g+=in(k,j)*in(i,k);
					}
					e[j]=g/h;
					f+=e(j)*in(i,j);
				}
				const nt1 hh=f/(h+h);
				for (j=0;j<i;j++) {
					f=in(i,j);
					g=e(j)-hh*f;
					e[j]=g;
					for (k=0;k<=j;k++)	in(j,k)-=f*e(k)+g*in(i,k);
				}
			}
		}else{
			e[i]=in(i,l);
		}
		d(i,i)=h;
	}
	d(0,0)=ZeroType<nt1>::zero();
	e(0)=ZeroType<nt1>::zero();

	for (i=0;i<n;i++) {
		if (d(i,i)!=0.0) {
			for (j=0;j<i;j++) {
				g=ZeroType<FloatType(nt1)>::zero();
				for (k=0;k<i;k++)  g += in(i,k)*in(k,j);
				for (k=0;k<i;k++)  in(k,j) -= g*in(k,i);
			}
		}
		d(i,i)=in(i,i);
		in(i,i)=OneType<nt1>::one();
		for (j=0;j<i;j++)  in(i,j)=in(j,i)=ZeroType<nt1>::zero();
	}
}

template<class nt1, class nt2>
OutType(nt1, nt2) pyth(nt1 a, nt2 b){
	OutType(nt1, nt2) aba=abs(a);
	OutType(nt1, nt2) abb=abs(b);
	if(aba > abb){
		OutType(nt1, nt2) adb=abb/aba;
		return aba*std::sqrt(1.+adb*adb);
	}else{
		OutType(nt1, nt2) adb=aba/abb;
		return (aba==ZeroType<OutType(nt1, nt2)>::zero() ? ZeroType<OutType(nt1, nt2)>::zero() : aba*std::sqrt(1+adb*adb));
	}
}


template<class nt1, class nt2>
inline OutType(nt1, nt2) sign(nt1 a, nt2 b) {
	return ((b) >= ZeroType<nt2>::zero() ? abs(a) : -abs(a));
}

template<class nt1>
void realtqli(	 _matrix<nt1, DiagonalMatrix> &d,
				  Vector<nt1> &e,
				 _matrix<nt1, FullMatrix> &U)
{
	int i,k,l,m;
	int n=d.rows();
	nt1 g,r,s,c,p,f,b,DD;
	for (i=1;i<n;i++) e(i-1)=e(i);
	e(n-1)=ZeroType<nt1>::zero();


	for (l=0;l<n;l++) {
		for (int iter=0;;iter++) {
			for (m=l;m<n-1;m++) {
				DD=abs(d(m,m))+abs(d(m+1,m+1));
				if ((abs(e(m))+DD)==DD) break;
			}
			if (m==l) break;

			if (iter>10*n) {
				BLEXCEPTION(" too many itterations perhaps strange input")
			}
			g=(d(l+1,l+1)-d(l,l))/(2*e(l));
			r=std::sqrt(g*g+1.0);
			g=d(m,m)-d(l,l)+e(l)/(g+sign(r,g));
			s=1.0;
			c=1.0;
			p=0.0;

			for (i=m-1;i>=l;i--) {
				f=s*e(i);
				b=c*e(i);
				if (abs(f)>=abs(g)) {
					c=g/f;
					r=std::sqrt(c*c+1.0);
					e(i+1)=f*r;
					s=1.0/r;
					c*=s;
				}else {
					s=f/g;
					r=std::sqrt(s*s+1.0);
					e(i+1)=g*r;
					c=1.0/r;
					s*=c;
				}
				g=d(i+1, i+1)-p;
				r=(d(i,i)-g)*s+2*c*b;
				p=s*r;
				d(i+1,i+1)=g+p;
				g=c*r-b;
				for (k=0;k<n;k++) {
					nt1 x=U(k, i+1);
					nt1 y=U(k,i);
					U(k,i+1)=s*y+c*x;
					U(k,i)=c*y-s*x;
				}
			}
			d(l,l)-=p;
			e(l)=g;
			e(m)=ZeroType<nt1>::zero();
		}
	}
}

template<class nt1, class st1>
void RealSymmetricDiag(const	 _matrix<nt1,st1> &inmat,
						  _matrix<nt1, DiagonalMatrix> &diag,
						  _matrix<nt1, FullMatrix> &U,
						  bool destory=false)
{
	int rs=inmat.rows(0);
	if(U.rows(0) != rs || U.cols(0) != rs) U.resize(rs,rs);
	if(diag.rows(0) != rs || diag.cols(0) != rs) diag.resize(rs,rs);
	diag.fill(0);
	if(inmat.type()==Mdiagonal  || inmat.type()==Midentity){
		diag=inmat;
		U.fill(ZeroType<nt1>::zero());
		for(int i=0;i<rs;++i) U(i,i)=OneType<nt1>::one();
		return;
	}

	if(inmat.type() != Msymmetric && inmat.type() != Mhermitian){
		BLEXCEPTION(" matrix to diagonalize MUST be a symmetric type")
	}

	Vector<nt1> e(rs, ZeroType<nt1>::zero());
	U=inmat;
	tred(U, diag, e);
	realtqli(diag, e, U);
	/*std::cout<<"U::"<<std::endl<<chop(U)<<std::endl
			<<"prop(U,diag)::"<<std::endl<<chop(prop(U,diag))<<std::endl
			<<"U*Udagar::"<<std::endl<<chop(U*transpose(U))<<std::endl
		<<"evals::"<<chop(diag)<<std::endl;*/
}

///*********************END REAL SYMMETRIC MATRIX DIAG************************

///*********************BEGIN COMPLEX HERMETIAN MATRIX DIAG*******************

//*************************************************************************
// reduces a complex Hermitian	matrix to a complex	**
// Hermitian tridiagonal form using the Housholder algorithm.		**


template<class nt1, class st1>
void cred(	 _matrix<nt1, st1> &in,
			  _matrix<nt1, FullMatrix> &U)
{
   int i,j,l;
   int n=in.rows();
   double s;			// norm of the vector a(l+1,l) to a(n-1,l)
   double sw;			// norm of the vextor a(l+2,l) to a(n-1,l)
   double t;

   Vector<nt1> w(n);

   /* create unit matrix */
   U.identity();
	/* Housholder reduction */
   for (l=0;l<n-2;l++) {

	double u=ZeroType<double>::zero();			// column scaling factor
	for (i=l+2;i<n;i++)
		u+=AbsNorm(in(i,l));

	if (u>ZeroType<double>::zero()) {
		const nt1 all=in(l+1,l);
		u+=AbsNorm(all);

		sw=ZeroType<double>::zero();
		for (i=l+2;i<n;i++)	sw+=square_norm(in(i,l)/u);

		s=sw+square_norm(all/u);
		s=u*std::sqrt(s);
		if (all!=ZeroType<double>::zero()) {
			t=norm(all);
			w(l+1)=all*(1.0+s/t);
		}else{
			w(l+1)= -s;
		}
		s=sw+square_norm(w(l+1)/u);
		s=std::sqrt(2/s)/u;

		for (i=l+2;i<n;i++) 	w(i)=in(i,l)*s;

		w(l+1)*=s;

		nt1 ss;
		/* transform matrix from the left */
		for (i=l;i<n;i++) {
			ss=0.0;
			for (j=l+1;j<n;j++) muladd_conj(ss,in(j,i),w(j));
			ss=-ss;

			for (j=l+1;j<n;j++) muladd(in(j,i),ss,w(j));
		}
		/* transform matrix from the right */
		for (i=l;i<n;i++) {
			ss=0.0;
			for (j=l+1;j<n;j++) muladd(ss,w(j),in(i,j));
			ss=-ss;
			for (j=l+1;j<n;j++) muladd_conj(in(i,j),ss,w(j));
       }

       /* accumulate transformations */
		for (i=0;i<n;i++) {
			ss=0.0;
			for (j=l+1;j<n;j++) muladd(ss,w(j),U(i,j));
			ss=-ss;
			for (j=l+1;j<n;j++) muladd_conj(U(i,j),ss,w(j));
		}
     }
  }
}





template<class nt1, class st1, class st2>
void rred(	 _matrix<nt1, st1> &in,
			  _matrix<nt1, st2> &ev)
{
	int i,l;
	int n=in.rows();
	double s;
	nt1 ss;

	for (l=1;l<n;l++) {
		s=norm(in(l,l-1));
		if (s!=0) {		// if there is something to be done
			ss=conj(in(l,l-1))/s;
			in(l,l-1)=in(l-1,l)=s;
			if (l<n-1) {		// transform
				in(l,l+1)*=ss;
				in(l+1,l)*=conj(ss);
			}
			for (i=0;i<n;i++)	ev(i,l)*=conj(ss);
		}
	}
}

template<class nt1, class st1>
void complextqli(	 _matrix<nt1, st1> &in,
					  _matrix<nt1, DiagonalMatrix> &d,
					  _matrix<nt1, FullMatrix> &U)
{

	int i,k,l,m,iter;
	int n=in.rows();
	float dd;
	float s,r,p,g,f,c,b;
	float x;
	nt1 xx,yy;

	Vector<float> e(n);

	for (i=0;i<n-1;i++)	e(i)=Re(in(i,i+1));

	for (i=n;i--;) d(i,i)=Re(in(i,i));

	e(n-1)=0.0;
	iter=0;
	for (l=0;l<n;l++) {		// loop for eigenvalues
		do {
			m=l;
			while (1) {		// look for an eigenvalue
				if (m==n-1) break;	// all eigenvalues found
				dd=abs(d(m,m))+abs(d(m+1, m+1));
				if (abs(e(m))+dd==dd) break;	// a new e.value found
				m+=1;
			}
			if (m!=l) {
				if (iter>10*n) {
					BLEXCEPTION(" too many itterations perhaps strange input")
				}
				// this happens almost surely only on erroneous input

				iter++;
				// form shift
				g=(Re(d(l+1, l+1)-d(l,l)))/(2.0*e(l));
				r=std::sqrt(g*g+1.0);
				x=e(l)/(g+sign(r,g));
				g=Re(d(m,m)-d(l,l))+x;
				s=1.0;
				c=1.0;
				p=0.0;
				// QL step
				for (i=m-1;i>=l;i--) {
					f=s*e(i);
					b=c*e(i);
					if (abs(f)>=abs(g)) {
						c=g/f;
						r=std::sqrt(c*c+1.0);
						e(i+1)=f*r;
						s=1.0/r;
						c*=s;
					}else {
						s=f/g;
						r=std::sqrt(s*s+1.0);
						e(i+1)=g*r;
						c=1.0/r;
						s*=c;
					}
					g=Re(d(i+1, i+1))-p;
					r=(Re(d(i,i))-g)*s+2.0*c*b;
					p=s*r;
					Re(d(i+1, i+1),g+p);
					g=c*r-b;
					// accumulate transforms
					for (k=n;k--;) {
						nt1& zi =U(k,i);
						nt1& zi1=U(k,i+1);
						xx=zi; xx*=s;
						yy=zi1; yy*=s;
						zi*=c; zi-=yy;
						zi1*=c; zi1+=xx;
					}
				}
				d(l,l)-=p;
				e(l)=g;
				e(m)=0.0;
			}
		} while (m!=l);		// all eigenvalues found
	}
}

template<class nt1, class st1>
void ComplexHermitianDiag(const _matrix<nt1,st1> &inmat,
						 _matrix<nt1, DiagonalMatrix> &diag,
						 _matrix<nt1, FullMatrix> &U,
						const bool destroy=false)
{
	int rs=inmat.rows(0);
	if(inmat.type()==Mdiagonal  || inmat.type()==Midentity){
		diag=inmat;
		U.identity();
		return;
	}

	if(inmat.type() != Mhermitian ){
		BLEXCEPTION(" matrix to diagonalize MUST be a hermitian type for complex numbers")
	}
	if(U.rows(0) != rs || U.cols(0) != rs) U.resize(rs,rs);
	if(diag.rows(0) != rs || diag.cols(0) != rs) diag.resize(rs,rs);
	diag.fill(0);

	U=inmat;
	if(!destroy){
		_matrix<nt1, FullMatrix> tmm=inmat;
		try{
			cred(tmm, U); 	
			rred(tmm, U);	
			complextqli(tmm,diag, U);
		}catch(BL_exception e){
			e.print(std::cerr);
			std::cerr<<"Original Matrix: "<<inmat;
			BLEXCEPTION("leaving diagonalize");
		}
	}else{
		_matrix<nt1, st1> &hold=const_cast<_matrix<nt1, st1> &>(inmat);
		try{
			cred(hold, U);
			rred(hold, U);
			complextqli(hold,diag, U);		
		}catch(BL_exception e){
			e.print(std::cerr); 
			std::cerr<<"Original Matrix: "<<inmat;
			BLEXCEPTION("leaving diagonalize");
		}
	}
/*	std::cout<<"U::"<<std::endl<<chop(U)<<std::endl
				<<"prop(U,diag)::"<<std::endl<<chop(prop(U,diag))<<std::endl
				<<"U*Udagar::"<<std::endl<<chop(U*adjoint(U))<<std::endl
		<<"evals::"<<chop(diag)<<std::endl;*/

}

//******************END HERMITIAN COMPLEX DIAG******************

//******************BEGIN GENERAL MATRIX DIAG*******************
/*
     this subroutine is a translation of a complex analogue of
     the algol procedure orthes, num. math. 12, 349-368(1968)
     by martin and wilkinson.
     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

     given a complex general matrix, this subroutine
     reduces a submatrix situated in rows and columns
     low through igh to upper hessenberg form by
     unitary similarity transformations.

     on input
        low and igh are integers determined by the balancing
          subroutine  cbal.  if  cbal  has not been used,
          set low=1, igh=n.

        a contains the complex input matrix.

     on output

        a of the hessenberg matrix.  information
          about the unitary transformations used in the reduction
          is stored in the remaining triangles under the
          hessenberg matrix.

        ort contains further information about the
          transformations.  only elements low through igh are used.

*/
//For a COMPLEX input 'A' matrix
//NOTE:: ORT MUST be complex
template<class st1, class Ctype>
void corth(	 _matrix<Complex<Ctype>, st1> &a,
			 Vector<Complex<Ctype> > &ort,
			int low,
			int igh)
{
  int i,j,m;
  double f,g,h,sc;
  Complex<Ctype> z;
  const int n=a.rows();

	for (m=low+1; m<igh; m++){
		h=0;
		sc=0;
		for (i=m; i<=igh; i++) sc+=AbsNorm( a(i,m-1) );

		if (sc>0.0) {
			for (i=igh;i>=m;i--) {
				ort(i)= a(i,m-1)/sc;
				h+=square_norm(ort(i));
			}
			g=std::sqrt(h);
			f=norm(ort(m));
			if (f==0) {
				ort(m)=g;
				Re(a(m,m-1),sc);
			}else{
				h+=f*g;
				g/=f;
				ort(m) *= (1.0+g);
			}

			for (j=m;j<n;j++) {
				z=0.0;
				for (i=igh;i>=m;--i)	muladd_conj(z,a(i,j),ort(i));
				z/=h;
				for (i=m;i<=igh;i++)	a(i,j)-=z*ort(i);
			}

			for (i=0;i<=igh;i++) {
				z=0.0;
				for (j=igh;j>=m;--j) muladd(z,ort(j),a(i,j));
				z/=h;
				for (j=m;j<=igh;j++) a(i,j) -= conj(ort(j))*z;
			}

			ort(m) *= sc;
			a(m,m-1) *= -g;
		}
	}
}

//For a REAL input 'A' matrix
//NOTE:: ORT Re also be complex the complex version is above
template<class nt1,class st1>
void corth(	 _matrix<nt1, st1> &a,
			 Vector<nt1> &ort,
			int low,
			int igh)
{
  int i,j,m;
  double f,g,h,sc;
  double z;
  const int n=a.rows();

	for (m=low+1; m<igh; m++){
		h=0;
		sc=0;
		for (i=m; i<=igh; i++) sc+=abs( a(i,m-1) );

		if (sc>0.0) {
			for (i=igh;i>=m;i--) {
				ort(i)= a(i,m-1)/sc;
				h+=ort(i)*ort(i);
			}
			g=std::sqrt(h);
			f=abs(ort(m));
			if (f==0) {
				ort(m)=g;
				a(m,m-1)=sc;
			}else{
				h+=f*g;
				g/=f;
				ort(m) *= (1.0+g);
			}

			for (j=m;j<n;j++) {
				z=0.0;
				for (i=igh;i>=m;--i)	z+=a(i,j)*ort(i);
				z/=h;
				for (i=m;i<=igh;i++)	a(i,j)-=z*ort(i);
			}

			for (i=0;i<=igh;i++) {
				z=0.0;
				for (j=igh;j>=m;--j) z+=ort(j)*a(i,j);
				z/=h;
				for (j=m;j<=igh;j++) a(i,j) -= ort(j)*z;
			}

			ort(m) *= sc;
			a(m,m-1) *= -g;
		}
	}
}

/*
     this subroutine is a translation of a unitary analogue of the
     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
     and wilkinson.
     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
     the unitary analogue substitutes the qr algorithm of francis
     (comp. jour. 4, 332-345(1962)) for the lr algorithm.

     this subroutine finds the eigenvalues and eigenvectors
     of a complex upper hessenberg matrix by the qr
     method.  the eigenvectors of a complex general matrix
     can also be found if  corth  has been used to reduce
     this general matrix to hessenberg form.

     on input

        low and igh are integers determined by the balancing
          subroutine  cbal.  if  cbal  has not been used,
          set low=1, igh=n.

        ort contains information about the unitary trans-
          formations used in the reduction by  corth, if performed.
          only elements low through igh are used.  if the eigenvectors
          of the hessenberg matrix are desired, set ortr(j) and
          orti(j) to 0.0d0 for these elements.

        h contains the complex upper hessenberg matrix.
          their lower triangles below the subdiagonal contain further
          information about the transformations which were used in the
          reduction by  corth, if performed.  if the eigenvectors of
          the hessenberg matrix are desired, these elements may be
          arbitrary.

     on output

        ort, and the upper hessenberg portions of h
          have been destroyed.

        w contains the eigenvalues.  if an error
          exit is made, the eigenvalues should be correct
          for indices ierr+1,...,n.

        z comtains the eigenvectors.  the eigenvectors
          are unnormalized.  if an error exit is made, none of
          the eigenvectors has been found.

        ierr is set to
          zero       for normal return,
          j          if the limit of 30*n iterations is exhausted
                     while the j-th eigenvalue is being sought.

*/

//this is the routine for a COMPLEX input matrix "a"
template<class st1, class Ctype>
 int comqr3( _matrix<Complex<Ctype>, st1> &a,int low,int igh,
		      Vector<Complex<Ctype> > &ort,
		      _matrix<Complex<Ctype>, DiagonalMatrix> &w,
		      _matrix<Complex<Ctype>, FullMatrix> & z,
		     int flag){
 	int n=a.rows();
	int i,j,k,l=0,m,en,ii,jj,ll,nn,its,iend,enm1,lp1, ip1;
	double u,s,ss,nor;
	Complex<Ctype> sc,x,y,t,zz;

	const double machep=1.0e-12;

//do if eignvectors are desired
	if(flag)	z.identity();

	iend=igh-low-1;
	if(iend>=0){
		for(ii=1;ii<=iend;ii++){
			i=igh-ii;
			if (ort[i-1]!=0.0){
				if (a(i-1,i-2)!=0.0){
					nor=Re(conj(a(i-1,i-2),ort[i-1]));
					ip1=i+1;
					for (k=ip1;k<=igh;k++)	ort[k-1]=a(k-1,i-2);
					if(flag) //do for evects only
					{
						for (j=i;j<=igh;j++){
							sc=0.0;
							for (k=i;k<=igh;k++)
							sc += conj(ort[k-1],z(k-1,j-1));
							sc/=nor;
						//do if eignvectors are desired
							for (k=i;k<=igh;k++)	z(k-1,j-1) +=sc*ort[k-1];
						}
					}
				}
			}
		}
      // Re subdiagonal elements
		l=low+1;
		for(i=l;i<=igh;i++){
			ll= ((i+1)<igh)?i+1:igh;
			if(Im(a(i-1,i-2))!=0.0){
				nor=norm(a(i-1,i-2));
				y=a(i-1,i-2)/nor;
				a(i-1,i-2)=nor;
				for (j=i;j<=n;j++)		a(i-1,j-1)=conj(y,a(i-1,j-1));
				for (j=1;j<=ll;j++)		a(j-1,i-1)*=y;
			//do if eignvectors are desired
				if(flag)	for (j=low;j<=igh;j++)	z(j-1,i-1)*=y;
			}
		}
	}

  // isolated roots
	for (i=1;i<=n;i++){
		if ((i<low)||(i>igh)){
			w(i-1,i-1)=a(i-1,i-1);
		}
	}


  en=igh;
  t=0.0;
//cout<<"A COMPLEX: "<<a<<endl;

  //iterate for eigenvalues
	while (en>=low){
		its=0;
		enm1=en-1;
		while (1){
			for (ll=low;ll<=en;ll++){
				l=en+low-ll;
				if (l==low) break;
				s=AbsNorm(a(l-2,l-2))+AbsNorm(a(l-1,l-1));
				if (s<1.0) s=1.0;
				if (fabs(Re(a(l-1,l-2)))<=machep*s) break;
			}
			if (l==en) break;
			if (its==30){
				BLEXCEPTION(" too many itterations perhaps strange input")
			}
		// form shift
			if ((its==10)||(its==20)){
				sc=fabs(Re(a(en-1,enm1-1)))+ fabs(Re(a(enm1-1,en-3)));
			}else{
				sc=a(en-1,en-1);
				x=a(enm1-1,en-1)*Re(a(en-1,enm1-1));
				if (x!=0.0){
					y=(a(enm1-1,enm1-1)-sc)/2.0;
					zz = std::sqrt(y*y+x);
					if ((Re(y)*Re(zz)+Im(y)*Im(zz))<0.0)	zz = -zz;
					zz = x/(y+zz);
					sc-=zz;
				}
			}
			for (i=low;i<=en;i++) a(i-1,i-1) -= sc;
			t+=sc;
			its=its+1;
			// QR decomposition
			lp1=l+1;
			for (i=lp1;i<=en;i++){
				s=Re(a(i-1,i-2));
				nor=std::sqrt( square_norm(a(i-2,i-2))+s*s );
				x=a(i-2,i-2)/nor;
				w(i-2,i-2)=x;
				a(i-2,i-2)=nor;
				s /= nor;
				Re(a(i-1,i-2),0.);
				Im(a(i-1,i-2),s);// a=complex(0,s);
				for (j=i;j<=n;j++){
					y=a(i-2,j-1);
					zz=a(i-1,j-1);
					a(i-2,j-1)= conj(x,y)+s*zz;
					a(i-1,j-1)= x*zz - s*y;
				}
			}
			const bool isreal=(Im(a(en-1,en-1))==0.);
			if(!isreal){
				nor=norm(a(en-1,en-1));
				sc=conj(a(en-1,en-1)/nor);
				a(en-1,en-1)=nor;
				if (en!=n){
					for (j=en+1;j<=n;j++)	a(en-1,j-1)*=sc;
				}
			}
		  // calculate RQ
			for (j=lp1;j<=en;j++){
				x=w(j-2,j-2);
				const complex cx=conj(x);
				ss=Im(a(j-1,j-2));
				for (i=1;i<=j;i++){
				y=Re(a(i-1,j-2));
				zz=a((i-1),j-1);
				if (i!=j){
					Im(y,Im(a((i-1),j-2)));
					Im(a((i-1),j-2), Re(x)*Im(y)+Im(x)*Re(y)+ss*Im(zz));
				}
				Re(a((i-1),j-2),Re(x)*Re(y)-Im(x)*Im(y)+ss*Re(zz));
				a((i-1),j-1)=cx*zz-ss*y;
			}
			if(flag) //do for evects only
			{
				for (i=low;i<=igh;i++){
					y=z((i-1),j-2);
					zz=z((i-1),j-1);
					z((i-1),j-2)=x*y+ss*zz;
					z((i-1),j-1)=cx*zz-ss*y;
				};
			}
		};
		sc=conj(sc);
		if (!isreal){
			for (i=1;i<=en;i++)		a((i-1),en-1)*=sc;
			for (i=low;i<=igh;i++)	z((i-1),en-1)*=sc;
		};
	};
     // a root found
		a((en-1),en-1)+=t;
		w(en-1,en-1)=a((en-1),en-1);
		en=enm1;

    };

//do if eignvectors are desired
	if(flag){
      // all roots found, backsubstitution of eigenvectors
		nor=0.0;
		for (i=1;i<=n;i++){
			for (j=i;j<=n;j++){
				nor+=AbsNorm(a((i-1),j-1));
			}
		}
		if ((n==1)||(nor==0.0)) return 0;
		for (nn=2;nn<=n;nn++){
			en=n+2-nn;
			x=w(en-1,en-1);
			a(en-1,en-1)=1;
			enm1=en-1;
			for(ii=1;ii<=enm1;ii++){
				i=en-ii;
				zz=a((i-1),en-1);
				if (i!=enm1){
					ip1=i+1;
					for (j=ip1;j<=enm1;j++)	muladd(zz,a((i-1),j-1),a((j-1),en-1)); //zz+= a((i-1),j-1)*a((j-1),en-1);
				}
				y=x-w(i-1,i-1);
				if ((fabs(Re(y))<machep)&&(fabs(Im(y))<machep))	Re(y,machep*nor);
				a((i-1),en-1)=zz/y;
			}
		}

		enm1=n-1;
	  // eigenvectors of isolated root
		for (i=1;i<=enm1;i++){
			if ((i<low)||(i>igh)){
				for (j=i+1;j<=n;j++)	z((i-1),j-1)=a((i-1),j-1);
			}
		}
	  // multiply for eigenbase
		for (jj=low;jj<=enm1;jj++){
			j=n+low-jj;
			m= ((j-1)<igh)?j-1:igh;
			for (i=low;i<=igh;i++){
				zz=z((i-1),j-1);
				for (k=low;k<=m;k++)	muladd(zz,z((i-1),k-1),a((k-1),j-1)); //zz+=z((i-1),k-1)*a((k-1),j-1);
				z((i-1),j-1)=zz;
			}
		}
	  //normalize vectors
		for (i=0;i<n;i++){
			u=0.0;
			for (j=0;j<n;j++)	u+=AbsNorm(z(j,i));
			s=0.0;
			for (j=0;j<n;j++)	s+=square_norm(z(j,i)/u);
			s=u*std::sqrt(s);
			for (j=0;j<n;j++)	z(j,i)/=s;
		}
	}
	return 0;
}


//this is the routine for a REAL input matrix "a"
// Ort is also Re
//the complex version is above...
// NOTE:: the matrix 'z' is still needed as a workspace
//  if NOT obtaining eignvectors
template<class nt1, class st1, class Ctype>
 int comqr3( _matrix<nt1, st1> &a,int low,int igh,
		      Vector<nt1> &ort,
		      _matrix<Complex<Ctype>, DiagonalMatrix> &w,
		      _matrix<Complex<Ctype>, FullMatrix> & z,
		     int flag){
 	int n=a.rows();
	int i,j,k,l=0,m,en,ii,jj,ll,nn,its,iend,enm1,lp1, ip1;
	double u,s,ss,nor,y;
	Complex<Ctype> sc,x,t,zz,yc;

	const double machep=1.0e-12;

//init the workspace for z
	if(flag) z.identity();

	iend=igh-low-1;
	if(iend>=0){
		for(ii=1;ii<=iend;ii++){
			i=igh-ii;
			if (ort[i-1]!=0.0){
				if (a(i-1,i-2)!=0.0){
					nor=a(i-1,i-2)*ort[i-1];
					ip1=i+1;
					for (k=ip1;k<=igh;k++)	ort[k-1]=a(k-1,i-2);

					if(flag) //do only if eigvects are desired
					{
						for (j=i;j<=igh;j++){
							sc=0.0;
							for (k=i;k<=igh;k++)
							sc += ort[k-1]*conj(z(k-1,j-1));
							sc/=nor;
							for (k=i;k<=igh;k++)	z(k-1,j-1) +=sc*ort[k-1];
						}
					}
				}
			}
		}
      // Re subdiagonal elements
	/*	l=low+1;
		for(i=l;i<=igh;i++){
			ll= ((i+1)<igh)?i+1:igh;
			if(Im(a(i-1,i-2))!=0.0){
				nor=abs(a(i-1,i-2));
				y=a(i-1,i-2)/nor;
				a(i-1,i-2)=nor;
				for (j=i;j<=n;j++)		a(i-1, j-1)*=y;//a(i-1,j-1)=conj(y,a(i-1,j-1));
				for (j=1;j<=ll;j++)		a(j-1,i-1)*=y;
			//do if eignvectors are desired
				if(flag)	for (j=low;j<=igh;j++)	z(j-1,i-1)*=y;
			}
		}
	*/
	}

  // isolated roots
	for (i=1;i<=n;i++){
		if ((i<low)||(i>igh)){
			w(i-1,i-1)=a(i-1,i-1);
		}
	}


  en=igh;
  t=0.0;

//do far everything up to here can be done using only REAL algebra...
//and everything was done internally on 'a' and 'ort'
//now, however, things must enter complex algebra as the roots of out
// 'polynomial' can be complex....
// becuase 'a' is Re and not able to hold complex numbers, we need
// to use a temporary workspace...i am not sure how much CPU time
// the 'Re' comqr3 saves v the 'complex comqr3 dues to this conversion

	_matrix<Complex<Ctype>, FullMatrix> A=a;
	//cout<<"A REAL: "<<A<<endl;
   //iterate for eigenvalues
  //iterate for eigenvalues
 	while (en>=low){
 		its=0;
 		enm1=en-1;
 		while (1){
 			for (ll=low;ll<=en;ll++){
 				l=en+low-ll;
 				if (l==low) break;
 				s=AbsNorm(A(l-2,l-2))+AbsNorm(A(l-1,l-1));
 				if (s<1.0) s=1.0;
 				if (fabs(Re(A(l-1,l-2)))<=machep*s) break;
 			}
 			if (l==en) break;
 			if (its==30){
				BLEXCEPTION(" too many itterations perhaps strange input")
 			}
 		// form shift
 			if ((its==10)||(its==20)){
 				sc=fabs(Re(A(en-1,enm1-1)))+ fabs(Re(A(enm1-1,en-3)));
 			}else{
 				sc=A(en-1,en-1);
 				x=A(enm1-1,en-1)*Re(A(en-1,enm1-1));
 				if (x!=0.0){
 					yc=(A(enm1-1,enm1-1)-sc)/2.0;
 					zz = std::sqrt(yc*yc+x);
 					if ((Re(yc)*Re(zz)+Im(yc)*Im(zz))<0.0)	zz = -zz;
 					zz = x/(yc+zz);
 					sc-=zz;
 				}
 			}
 			for (i=low;i<=en;i++) A(i-1,i-1) -= sc;
 			t+=sc;
 			its=its+1;
 			// QR decomposition
 			lp1=l+1;
 			for (i=lp1;i<=en;i++){
 				s=Re(A(i-1,i-2));
 				nor=std::sqrt( square_norm(A(i-2,i-2))+s*s );
 				x=A(i-2,i-2)/nor;
 				w(i-2,i-2)=x;
 				A(i-2,i-2)=nor;
 				s /= nor;
 				Re(A(i-1,i-2),0.);
 				Im(A(i-1,i-2),s);// a=complex(0,s);
 				for (j=i;j<=n;j++){
 					yc=A(i-2,j-1);
 					zz=A(i-1,j-1);
 					A(i-2,j-1)= conj(x,yc)+s*zz;
 					A(i-1,j-1)= x*zz - s*yc;
 				}
 			}
 			const bool isreal=(Im(A(en-1,en-1))==0.);
 			if(!isreal){
 				nor=norm(A(en-1,en-1));
 				sc=conj(A(en-1,en-1)/nor);
 				A(en-1,en-1)=nor;
 				if (en!=n){
 					for (j=en+1;j<=n;j++)	A(en-1,j-1)*=sc;
 				}
 			}
 		  // calculate RQ
 			for (j=lp1;j<=en;j++){
 				x=w(j-2,j-2);
 				const complex cx=conj(x);
 				ss=Im(A(j-1,j-2));
 				for (i=1;i<=j;i++){
 				yc=Re(A(i-1,j-2));
 				zz=A((i-1),j-1);
 				if (i!=j){
 					Im(yc,Im(A((i-1),j-2)));
 					Im(A((i-1),j-2), Re(x)*Im(yc)+Im(x)*Re(yc)+ss*Im(zz));
 				}
 				Re(A((i-1),j-2),Re(x)*Re(yc)-Im(x)*Im(yc)+ss*Re(zz));
 				A((i-1),j-1)=cx*zz-ss*yc;
 			}
 			if(flag) //do for evects onlyc
 			{
 				for (i=low;i<=igh;i++){
 					yc=z((i-1),j-2);
 					zz=z((i-1),j-1);
 					z((i-1),j-2)=x*yc+ss*zz;
 					z((i-1),j-1)=cx*zz-ss*yc;
 				};
 			}
 		};
 		sc=conj(sc);
 		if (!isreal){
 			for (i=1;i<=en;i++)		A((i-1),en-1)*=sc;
 			for (i=low;i<=igh;i++)	z((i-1),en-1)*=sc;
 		};
 	};
      // a root found
 		A((en-1),en-1)+=t;
 		w(en-1,en-1)=A((en-1),en-1);
 		en=enm1;

    };

//do if eignvectors are desired
	if(flag){
      // all roots found, backsubstitution of eigenvectors
		nor=0.0;
		for (i=1;i<=n;i++){
			for (j=i;j<=n;j++){
				nor+=AbsNorm(A((i-1),j-1));
			}
		}
		if ((n==1)||(nor==0.0)) return 0;
		for (nn=2;nn<=n;nn++){
			en=n+2-nn;
			x=w(en-1,en-1);
			A(en-1,en-1)=1;
			enm1=en-1;
			for(ii=1;ii<=enm1;ii++){
				i=en-ii;
				zz=A((i-1),en-1);
				if (i!=enm1){
					ip1=i+1;
					for (j=ip1;j<=enm1;j++)	muladd(zz,A((i-1),j-1),A((j-1),en-1)); //zz+= a((i-1),j-1)*a((j-1),en-1);
				}
				yc=x-w(i-1,i-1);
				if ((fabs(Re(yc))<machep)&&(fabs(Im(yc))<machep))	Re(yc,machep*nor);
				A((i-1),en-1)=zz/yc;
			}
		}

		enm1=n-1;
	  // eigenvectors of isolated root
		for (i=1;i<=enm1;i++){
			if ((i<low)||(i>igh)){
				for (j=i+1;j<=n;j++)	z((i-1),j-1)=A((i-1),j-1);
			}
		}
	  // multiply for eigenbase
		for (jj=low;jj<=enm1;jj++){
			j=n+low-jj;
			m= ((j-1)<igh)?j-1:igh;
			for (i=low;i<=igh;i++){
				zz=z((i-1),j-1);
				for (k=low;k<=m;k++)	muladd(zz,z((i-1),k-1),A((k-1),j-1)); //zz+=z((i-1),k-1)*a((k-1),j-1);
				z((i-1),j-1)=zz;
			}
		}
		  //normalize vectors
		for (i=0;i<n;i++){
			u=0.0;
			for (j=0;j<n;j++)	u+=AbsNorm(z(j,i));
			s=0.0;
			for (j=0;j<n;j++)	s+=square_norm(z(j,i)/u);
			s=u*std::sqrt(s);
			for (j=0;j<n;j++)	z(j,i)/=s;
		}

	}
	return 0;
}


//Routine to get the Eigen System (values and vectors)
// of a Re or complex matrix
//NOTE:: the diagonal eigenvalue matrix must be complex
//NOTE:: the Unitary eigenvector matrix must be complex...

template<class nt1, class st1, class CType>
void ComplexFullDiag(const	 _matrix<nt1, st1> &in,
						  _matrix<Complex<CType>,DiagonalMatrix> &diag,
						  _matrix<Complex<CType>, FullMatrix> &U,
						const bool destroy=false)
{
	int rs=in.rows(0);
 	if(in.type()==Mdiagonal  || in.type()==Midentity){
 		diag=in;
 		U.identity();
 		return;
 	}

 	/*if(in.type() == Mhermitian ){
 		ComplexHermitianDiag(in, diag, U);
 		return;
 	}else if(in.type() == Msymmetric){
		RealSymmetricDiag(in,diag, U);
		return;
	}*/

 	if(U.rows(0) != rs || U.cols(0) != rs) U.resize(rs,rs);
 	if(diag.rows(0) != rs || diag.cols(0) != rs) diag.resize(rs,rs);
	diag.fill(0);
	//U.fill(0);
	Vector<nt1> tmp(rs,0);
	if(!destroy){
		_matrix<nt1, FullMatrix> workspace=in;
 		corth(workspace,tmp,0,rs-1);	// To upper Hessenberg form
 		comqr3(workspace,1, rs, tmp,diag,U,1); // To diagonal form
	}else{
		_matrix<nt1, FullMatrix> &hold=const_cast<_matrix<nt1, FullMatrix> &>(in);
		corth(hold,tmp,0,rs-1);	// To upper Hessenberg form
 		comqr3(hold,1, rs, tmp,diag,U,1); // To diagonal form
	}
 	/*std::cout<<"U::"<<std::endl<<chop(U)<<std::endl
		<<"prop(U,diag)-in::"<<std::endl<<chop(U*diag*adjoint(U))-in<<std::endl
		<<"U*Udagar::"<<std::endl<<chop(U*adjoint(U))<<std::endl
		<<"evals::"<<chop(diag)<<std::endl;*/
}

//Routine to get the EIGENVALUES ONLY of a Re or complex matrix
//NOTE:: the diagonal eigenvalue matrix BETTER be complex
template<class nt1, class st1, class Ctype>
void ComplexFullEigenValues(const _matrix<nt1, st1> &in,
						  _matrix<Complex<Ctype>,DiagonalMatrix> &diag,
						const bool destroy=false)
{
	int rs=in.rows(0);
 	if(in.type()==Mdiagonal  || in.type()==Midentity){
 		diag=in;
 		return;
 	}

 	/*if(in.type() == Mhermitian ){
 		ComplexHermitianDiag(in, diag, U);
 		return;
 	}else if(in.type() == Msymmetric){
		RealSymmetricDiag(in,diag, U);
		return;
	}*/
	static _matrix<Complex<Ctype>, FullMatrix> dummyU;
 	if(diag.rows(0) != rs || diag.cols(0) != rs) diag.resize(rs,rs);
	diag.fill(0);
	//U.fill(0);
	Vector<nt1> tmp(rs,0);
	if(!destroy){
		_matrix<nt1, FullMatrix> workspace=in;
		try{
			corth(workspace,tmp,0,rs-1);	// To upper Hessenberg form
 			comqr3(workspace,1, rs, tmp,diag,dummyU,0);
		}catch(BL_exception e){
			e.print(std::cerr);
			std::cerr<<"Original Matrix: "<<in;
			BLEXCEPTION("leaving diagonalize");
		}
		
	}else{
		try{
			corth(in,tmp,0,rs-1);		// To upper Hessenberg form
			comqr3(in,1, rs, tmp,diag,dummyU,0); 	// To diagonal form (the '0' tells me NOT to solver for eigvetors)

		}catch(BL_exception e){
			e.print(std::cerr);
			std::cerr<<"Original Matrix: "<<in;
			BLEXCEPTION("leaving diagonalize");
		}
	}
		/*std::cout<<"U::"<<std::endl<<chop(U)<<std::endl
		<<"prop(U,diag)-in::"<<std::endl<<chop(U*diag*adjoint(U))-in<<std::endl
		<<"U*Udagar::"<<std::endl<<chop(U*adjoint(U))<<std::endl
		<<"evals::"<<chop(diag)<<std::endl;*/
}

//**********master eigenvalue functions
template<class nt1, class Ctype>
void eigenvalues(const	 _matrix<nt1,FullMatrix> &inmat,
					  _matrix<Complex<Ctype>, DiagonalMatrix> &dia,
					const bool destroy=false)
{
	if(inmat.rows() != inmat.cols()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	ComplexFullEigenValues(inmat,dia,destroy);
}


template<class T, class INstructure>
template<class nt1>
void _matrix<T, INstructure>::eigenvalues(_matrix<nt1, DiagonalMatrix> &dia)
{
	std::eigenvalues((*this), dia, false);
}

//***********master diag function

template<class expr1, class nt1>
void diag(const	 _matrixExpr< expr1 > in,
			 _matrix<nt1, DiagonalMatrix> &dia,
			 _matrix<nt1, FullMatrix> &U,
			const bool destroy=false)
{
	_matrix<nt1, typename expr1::structure> inmat(in);
	diag(inmat, dia, U, true);
}

template<class nt1, class st1>
void Diagonalize(const	 _matrix<nt1,st1> &inmat,
					  _matrix<nt1, DiagonalMatrix> &dia,
					  _matrix<nt1, FullMatrix> &U,
					const bool destroy=false)
{
	diag(inmat,dia,U,destroy);
}



//special cases for the diag function given number type AND matrix structure
//full complex
template<class Ctype>
void diag(const	 _matrix<Complex<Ctype>,FullMatrix> &inmat,
			 _matrix<Complex<Ctype>, DiagonalMatrix> &dia,
			 _matrix<Complex<Ctype>, FullMatrix> &U,
			 bool des=false)
{
	if(!inmat.issquare()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	ComplexFullDiag(inmat, dia, U, des);
}

/*
//full complex
template<class Ctype>
void diag(const	 _matrix<Complex<Ctype>,FullMatrix> &inmat,
			 _matrix<Complex<Ctype>, DiagonalMatrix> &dia,
			 _matrix<Complex<Ctype>, FullMatrix> &U,
			 bool des=false)
{
	diag(inmat, dia, U);
}
*/
//Re full (not possible)
template<class nt1>
void diag(const	 _matrix<nt1,FullMatrix> &inmat,
			 _matrix<nt1, DiagonalMatrix> &dia,
			 _matrix<nt1, FullMatrix> &U,
			 bool des=false)
{

	if(inmat.rows() != inmat.cols()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	BLEXCEPTION(" Full matrix should MUST be complex for general solution...")
}

//complex hermitian
template<class Ctype>
void diag(const	 _matrix<Complex<Ctype>,HermitianMatrix> &inmat,
			 _matrix<Complex<Ctype>, DiagonalMatrix> &dia,
			 _matrix<Complex<Ctype>, FullMatrix> &U,
			const bool des=false)
{
	if(!inmat.issquare()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	ComplexHermitianDiag(inmat, dia, U,des);
}

//Re Hermitian matrix
template<class nt1>
void diag(const _matrix<nt1,HermitianMatrix> &inmat,
			_matrix<nt1, DiagonalMatrix> &dia,
			_matrix<nt1, FullMatrix> &U,
			const bool des=false)
{
	if(inmat.rows() != inmat.cols()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	RealSymmetricDiag(inmat,dia,U,des);
}

//complex Symmetric (Not possible....)
template<class Ctype>
void diag(const	_matrix<Complex<Ctype>,SymmetricMatrix> &inmat,
			_matrix<Complex<Ctype>, DiagonalMatrix> &dia,
			_matrix<Complex<Ctype>, FullMatrix> &U,
			const bool des=false)
{
	if(!inmat.issquare()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	BLEXCEPTION(" Symmetric matrix should not be of complex type...use Hermitian instead")
}


// Re Symmetric
template<class nt1>
void diag(const	 _matrix<nt1,SymmetricMatrix> &inmat,
			  _matrix<nt1, DiagonalMatrix> &dia,
			  _matrix<nt1, FullMatrix> &U,
			const bool des=false)
{
	if(inmat.rows() != inmat.cols()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	RealSymmetricDiag(inmat,dia,U,des);
}

// diagonal
template<class nt1>
void diag(const	_matrix<nt1,DiagonalMatrix> &inmat,
			_matrix<nt1, DiagonalMatrix> &dia,
			_matrix<nt1, FullMatrix> &U,
			const bool des=false)
{
	if(inmat.rows() != inmat.cols()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	U.identity(inmat.rows());
	dia=inmat;
}

// identity
template<class nt1>
void diag(const 	_matrix<nt1,IdentityMatrix> &inmat,
			_matrix<nt1, DiagonalMatrix> &dia,
			_matrix<nt1, FullMatrix> &U,
			const bool des=false)
{
	if(!inmat.issquare()){
		BLEXCEPTION(" matrix must be square to diagonalize...")
	}
	U.identity(inmat.rows(), inmat.cols());
	dia=inmat;
}

//FOR SOME REASON MW Bitches about the default argument being BAD
template<class T, class INstructure>
template<class nt1>
void _matrix<T, INstructure>::diag( _matrix<nt1, DiagonalMatrix> &dia,  _matrix<nt1, FullMatrix> &U, const bool des)
{
	BL_NAMESPACE::diag((*this), dia, U, des);
}


template<class T, class INstructure>
template<class nt1>
void _matrix<T, INstructure>::diagonalize( _matrix<nt1, DiagonalMatrix> &dia,  _matrix<nt1, FullMatrix> &U, const bool des){
	BL_NAMESPACE::diag((*this), dia, U, des);
}


END_BL_NAMESPACE



#endif

