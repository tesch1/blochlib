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

/*performs templated versions of the mighty
LU decomposistion stuff

*/

#ifndef _LU_h_
#define _LU_h_ 1

#include "container/Vector/Vector.h"
#include "container/rankType.h"
#include<string>


BEGIN_BL_NAMESPACE


//**********************************************************
//**********************************************************
//					EXTERNAL TO MATRIC CLASS
//**********************************************************
//**********************************************************

/*the LU decomposition
solves
	L U=A

	where L is a lower triangular matrix
	and U is an upper triangular matrix

	a key step to solveing linear equations or inverting matrices
*/

//modifies the matrix A in such that the top is the U matrix and bottom tri is the L matrix
//nicely saves storage space and makes for a faster algo
//	'a' a matrix to be decomped...MUST BE SQUARE
// 'rowperm' is a vector that stores which rows were flipped around
// 'd' ==+/-1 if even number of row flips or odd

//return 1--happy, -1--fail
//this constant is found in 'constants.cc'
extern bool LUwarn;

template<class T, class INstructure>
int _matrix<T, INstructure>::LUdecomp_( Vector<int> &rowperm, float &d){
	//test squareness
	static const float SMALL=1.e-20;
	int i,j,k, len=rows();
	T big=0, tm=0;
	Vector<T> wt(rows(), ZeroType<T>::zero()); //row scaling factors
	d=0.;

	T ss, ttmm;
	int maxi=0;

	for(i=0;i<len;i++){
		big=0;
		for(j=0;j<len;j++){
			tm=abs(this->get(i,j));
			if(tm>big) big=tm;
		}
		if(big==0){
			if(LUwarn)
			{	std::cerr<<"LUdecomp error::"<<std::endl;
				std::cerr<<" Signular matrix input"<<std::endl;
			}
			return -1;
		}
		wt[i]=1./big;


	}
	for(j=0;j<len;j++){
		for(i=0;i<j;i++){
			ss=this->get(i,j);
			for(k=0;k<i;k++)	 ss -= this->get(i,k)* this->get(k,j);

			this->put(i,j,ss);
		}
		big=0;
		for(i=j;i<len;i++){
			ss=this->get(i,j);
			for(k=0;k<j;k++) ss -= this->get(i,k)* this->get(k,j);
			this->put(i,j,ss);

			ttmm=wt[i]*abs(ss);
			if(ttmm>=big){
				big=ttmm;
				maxi=i;
			}
		}
		if(j!=maxi){
			for(k=0;k<len;k++){
				ttmm=this->get(maxi,k);
				this->put(maxi,k, this->get(j,k));
				this->put(j,k,ttmm);
			}
			d=-d;		//even or odd??
			wt[maxi]=wt[j];
		}
		rowperm[j]=maxi;
		if(this->get(j,j)==0) this->put(j,j,SMALL);		//to avoid mess with signularity later in on other algos

		if(j!=(len-1)){
			ttmm=1./this->get(j,j);
			for(i=j+1;i<len;i++) this->put(i,j, this->get(i,j)*ttmm);
		}
	}
	return 1;	//happy
}

template<class T, class INstructure>
int _matrix<T, INstructure>::LUdecomp_( Vector<int> &rowperm, float &d, _matrix<T, FullMatrix> &LU){
	LU=(*this);
	return LU.LUdecomp_(rowperm, d);
}

template<class T, class INstructure>
int _matrix<T, INstructure>::LUdecomp_(
			Vector<int> &rowperm, float &d,
			_matrix<T, FullMatrix> &L,
			_matrix<T, FullMatrix> &U
)
{
	//test squareness
	L.identity(this->rows());
	U=(*this);
	int i,j;
	U.LUdecomp_(rowperm, d);
	for(i=0;i<U.rows();++i){
		for(j=0;j<i;++j){
			L(i,j)=U(i,j);
			U(i,j)=ZeroType<T>::zero();
		}
	}

	return 1;	//happy
}


//the Linear System Solve for a TriDiagonalMatrix.....
// Solves a*sols=r
template<class T>
bool solve(_matrix<T, TriDiagonalMatrix> &a, const Vector<T> &r, Vector<T> &sols)
{
	int j, n=a.rows();
	if(r.size() != n)
	{
		std::cerr<<std::endl<<" Error:: solve(TriDiagMatrix a, Vector r, Vector sols)"<<std::endl;
		std::cerr<<" vector 'r' is NOT the same length as 'a' has rows "<<std::endl;
		return false;
	}


	if(sols.size() != r.size()) sols.resize(n);

	T bet;
	Vector<T> gam(n,0);
	if (a(0,0) == 0.0)
	{
		std::cerr<<std::endl<<" Error:: solve(TriDiagMatrix, Vector, Vector)"<<std::endl;
		std::cerr<<" First element in Tridiagonal Matrix is 0...need to recast your "<<std::endl;
		std::cerr<<" Problem...(i.e. eliminate one of the variables)"<<std::endl;
		return false;
	}
	bet=a(0,0);
	sols[0]=r[0]/bet;
	for (j=1;j<n;j++) { //Decomposition and forward substitution.
		gam[j]=a(j-1,j)/bet;
		bet=a(j,j)-a(j, j-1)*gam[j];
		if (bet == 0.0)
		{
			std::cerr<<std::endl<<" Error:: solve(TriDiagMatrix, Vector, Vector)"<<std::endl;
			std::cerr<<" Saddly...you either have a Singular matrix "<<std::endl;
			std::cerr<<" or this algorithm is not suited for none 'diagonal-dominant' problems."<<std::endl;
			return false;
		}
		sols[j]=(r[j]-a(j,j-1)*sols[j-1])/bet;
	}
	for (j=(n-2);j>=0;j--)	sols[j] -= gam[j+1]*sols[j+1];// Backsubstitution.
	return true;
}

/* Now that we have our U and L we can back subsitute to get a
a solution to linear equations
	Ax=b
	L U x=b

	here b gets replaced with the solutions 'x'
*/
template<class T, class INstructure>
template<class T1>
int _matrix<T, INstructure>::LUbackSub_(Vector<int> &rowperm, Vector<T1> &b){
	//test squareness
	if(!issquare()){
		std::cerr<<"error: LU Back Subsitution::"<<std::endl;
		std::cerr<<" input matrix MUST be square"<<std::endl;
		return -1;
	}
	int i,j, len=rows(), indx;
	T ss;
	int itest=-1;

	for(i=0;i<len;i++){
		indx=rowperm[i];
		ss=b[indx];
		b[indx]=b[i];
		if(itest!=-1){
			for(j=itest;j<i;j++) ss -= this->get(i,j)*b[j];
		}else if(ss!=0){
			itest=i;
		}
		b[i]=ss;
	}
	for(i=len-1;i>=0;i--){
		ss=b[i];
		for(j=i+1;j<len;j++) ss -= this->get(i,j)*b[j];
		b[i]=ss/this->get(i,i);
	}
	return 1;
}

//special case for a Vector<coord<T, N> >
// NOTE:: matrix size must be N*vector.size()
template<class T, class INstructure>
template<int N>
int _matrix<T, INstructure>::LUbackSub_(Vector<int> &rowperm, Vector<coord<T, N> > &bi){
	//test squareness
	if(!issquare()){
		std::cerr<<"error: LU Back Subsitution::"<<std::endl;
		std::cerr<<" input matrix MUST be square"<<std::endl;
		return -1;
	}
	int i,j, len=rows(), indx;
	T ss;
	int itest=-1;
	int vecRo=0;

	//need to create a point of T to the elements in the VEctor
	double **b; b=new T*[len];
	for(i=0;i<bi.size();++i)
		for(j=0;j<N;++j)
			b[vecRo++]=&bi[i][j];

	for(i=0;i<len;i++){
		indx=rowperm[i];
		ss=(*b[indx]);
		(*b[indx])=(*b[i]);
		if(itest!=-1){
			for(j=itest;j<i;j++) ss -= this->get(i,j)*(*b[j]);
		}else if(ss!=0){
			itest=i;
		}
		(*b[i])=ss;

	}
	for(i=len-1;i>=0;i--){
		ss=(*b[i]);
		for(j=i+1;j<len;j++)	ss -= this->get(i,j)*(*b[j]);

	}
	delete [] b;
	return 1;
}

/*
here is the full monty guy...you should ONLY call this
when you need to solve Ax=b ONLY ONCE!!! otherwise if you are
using the same matrix 'A', you will do LUdecomp over and
over on the same matrix...a silly waste of time!

***here 'a' IS NOT modified!!
but 'b' IS converted to the solution
*/

template<class nt1, class st1>
template<class T1>
int _matrix<nt1, st1>::LUsolve_( Vector<T1> &b){
	if(!issquare()){
		std::cerr<<"error: LU solve::"<<std::endl;
		std::cerr<<" input matrix MUST be square"<<std::endl;
		return -1;
	}
	_matrix<nt1, st1> tt(*this);

	Vector<int> ii(tt.rows(), 0.);
	float d=0;
	if(LUdecomp(tt, ii, d)==-1 && LUwarn){
		std::cerr<<"Error: LUsolve::"<<std::endl;
		std::cerr<<" can go no further"<<std::endl;
		return -1;
	}

	int err=tt.LUbackSub_(ii, b);

	return err;
}

template<class T, class st1>
bool solve( _matrix<T, st1> &a, const Vector<T> &r, Vector<T> &sols){
	Vector<T> sol;
	sol=r;
	LUsolve_(sol);

	return true;
}

template<class nt1, class st1>
_matrix<nt1, FullMatrix> inv(const _matrix<nt1, st1> &a)
{
	return a.inv();
}

template<class expr>
_matrix<typename expr::numtype, FullMatrix> inv(const _matrixExpr<expr> &a)
{
	_matrix<typename expr::numtype, typename expr::structure> tm(a);
	return tm.inv();
}


//**********************************************************
//**********************************************************
//					INTERNAL TO MATRIC CLASS
//**********************************************************
//**********************************************************

//DESCRUCTIVE
template<class T, class INstructure>
int _matrix<T, INstructure>::LU(){
	Vector<int> ii(rows(), 0.); //hold row permutations
	float d=0;
	int ret = LUdecomp_( ii,d);
	if(ret==-1 && LUwarn){
std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUdecomp()"<<std::endl;
std::cerr<<" the pain........"<<std::endl;
	}
	return ret;
}

template<class T, class INstructure>
Vector<int> _matrix<T, INstructure>::LUdecomp(){
	Vector<int> ii(rows(), 0.); //hold row permutations
	float d=0;
	if(LUdecomp_( ii,d)==-1 && LUwarn){
		std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUdecomp()"<<std::endl;
		std::cerr<<" the pain........"<<std::endl;
	}
	return ii;
}

//Simple NON destructive LU the returns the L & U in the same matrix
template<class T, class INstructure>
int _matrix<T, INstructure>::LU(_matrix<T, FullMatrix> &LU){
	Vector<int> ii(rows(), 0.); //hold row permutations
	float d=0;
	int err=LUdecomp_( ii,d, LU);
	if(err==-1) {
		if(LUwarn){
			std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUdecomp()"<<std::endl;
			std::cerr<<" the pain........"<<std::endl;
		}
		return 0;
	}
	return 1;
}

//Simple NON destructive LU the returns the L & U in the separate matriices
template<class T, class INstructure>
int _matrix<T, INstructure>::LU(_matrix<T, FullMatrix> &L,_matrix<T, FullMatrix> &U){
	Vector<int> ii(rows(), 0.); //hold row permutations
	float d=0;
	int err=LUdecomp_( ii,d,L,U);
	if(err==-1) {
		if(LUwarn){
			std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LU()"<<std::endl;
			std::cerr<<" the pain........"<<std::endl;
		}
		return 0;
	}
	return 1;
}

template<class T, class INstructure>
int _matrix<T, INstructure>::LUdecomp(Vector<int> &rowperm){
	rowperm.resize(rows());
/*	if(rowperm.size()!=rows()){
		std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUdecomp(Vector<int>)"<<std::endl;
		std::cerr<<" number of columns in the input vector"<<std::endl;
		std::cerr<<" must be the same as the number of rows of the matrix"<<std::endl;
		std::cerr<<" hurts don't it?........"<<std::endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}
*/
	float d=0;
	int err=LUdecomp_( rowperm,d);
	if(err==-1) {
		if(LUwarn){
			std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUdecomp()"<<std::endl;
			std::cerr<<" the pain........"<<std::endl;
		}
		return 0;
	}
	return 1;
}

template<class T, class INstructure>
int _matrix<T, INstructure>::LUdecomp(Vector<int> &rowperm, float &d){
	rowperm.resize(rows());
/*	if(rowperm.size()!=rows()){
		std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUdecomp(Vector<int>)"<<std::endl;
		std::cerr<<" number of elements in the input vector"<<std::endl;
		std::cerr<<" must be the same as the number of rows of the matrix"<<std::endl;
		std::cerr<<" hurts don't it?........"<<std::endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}
*/
	int err=LUdecomp_( rowperm,d);
	if(err==-1) {
		if(LUwarn){
			std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUdecomp()"<<std::endl;
			std::cerr<<" the pain........"<<std::endl;
		}
		return 0;
	}
	return 1;
}

template<class T, class INstructure>
template<class T1>
void _matrix<T, INstructure>::LUbackSub(Vector<int> &rowperm, Vector<T1> &b){
	if(rowperm.size()!=rows() || b.size()!=rows()){
		BLEXCEPTION(std::string(" number of elements in the input vectors")+
				std::string("\n must be the same as the number of rows of the matrix")+
				std::string("\n  hurts don't it?........"))
	}
	if(LUbackSub_( rowperm, b)==-1  && LUwarn){
		std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUbackSub()"<<std::endl;
		std::cerr<<" the pain........"<<std::endl;
	}
}

template<class T, class INstructure>
template<int N>
void _matrix<T, INstructure>::LUbackSub(Vector<int> &rowperm, Vector<coord<T, N> > &b){
	if(rowperm.size()!=rows() || b.size()!=rows()/N){
		BLEXCEPTION(std::string(" number of elements in the input vectors")+
				std::string("\n must be the same as the number of rows of the matrix")+
				std::string("\n  hurts don't it?........"))
	}
	if(LUbackSub_( rowperm, b)==-1  && LUwarn){
		std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUbackSub()"<<std::endl;
		std::cerr<<" the pain........"<<std::endl;
	}
}


template<class T, class INstructure>
template<class T1>
void _matrix<T, INstructure>::LUsolve(Vector<T1> &b){
	if(b.size()!=rows()){
		BLEXCEPTION(std::string(" number of elements in the input vectors")+
				std::string("\n must be the same as the number of rows of the matrix")+
				std::string("\n  hurts don't it?........"))
	}
	int err=LUsolve_(b);
	if(err==-1  && LUwarn){
		std::cerr<<std::endl<<"Error:  _matrix<nt1, st1>::LUsolve()"<<std::endl;
		std::cerr<<" the pain........"<<std::endl;
	}
	return err;
}


template<class T, class INstructure>
Vector<T> _matrix<T, INstructure>::solve(const Vector<T> &r)
{
	Vector<T> sol;
	if(!std::solve(*this, r, sol)  && LUwarn)
	{
		std::cerr<<std::endl<<"Error:  _matrix<T, TriDiagonalMatrix>::solve(Vector<T>)"<<std::endl;
		std::cerr<<" solve Failed..."<<std::endl;
	}
	return sol;
}





template<class T, class INstructure>
_matrix<T, typename _matrix<T, INstructure>::invers_structure> _matrix<T, INstructure>::inv() const
{
	int nr = rows();
	_matrix<T, invers_structure> theinv(nr,nr);
	Vector<int> indx(cols());
	_matrix<T,invers_structure> ALU(*this);
	ALU.LUdecomp(indx);
	Vector<T> Ii(nr, ZeroType<T>::zero());
	Vector<T> Ainvi(nr, ZeroType<T>::zero());
	int i=0, j=0;
	for(i=0; i<nr; i++)
	{
		Ii[i] = 1;
		Ainvi = Ii;
		ALU.LUbackSub(indx, Ainvi);
		for(j=0; j<nr; j++)	theinv(j,i) = Ainvi[j];
		Ii[i] = 0;
	}
	return theinv;
}

/*template<class nt1>
_matrix<nt1, DiagonalMatrix> _matrix<nt1, DiagonalMatrix>::inv() const
{
	_matrix<nt1, DiagonalMatrix> theinv(rows(), cols());
	if(!issquare()){
		std::cerr<<"error: inv()::"<<std::endl;
		std::cerr<<" input matrix MUST be square"<<std::endl;
		//return theinv;
		    BLEXCEPTION(__FILE__,__LINE__)
	}

	_matrix<nt1, DiagonalMatrix> theinv(rows(), cols());
	for(int i=0;i<rows();i++){
		theinv(i,i)=inv(get(i,i));
	}
	return theinv;
}

template<class nt1>
_matrix<nt1, IdentityMatrix> _matrix<nt1, IdentityMatrix>::inv() const
{
	_matrix<nt1, IdentityMatrix> theinv(rows(), cols());
	return theinv;
}

*/

END_BL_NAMESPACE


#endif






