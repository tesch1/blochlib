/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-08-02
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
 matmat.h

 performs various specialized matrix matrix operations

 1) cross(Matrix, Matrix)=tensor_product(Matrix, Matrix) --> _matrixExpr
 2) transpose(matrix) --> _matrixExpr
 3) adjoint(matrix) --> _matrixExpr
 4) prop(matrix, matrix)=propogate(matrix, matrix) --> _matrixExpr

 should be included from '_matrix.h'

 **************************************************/


 #ifndef _matmat_h_
 #define _matmat_h_ 1

#include "blochconfig.h"

BEGIN_BL_NAMESPACE





 //*********cross Matrix Matrix class

template<class m1, class m2>
class crossMatMatBinOp{
	private:
		m1 expr1_;
		m2 expr2_;


	public:

		typedef m1 expr1;
		typedef m2 expr2;
		typedef typename expr1::numtype numtype1;
		typedef typename expr2::numtype numtype2;
		typedef SumType(OutType(numtype1, numtype2)) numtype;
		typedef typename MatOutType(expr1, expr2) structure;
		typedef typename MatOutTypeS(typename m1::invers_structure, typename m2::invers_structure) invers_structure;

		_matrix<numtype, structure> dd;

		crossMatMatBinOp(const m1 &A, const m2 &B): expr1_(A), expr2_(B),
			dd(expr1_.rows(0)*expr2_.rows(0), expr1_.cols(0)*expr2_.cols(0))
		{
			dd.fill(ZeroType<numtype>::zero());
			calccross();
		}

		void calccross(){
			//cout<<"rs:"<<rows(expr1_.rows(0)*expr2_.rows(0))<<endl<<"cl:"<<cols(expr1_.cols(0)*expr2_.cols(0))<<endl;
			typename expr1::structure::iterator iter1(expr1_.rows(0), expr1_.cols(0));
			int rr=expr2_.rows(0);
			int cc=expr2_.cols(0);
			while(iter1){
				//typename expr2::structure::iterator iter2(expr2_.rows(0), expr2_.cols(0));
				FullMatrix::iterator iter2(expr2_.rows(0), expr2_.cols(0));
				while(iter2){
					dd.put(iter1.row()*rr+iter2.row(), iter1.col()*cc+iter2.col(),
						expr1_(iter1.row(), iter1.col()) * expr2_(iter2.row(), iter2.col()));
					++iter2;
				}
				++iter1;
			}
		}



		numtype operator()(int i, int j) const{
			return dd(i,j);
		}

		int rows(int guessCols)const{
			return dd.rows();
		}

		int cols(int guessCols)const{
			return dd.cols();
		}
 };

 //**********the cross operations******

template<class nt1, class st1, class nt2, class st2>
inline _matrixExpr< crossMatMatBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2> > >
cross(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2){
	typedef crossMatMatBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2> > expr;
	return _matrixExpr<expr>(expr(m1.getRef(), m2.getRef()));
}

template<class ex1, class nt2, class st2>
inline _matrixExpr< crossMatMatBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2> > >
cross(const _matrixExpr<ex1> &m1, const _matrix<nt2, st2> &m2){
	 typedef crossMatMatBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2> > expr;
	 return _matrixExpr<expr>(expr(m1, m2.getRef()));
}

template<class nt1, class st1, class ex2 >
inline _matrixExpr< crossMatMatBinOp<  _matrixReference<nt1, st1> , _matrixExpr<ex2> > >
cross(const _matrix<nt1, st1> &m1, const _matrixExpr<ex2> &m2){
	 typedef crossMatMatBinOp<  _matrixReference<nt1, st1> , _matrixExpr<ex2> > expr;
	 return _matrixExpr<expr>(expr(m1.getRef(), m2));
}

template<class ex1, class ex2 >
inline _matrixExpr< crossMatMatBinOp<  _matrixExpr<ex1> , _matrixExpr<ex2> > >
cross(const _matrixExpr<ex1> &m1, const _matrixExpr<ex2> &m2){
	 typedef crossMatMatBinOp<  _matrixExpr<ex1> , _matrixExpr<ex2> > expr;
	 return _matrixExpr<expr>(expr(m1, m2));
}

template<class nt1, class st1, class nt2, class st2>
inline _matrixExpr< crossMatMatBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2> > >
tensor_product(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2){
	return cross(m1, m2);
}

template<class ex1, class nt2, class st2>
inline _matrixExpr< crossMatMatBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2> > >
tensor_product(const _matrixExpr<ex1> &m1, const _matrix<nt2, st2> &m2){
	return cross(m1, m2);
}

template<class nt1, class st1, class ex2 >
inline _matrixExpr< crossMatMatBinOp<  _matrixReference<nt1, st1> , _matrixExpr<ex2> > >
tensor_product(const _matrix<nt1, st1> &m1, const _matrixExpr<ex2> &m2){
	return cross(m1, m2);
}

template<class ex1, class ex2 >
inline _matrixExpr< crossMatMatBinOp<  _matrixExpr<ex1> , _matrixExpr<ex2> > >
tensor_product(const _matrixExpr<ex1> &m1, const _matrixExpr<ex2> &m2){
	return cross(m1, m2);
}

//***************************************************
//	TRANSPOSE operation class...if request for i,j return j,i element

template<class m1>
class transposeMatUnOp{
	private:
		m1 expr1_;

	public:

		typedef m1 expr1;
		typedef typename expr1::numtype numtype;
		typedef typename m1::structure structure;
		typedef typename m1::invers_structure invers_structure;

		transposeMatUnOp(const m1 &A): expr1_(A){}

		numtype operator()(int i, int j) const{
			return expr1_(j,i);
		}

		int rows(int guessCols)const{
			return expr1_.rows(guessCols);
		}

		int cols(int guessCols)const{
			return expr1_.cols(guessCols);
		}
 };


//****************the transpose operations
template<class nt1, class st1>
inline _matrixExpr< transposeMatUnOp< _matrixReference<nt1, st1> > >
transpose(const _matrix<nt1, st1> &m1){
	typedef transposeMatUnOp< _matrixReference<nt1, st1> > expr;
	return _matrixExpr<expr>(expr(m1.getRef()));
}

template<class ex1>
inline _matrixExpr< transposeMatUnOp< _matrixExpr<ex1> > >
transpose(const _matrixExpr<ex1> &m1){
	 typedef transposeMatUnOp< _matrixExpr<ex1> > expr;
	 return _matrixExpr<expr>(expr(m1));
}

/**** INTERNAL CLASS Transpose ****/
//solves the problem assosicated with A=adjoint(A)
namespace BL_Transpose{

template<class nt1, class st1>
inline _matrixExpr< transposeMatUnOp< _matrixReference<nt1, st1> > >
transpose(const _matrix<nt1, st1> &m1){
	typedef transposeMatUnOp< _matrixReference<nt1, st1> > expr;
	return _matrixExpr<expr>(expr(m1.getRef()));
}
}

template<class T, class INstructure>
void _matrix<T, INstructure>::transpose()
{
	_matrix<T, INstructure> B=(*this);
	(*this)=BL_Transpose::transpose(B);
}

//******************************************************
//		ADJOINT operation --> transpose(conj(m))

template<class m1>
class adjointMatUnOp{
	private:
		m1 expr1_;

	public:

		typedef m1 expr1;
		typedef typename expr1::numtype numtype;
		typedef typename m1::structure structure;
		typedef typename m1::invers_structure invers_structure;


		adjointMatUnOp(const m1 &A): expr1_(A){}

		numtype operator()(int i, int j) const{
			return conj(expr1_(j,i));
		}

		int rows(int guessCols)const{
			return expr1_.rows(guessCols);
		}

		int cols(int guessCols)const{
			return expr1_.cols(guessCols);
		}
 };


//*****************adjoit operation==conjugate(transpose(m1))

template<class nt1, class st1>
inline _matrixExpr< adjointMatUnOp< _matrixReference<nt1, st1> > >
adjoint(const _matrix<nt1, st1> &m1){
	typedef adjointMatUnOp< _matrixReference<nt1, st1> > expr;
	return _matrixExpr<expr>(expr(m1.getRef()));
}

template<class ex1>
inline _matrixExpr< adjointMatUnOp< _matrixExpr<ex1> > >
adjoint(const _matrixExpr<ex1> &m1){
	 typedef adjointMatUnOp< _matrixExpr<ex1> > expr;
	 return _matrixExpr<expr>(expr(m1));
}



/**** INTERNAL CLASS AJOINT ****/
//solves the problem assosicated with A=adjoint(A)

//need to provide this namespace scope for the
// refing the adjoint external to the class
namespace BL_Adjoint{

template<class nt1, class st1>
inline _matrixExpr< adjointMatUnOp< _matrixReference<nt1, st1> > >
adjoint(const _matrix<nt1, st1> &m1){
	typedef adjointMatUnOp< _matrixReference<nt1, st1> > expr;
	return _matrixExpr<expr>(expr(m1.getRef()));
}

}

template<class T, class INstructure>
void _matrix<T, INstructure>::adjoint()
{
	_matrix<T, INstructure> B(*this);
	(*this)=BL_Adjoint::adjoint(B); //point to the 'outside class' function
}


/*
 // *********propogation Matrix Matrix class
// performs a normal QM matrix propogation
// <out> = adjoint(<U>)*<in>*<U>
// The first matrix input is U
// i.e. prop(<U>, <Mat>)
// the structure of <out> will be the structure of <in>
// !!!NOTE::!!! ************** the <U> best be UNITARY
template<class m1, class m2>
class propMatMatBinOp{
	private:
		m1 expr1_;
		m2 expr2_;
		void lenerr() const{
			std::cerr<<std::endl<<"Error: prop(U, Mat)"<<std::endl;
			std::cerr<<" U and Mat must be the same size"<<std::endl;
			    BLEXCEPTION(__FILE__,__LINE__)
		}

		propMatMatBinOp(){};

	public:

		typedef m1 expr1;
		typedef m2 expr2;
		typedef typename expr1::numtype numtype1;
		typedef typename expr2::numtype numtype2;
		typedef SumType(OutType(numtype1, numtype2)) numtype;
		typedef typename MatOutType(typename m2::structure, HermitianMatrix) structure; //seting <out> structure to the <in> sturcture
		typedef typename MatOutTypeS(typename m1::invers_structure, typename m2::invers_structure) invers_structure;


		_matrix<numtype, structure> dd;

		propMatMatBinOp( const propMatMatBinOp &in):
			expr1_(in.expr1_), expr2_(in.expr2_),
			dd(in.dd.getRef())
		{
		};

		propMatMatBinOp(const m1 &A, const m2 &B): expr1_(A), expr2_(B),
			dd(rows(0), cols(0)){
			calcprop();
		}

	// this prop is broken up into
	// dd_ij = sum_l [adjoint(U)_il sum_k[ in_lk U_kj] ]
	// which can be a bit simplified a bit, by removing the 'adjoint' op by
	// dd_ij = sum_l [conj(U_li) * sum_k[ in_lk * U_kj] ]
		void calcprop(){
			typename structure::iterator iter1(expr1_.rows(0), expr1_.cols(0));  //dd_ij
			//int rr=expr2_.rows(0);
			//int cc=expr2_.cols(0);
			numtype ddij;
			while(iter1){
				typename expr1::structure::iterator iter2(expr2_.rows(0), expr2_.cols(0)); //in_lk
				ddij=ZeroType<numtype>::zero();
				while(iter2){
					ddij+=expr1_(iter1.row(),iter2.row())*  //conj(U_li)
								expr2_(iter2.row(), iter2.col())*				//<in>_lk
								conj(expr1_(iter1.col(),iter2.col()));				//U_kj
					++iter2;
				}
				dd.put(iter1.row(), iter1.col(), ddij);
				++iter1;
			}
		}



		numtype operator()(int i, int j) const{
			return dd(i,j);
		}

		int rows(int guessRows)const{
#ifdef MatDeBug
			if(expr1_.rows(guessRows) != expr2_.rows(guessRows)) lenerr();
#endif
			return expr1_.rows(guessRows);
		}

		int cols(int guessCols)const{
#ifdef MatDeBug
			if(expr1_.cols(guessCols) != expr2_.cols(guessCols)) lenerr();
#endif
			return expr1_.cols(guessCols);
		}
 };



// **************the Propogator operations
template<class nt1, class st1, class nt2, class st2>
inline _matrixExpr< propMatMatBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2> > >
prop(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2){
	typedef propMatMatBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2> > expr;
	return _matrixExpr<expr>(expr(m1.getRef(), m2.getRef()));
}

template<class ex1, class nt2, class st2>
inline _matrixExpr< propMatMatBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2> > >
prop(const _matrixExpr<ex1> &m1, const _matrix<nt2, st2> &m2){
	 typedef propMatMatBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2> > expr;
	 return _matrixExpr<expr>(expr(m1, m2.getRef()));
}

template<class nt1, class st1, class ex2 >
inline _matrixExpr< propMatMatBinOp<  _matrixReference<nt1, st1> , _matrixExpr<ex2> > >
prop(const _matrix<nt1, st1> &m1, const _matrixExpr<ex2> &m2){
	 typedef propMatMatBinOp<  _matrixReference<nt1, st1> , _matrixExpr<ex2> > expr;
	 return _matrixExpr<expr>(expr(m1.getRef(), m2));
}

template<class ex1, class ex2 >
inline _matrixExpr< propMatMatBinOp<  _matrixExpr<ex1> , _matrixExpr<ex2> > >
prop(const _matrixExpr<ex1> &m1, const _matrixExpr<ex2> &m2){
	 typedef propMatMatBinOp<  _matrixExpr<ex1> , _matrixExpr<ex2> > expr;
	 return _matrixExpr<expr>(expr(m1, m2));
}

*/
template<class nt1, class st1, class nt2, class st2>
inline  MUL_STRUCTURE3( _matrix<nt1, st1>, _matrix<nt2, st2>)
	prop(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2)
	{	return (m1)*m2*adjoint(m1);	}

//template<class Ctype>
inline _matrix<Complex<float>, FullMatrix>
	prop(const _matrix<Complex<float>, FullMatrix> &m1, const _matrix<Complex<float>, FullMatrix> &m2)
{

#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef Complex<float> numtype;
	if(m2.rows() != m2.cols()) m1.PropErr();

	_matrix<numtype, structure> dd(m1.rows(), m1.cols(), 0.0);
	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());

	const numtype alp(1.0,0.0), beta(0.0,0.0);

	int M=m2.rows(), N=m1.cols(), K=m2.cols();
	numtype *TMMmat; TMMmat=new numtype[M*N];

	//tmm=m2*adjoint(m1)..tmm-->[M x N]=[M x M]([N x M]')
	CGEMM("N","C", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, TMMmat, &M);
	M=m1.cols(), N=m2.cols(), K=m1.rows();
	//dd=m1*tmm-->[N x N]=[N x M][M x N]
	CGEMM("N","N", &M, &N, &K,&alp, aD, &M,TMMmat, &K,&beta, ddD, &M);

	delete [] TMMmat;
	return dd;
#else
	return (m1)*m2*adjoint(m1);
#endif

}

//template<class Ctype>
inline _matrix<Complex<double>, FullMatrix>
	prop(const _matrix<Complex<double>, FullMatrix> &m1,
	     const _matrix<Complex<double>, FullMatrix> &m2)
{
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef Complex<double> numtype;
	if(m2.rows() != m2.cols()) m1.PropErr();

	_matrix<numtype, structure> dd(m1.rows(), m1.cols(), 0.0);
	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());

	const numtype alp(1.0,0.0), beta(0.0,0.0);

	int M=m2.rows(), N=m1.cols(), K=m2.cols();
	numtype *TMMmat; TMMmat=new numtype[M*N];

	//tmm=m2*adjoint(m1)..tmm-->[M x N]=[M x M]([N x M]')
	ZGEMM("N","C", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, TMMmat, &M);
	M=m1.cols(), N=m2.cols(), K=m1.rows();
	//dd=m1*tmm-->[N x N]=[N x M][M x N]
	ZGEMM("N","N", &M, &N, &K,&alp, aD, &M,TMMmat, &K,&beta, ddD, &M);

	delete [] TMMmat;
	return dd;
#else
	return (m1)*m2*adjoint(m1);
#endif
}

inline _matrix<double, FullMatrix>
	prop(const _matrix<double, FullMatrix> &m1,
	     const _matrix<double, FullMatrix> &m2)
{
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef double numtype;
	if(m2.rows() != m2.cols()) m1.PropErr();

	_matrix<numtype, structure> dd(m1.rows(), m1.cols());
	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());

	numtype alp(1.0), beta(0.0);

	int M=m2.rows(), N=m1.cols(), K=m2.cols();
	numtype *TMMmat; TMMmat=new numtype[M*N];

	//tmm=m2*adjoint(m1)..tmm-->[M x N]=[M x M]([N x M]')
	DGEMM("N","C", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, TMMmat, &M);
	M=m1.cols(), N=m2.cols(), K=m1.rows();
	//dd=m1*tmm-->[N x N]=[N x M][M x N]
	DGEMM("N","N", &M, &N, &K,&alp, aD, &M,TMMmat, &K,&beta, ddD, &M);

	delete [] TMMmat;
	return dd;
#else
	return (m1)*m2*adjoint(m1);
#endif
}

inline _matrix<float, FullMatrix>
	prop(const _matrix<float, FullMatrix> &m1,
	     const _matrix<float, FullMatrix> &m2)
{
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef float numtype;
	if(m2.rows() != m2.cols()) m1.PropErr();

	_matrix<numtype, structure> dd(m1.rows(), m1.cols());
	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());

	numtype alp(1.0), beta(0.0);

	int M=m2.rows(), N=m1.cols(), K=m2.cols();
	numtype *TMMmat; TMMmat=new numtype[M*N];

	//tmm=m2*adjoint(m1)..tmm-->[M x N]=[M x M]([N x M]')
	SGEMM("N","C", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, TMMmat, &M);
	M=m1.cols(), N=m2.cols(), K=m1.rows();
	//dd=m1*tmm-->[N x N]=[N x M][M x N]
	SGEMM("N","N", &M, &N, &K,&alp, aD, &M,TMMmat, &K,&beta, ddD, &M);

	delete [] TMMmat;
	return dd;
#else
	return (m1)*m2*adjoint(m1);
#endif
}

template<class ex1, class nt2, class st2>
inline MUL_STRUCTURE2(_matrix<nt2, st2>,_matrixExpr<ex1>)
	prop(const _matrixExpr<ex1> &m1, const _matrix<nt2, st2> &m2)
	{	return (m1)*m2*adjoint(m1);	}

template<class nt1, class st1, class ex2 >
inline MUL_STRUCTURE2(_matrix<nt1, st1>,_matrixExpr<ex2>)
	prop(const _matrix<nt1, st1> &m1, const _matrixExpr<ex2> &m2)
	{	return (m1)*m2*adjoint(m1);	}

template<class ex1, class ex2 >
inline  MUL_STRUCTURE(_matrixExpr<ex1>,_matrixExpr<ex2>)
	prop(const _matrixExpr<ex1> &m1, const _matrixExpr<ex2> &m2)
	{	return (m1)*m2*adjoint(m1);	}


/**** INTERNAL CLASS PROPAGATORS ****/

template<class T, class INstructure>
template<class Ctype>
void _matrix<T, INstructure>::simTrans(const _matrix<Complex<Ctype>, FullMatrix> &A)
{
#ifdef AUX_BLAS_LIB
	//this=A*this*adjoint(A)
	typedef FullMatrix structure;
	typedef Complex<Ctype> numtype;
	if(rows() != cols()) PropErr();

	numtype *aD=const_cast<numtype *>(A.data());
	numtype *bD=const_cast<numtype *>(data());

	const Complex<Ctype> alp(1.0), beta(0.0);

	int M=rows(), N=A.cols(), K=cols();
	numtype *TMMmat; TMMmat=new numtype[M*N];

	//this=this*adjoint(A)..tmm-->[M x N]=[M x M]([N x M]')
	DGEMM("N","C", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, TMMmat, &M);
	M=A.cols(), N=cols(), K=rows();
	//this=A*this-->[N x N]=[N x M][M x N]
	DGEMM("N","N", &M, &N, &K,&alp, aD, &M,TMMmat, &K,&beta, bD, &M);
	delete TMMmat;
#else
	(*this)= (A)*(*this)*adjoint(A);
#endif
}

template<class T, class INstructure>
void _matrix<T, INstructure>::simTrans(const _matrix<double, FullMatrix> &A)
{
	//this=A*this*adjoint(A)
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef double numtype;
	if(rows() != cols()) PropErr();

	numtype *aD=const_cast<numtype *>(A.data());
	numtype *bD=const_cast<numtype *>(data());

	double alp(1.0), beta(0.0);

	int M=rows(), N=A.cols(), K=cols();
//	numtype *TMMmat; TMMmat=new numtype[M*N];

	//this=this*adjoint(A)..tmm-->[M x N]=[M x M]([N x M]')
	DGEMM("N","C", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, bD, &M);
	M=A.cols(), N=cols(), K=rows();
	//this=A*this-->[N x N]=[N x M][M x N]
	DGEMM("N","N", &M, &N, &K,&alp, aD, &M,bD, &K,&beta, bD, &M);
#else
	(*this)= (A)*(*this)*adjoint(A);
#endif
}

template<class T, class INstructure>
void _matrix<T, INstructure>::simTrans(const _matrix<float, FullMatrix> &A)
{
	//this=A*this*adjoint(A)
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef float numtype;
	if(rows() != cols()) PropErr();

	numtype *aD=const_cast<numtype *>(A.data());
	numtype *bD=const_cast<numtype *>(data());

	numtype alp(1.0), beta(0.0);

	int M=rows(), N=A.cols(), K=cols();
//	numtype *TMMmat; TMMmat=new numtype[M*N];

	//this=this*adjoint(A)..tmm-->[M x N]=[M x M]([N x M]')
	SGEMM("N","C", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, bD, &M);
	M=A.cols(), N=cols(), K=rows();
	//this=A*this-->[N x N]=[N x M][M x N]
	SGEMM("N","N", &M, &N, &K,&alp, aD, &M,bD, &K,&beta, bD, &M);
#else
	(*this)= (A)*(*this)*adjoint(A);
#endif
}

template<class T, class INstructure>
template<class nt1, class st1>
void _matrix<T, INstructure>::simTrans(const _matrix<nt1, st1> &A)
{	simTrans(_matrix<complex, FullMatrix>(A));	}

 //same as above
template<class T, class INstructure>
template<class nt1, class st1>
void _matrix<T, INstructure>::prop(const _matrix<nt1, st1> &A)
{	simTrans(A);	}

 //same as above
template<class T, class INstructure>
template<class nt1>
void _matrix<T, INstructure>::prop(const _matrix<nt1, FullMatrix> &A)
{	simTrans(A);	}

/*

template<class nt1, class st1, class nt2, class st2>
inline  MUL_STRUCTURE3( _matrix<nt1, st1>, _matrix<nt2, st2>)
	prop(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2)
{
	typename st1::iterator iter1(m1.rows(0), m1.cols(0));  //dd_ij
	//int rr=expr2_.rows(0);
	//int cc=expr2_.cols(0);
	MUL_STRUCTURE3( _matrix<nt1, st1>, _matrix<nt2, st2>) dd(m1.rows(0), m1.cols(0));
	typedef OutType(nt1, nt2) numtype ;
	numtype ddij;
	while(iter1){
		typename st2::iterator iter2(m2.rows(0), m2.cols(0)); //in_lk
		ddij=ZeroType<OutType(nt1, nt2)>::zero();
		while(iter2){
			ddij+=m1(iter1.row(),iter2.row())*  //conj(U_li)
						m2(iter2.row(), iter2.col())*				//<in>_lk
						conj(m1(iter1.col(),iter2.col()));				//U_kj
			++iter2;
		}
		dd.put(iter1.row(), iter1.col(), ddij);
		++iter1;
	}
	return dd;
//	return m1*m2*adjoint(m1);
}

template<class ex1, class nt2, class st2>
inline MUL_STRUCTURE2(_matrix<nt2, st2>,_matrixExpr<ex1>)
	prop(const _matrixExpr<ex1> &m1, const _matrix<nt2, st2> &m2)
{
	typename _matrixExpr<ex1>::structure iter1(m1.rows(0), m1.cols(0));  //dd_ij
	//int rr=expr2_.rows(0);
	//int cc=expr2_.cols(0);
	MUL_STRUCTURE2(_matrix<nt2, st2>,_matrixExpr<ex1>) dd(m1.rows(0), m1.cols(0));
	typedef OutType(typename _matrixExpr<ex1>::numtype, nt2) numtype;
	numtype ddij;
	while(iter1){
		typename st2::structure::iterator iter2(m2.rows(0), m2.cols(0)); //in_lk
		ddij=ZeroType<numtype>::zero();
		while(iter2){
			ddij+=m1(iter1.row(),iter2.row())*  //conj(U_li)
						m2(iter2.row(), iter2.col())*				//<in>_lk
						conj(m1(iter1.col(),iter2.col()));				//U_kj
			++iter2;
		}
		dd.put(iter1.row(), iter1.col(), ddij);
		++iter1;
	}
	return dd;

	//return m1*m2*adjoint(m1);
}

template<class nt1, class st1, class ex2 >
inline MUL_STRUCTURE2(_matrix<nt1, st1>,_matrixExpr<ex2>)
	prop(const _matrix<nt1, st1> &m1, const _matrixExpr<ex2> &m2)
{
	typename st1::iterator iter1(m1.rows(0), m1.cols(0));  //dd_ij
	//int rr=expr2_.rows(0);
	//int cc=expr2_.cols(0);
	MUL_STRUCTURE2(_matrix<nt1, st1>,_matrixExpr<ex2>) dd(m1.rows(0), m1.cols(0));
	typedef OutType(nt1, typename _matrixExpr<ex2>::numtype) numtype;
	numtype ddij;
	while(iter1){
		typename _matrixExpr<ex2>::structure::iterator iter2(m2.rows(0), m2.cols(0)); //in_lk
		ddij=ZeroType<numtype>::zero();
		while(iter2){
			ddij+=m1(iter1.row(),iter2.row())*  //conj(U_li)
						m2(iter2.row(), iter2.col())*				//<in>_lk
						conj(m1(iter1.col(),iter2.col()));				//U_kj
			++iter2;
		}
		dd.put(iter1.row(), iter1.col(), ddij);
		++iter1;
	}
	return dd;
	//return m1*m2*adjoint(m1);
}

template<class ex1, class ex2 >
inline  MUL_STRUCTURE(_matrixExpr<ex1>,_matrixExpr<ex2>)
	prop(const _matrixExpr<ex1> &m1, const _matrixExpr<ex2> &m2)
{
	typename _matrixExpr<ex1>::iterator iter1(m1.rows(0), m1.cols(0));  //dd_ij
	//int rr=expr2_.rows(0);
	//int cc=expr2_.cols(0);
	MUL_STRUCTURE(_matrixExpr<ex1>,_matrixExpr<ex2>) dd(m1.rows(0), m1.cols(0));
	typedef OutType( typename _matrixExpr<ex1>::numtype, typename _matrixExpr<ex2>::numtype) numtype;
	 numtype ddij;
	while(iter1){
		typename _matrixExpr<ex2>::structure::iterator iter2(m2.rows(0), m2.cols(0)); //in_lk
		ddij=ZeroType<numtype>::zero();
		while(iter2){
			ddij+=m1(iter1.row(),iter2.row())*  //conj(U_li)
						m2(iter2.row(), iter2.col())*				//<in>_lk
						conj(m1(iter1.col(),iter2.col()));				//U_kj
			++iter2;
		}
		dd.put(iter1.row(), iter1.col(), ddij);
		++iter1;
	}
	return dd;

	//return m1*m2*adjoint(m1);
}

*/
/*


// *********Adjoint--propogation Matrix Matrix class
// performs a normal QM matrix propogation
// <out> = <U>*<in>*adjoint(<U>)
// The first matrix input is U
// i.e. prop(<U>, <Mat>)
// the structure of <out> will be the structure of <in>
// !!!NOTE::!!! ************** the <U> best be UNITARY
template<class m1, class m2>
class adjpropMatMatBinOp{
	private:
		m1 expr1_;
		m2 expr2_;
		void lenerr() const{
			std::cerr<<std::endl<<"Error: prop(U, Mat)"<<std::endl;
			std::cerr<<" U and Mat must be the same size"<<std::endl;
			    BLEXCEPTION(__FILE__,__LINE__)
		}

	public:

		typedef m1 expr1;
		typedef m2 expr2;
		typedef typename expr1::numtype numtype1;
		typedef typename expr2::numtype numtype2;
		typedef SumType(OutType(numtype1, numtype2)) numtype;
		typedef typename MatOutType(typename m2::structure, HermitianMatrix) structure; //seting <out> structure to the <in> sturcture
		typedef typename MatOutTypeS(typename m1::invers_structure, typename m2::invers_structure) invers_structure;


		_matrix<numtype, structure> dd;

		adjpropMatMatBinOp( const adjpropMatMatBinOp &in):
			expr1_(in.expr1_), expr2_(in.expr2_), dd(in.dd.getRef())
		{};

		adjpropMatMatBinOp(const m1 &A, const m2 &B): expr1_(A), expr2_(B),
			dd(rows(0), cols(0)){
			calcprop();
		}

	// this prop is broken up into
	// dd_ij = sum_l [adjoint(U)_il sum_k[ in_lk U_kj] ]
	// which can be a bit simplified a bit, by removing the 'adjoint' op by
	// dd_ij = sum_l [conj(U_li) * sum_k[ in_lk * U_kj] ]
		void calcprop(){
			typename structure::iterator iter1(expr1_.rows(0), expr1_.cols(0));  //dd_ij
			//int rr=expr2_.rows(0);
			//int cc=expr2_.cols(0);
			numtype ddij;
			while(iter1){
				typename expr1::structure::iterator iter2(expr2_.rows(0), expr2_.cols(0)); //in_lk
				ddij=ZeroType<numtype>::zero();
				while(iter2){
					ddij+=conj(expr1_(iter1.row(),iter2.row()))*  //conj(U_li)
								expr2_(iter2.row(), iter2.col())*				//<in>_lk
								expr1_(iter1.col(),iter2.col());				//U_kj
					++iter2;
				}
				dd.put(iter1.row(), iter1.col(), ddij);
				++iter1;
			}
		}



		numtype operator()(int i, int j) const{
			return dd(i,j);
		}

		int rows(int guessRows)const{
#ifdef MatDeBug
			if(expr1_.rows(guessRows) != expr2_.rows(guessRows)) lenerr();
#endif
			return expr1_.rows(guessRows);
		}

		int cols(int guessCols)const{
#ifdef MatDeBug
			if(expr1_.cols(guessCols) != expr2_.cols(guessCols)) lenerr();
#endif
			return expr1_.cols(guessCols);
		}
 };

// **************the Propogator operations
template<class nt1, class st1, class nt2, class st2>
inline _matrixExpr< adjpropMatMatBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2> > >
adjprop(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2){
	typedef adjpropMatMatBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2> > expr;
	return _matrixExpr<expr>(expr(m1.getRef(), m2.getRef()));
}

template<class ex1, class nt2, class st2>
inline _matrixExpr< adjpropMatMatBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2> > >
adjprop(const _matrixExpr<ex1> &m1, const _matrix<nt2, st2> &m2){
	 typedef adjpropMatMatBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2> > expr;
	 return _matrixExpr<expr>(expr(m1, m2.getRef()));
}

template<class nt1, class st1, class ex2 >
inline _matrixExpr< adjpropMatMatBinOp<  _matrixReference<nt1, st1> , _matrixExpr<ex2> > >
adjprop(const _matrix<nt1, st1> &m1, const _matrixExpr<ex2> &m2){
	 typedef adjpropMatMatBinOp<  _matrixReference<nt1, st1> , _matrixExpr<ex2> > expr;
	 return _matrixExpr<expr>(expr(m1.getRef(), m2));
}

template<class ex1, class ex2 >
inline _matrixExpr< adjpropMatMatBinOp<  _matrixExpr<ex1> , _matrixExpr<ex2> > >
adjprop(const _matrixExpr<ex1> &m1, const _matrixExpr<ex2> &m2){
	 typedef adjpropMatMatBinOp<  _matrixExpr<ex1> , _matrixExpr<ex2> > expr;
	 return _matrixExpr<expr>(expr(m1, m2));
}
*/
/*template<class nt1,  class nt2>
inline  MUL_STRUCTURE3( _matrix<nt1, FullMatrix>, _matrix<nt2, FullMatrix>)
	adjprop(const _matrix<nt1, FullMatrix> &m1, const _matrix<nt2, FullMatrix> &m2)
{
	typedef OutType(nt1, nt2) numtype;
	if(m1.rows() != m2.cols()) m1.MulErr();
	_matrix<numtype, FullMatrix> dd(m1.rows(), m1.cols());

	static int i,j,k,p;
	for(i=0;i<dd.rows();++i){
		for(j=0;j<dd.cols();++j){
			for(k=0; k<dd.rows();++k){
				dd(i,j)=conj(m1(k,i)) * m2(k,0)* m1(0, j);
				for(p=1;p<dd.cols();++k){
					dd(i,j)+=conj(m1(k,i)) * m2(k,p)* m1(p, j);
				}
			}
		}
	}

	return dd;

}
*/
template<class nt1, class st1, class nt2, class st2>
inline  MUL_STRUCTURE3( _matrix<nt1, st1>, _matrix<nt2, st2>)
	adjprop(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2)
	{	return adjoint(m1)*m2*(m1);	}

template<class ex1, class nt2, class st2>
inline MUL_STRUCTURE2(_matrix<nt2, st2>,_matrixExpr<ex1>)
	adjprop(const _matrixExpr<ex1> &m1, const _matrix<nt2, st2> &m2)
	{	return adjoint(m1)*m2*(m1);	}

template<class Ctype>
inline _matrix<Complex<Ctype>, FullMatrix>
	adjprop(const _matrix<Complex<Ctype>, FullMatrix> &m1,
	        const _matrix<Complex<Ctype>, FullMatrix> &m2)
{
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef Complex<Ctype> numtype;
	if(m2.rows() != m2.cols()) m1.PropErr();

	_matrix<Complex<Ctype>, structure> dd(m1.rows(), m1.cols());
	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());

	Complex<Ctype> alp(1.0,0.0), beta(0.0,0.0);

	int M=m2.rows(), N=m1.cols(), K=m2.cols();
	numtype *TMMmat; TMMmat=new numtype[M*N];

	//tmm=m2*m1..tmm-->[M x N]=[M x M]([M x N])
	DGEMM("N","N", &M, &N, &K,&alp, bD, &M,	aD, &K,&beta, TMMmat, &M);
	M=m1.cols(), N=m2.cols(), K=m1.rows();

	//dd=adjoint(m1)*tmm-->[N x N]=[N x M][M x N]
	DGEMM("C","N", &M, &N, &K,&alp, aD, &K,TMMmat, &K,&beta, ddD, &M);

	delete [] TMMmat;
	return dd;
#else
	return adjoint(m1)*m2*(m1);
#endif
}

inline _matrix<double, FullMatrix>
	adjprop(const _matrix<double, FullMatrix> &m1,
	        const _matrix<double, FullMatrix> &m2)
{
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef double numtype;
	if(m2.rows() != m2.cols()) m1.PropErr();

	_matrix<numtype, structure> dd(m1.rows(), m1.cols());
	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());

	double alp=1.0, beta=0.0;

	int M=m2.rows(), N=m1.cols(), K=m2.cols();
	numtype *TMMmat; TMMmat=new numtype[M*N];

	//tmm=m2*adjoint(m1)..tmm-->[M x N]=[M x M]([M x N])
	DGEMM("N","N", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, TMMmat, &M);
	M=m1.cols(), N=m2.cols(), K=m1.rows();
	//dd=m1*tmm-->[N x N]=([N x M]')[M x N]
	DGEMM("C","N", &M, &N, &K,&alp, aD, &M,TMMmat, &K,&beta, ddD, &M);

	delete [] TMMmat;
	return dd;
#else
	return adjoint(m1)*m2*(m1);
#endif
}

inline _matrix<float, FullMatrix>
	adjprop(const _matrix<float, FullMatrix> &m1,
	        const _matrix<float, FullMatrix> &m2)
{
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef float numtype;
	if(m2.rows() != m2.cols()) m1.PropErr();

	_matrix<numtype, structure> dd(m1.rows(), m1.cols());
	numtype *ddD=const_cast<numtype *>(dd.data());
	numtype *aD=const_cast<numtype *>(m1.data());
	numtype *bD=const_cast<numtype *>(m2.data());

	float alp=1.0, beta=0.0;

	int M=m2.rows(), N=m1.cols(), K=m2.cols();
	numtype *TMMmat; TMMmat=new numtype[M*N];

	//tmm=m2*adjoint(m1)..tmm-->[M x N]=[M x M]([M x N])
	SGEMM("N","N", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, TMMmat, &M);
	M=m1.cols(), N=m2.cols(), K=m1.rows();
	//dd=m1*tmm-->[N x N]=([N x M]')[M x N]
	SGEMM("C","N", &M, &N, &K,&alp, aD, &M,TMMmat, &K,&beta, ddD, &M);

	delete [] TMMmat;
	return dd;
#else
	return adjoint(m1)*m2*(m1);
#endif
}

template<class nt1, class st1, class ex2 >
inline MUL_STRUCTURE2(_matrix<nt1, st1>,_matrixExpr<ex2>)
	adjprop(const _matrix<nt1, st1> &m1, const _matrixExpr<ex2> &m2)
	{	return adjoint(m1)*m2*(m1);	}

template<class ex1, class ex2 >
inline  MUL_STRUCTURE(_matrixExpr<ex1>,_matrixExpr<ex2>)
	adjprop(const _matrixExpr<ex1> &m1, const _matrixExpr<ex2> &m2)
	{	return adjoint(m1)*m2*(m1);	}

/******* INTERNAL MATRIX BITS...these save a memory copy and
* should be a bit faster then the obove ones
*/
//this is B=adjoint(A)*B*A-->similarity tranformations
template<class T, class INstructure>
template<class nt1>
void _matrix<T, INstructure>::adjprop(const _matrix<nt1, FullMatrix> &A)
{
#ifdef AUX_BLAS_LIB
	typedef FullMatrix structure;
	typedef double numtype;
	if(rows() != cols()) PropErr();

	numtype *aD=const_cast<numtype *>(A.data());
	numtype *bD=const_cast<numtype *>(data());

	static double alp(1.0), beta(0.0);

	int M=rows(), N=A.cols(), K=cols();
	//numtype *TMMmat; TMMmat=new numtype[M*N];

	//tmm=m2*adjoint(m1)..tmm-->[M x N]=[M x M]([N x M]')
	DGEMM("N","N", &M, &N, &K,&alp, bD, &M,	aD, &N,&beta, bD, &M);
	M=cols(), N=A.cols(), K=rows();
	//dd=m1*tmm-->[N x N]=([N x M]')[M x N]
	DGEMM("C","N", &M, &N, &K,&alp, aD, &M,bD, &K,&beta, bD, &M);
#else
	(*this)=adjoint(A)*(*this)*(m1);
#endif
}

/*

template<class nt1, class st1, class nt2, class st2>
inline  MUL_STRUCTURE3( _matrix<nt1, st1>, _matrix<nt2, st2>)
	adjprop(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2)
{
	typename st1::iterator iter1(m1.rows(0), m1.cols(0));  //dd_ij
	//int rr=expr2_.rows(0);
	//int cc=expr2_.cols(0);
	MUL_STRUCTURE3( _matrix<nt1, st1>, _matrix<nt2, st2>) dd(m1.rows(0), m1.cols(0));
	typedef OutType(nt1, nt2) numtype ;
	numtype ddij;
	while(iter1){
		typename st2::iterator iter2(m2.rows(0), m2.cols(0)); //in_lk
		ddij=ZeroType<OutType(nt1, nt2)>::zero();
		while(iter2){
			ddij+=conj(m1(iter1.row(),iter2.row()))*  //conj(U_li)
						m2(iter2.row(), iter2.col())*				//<in>_lk
						(m1(iter1.col(),iter2.col()));				//U_kj
			++iter2;
		}
		dd.put(iter1.row(), iter1.col(), ddij);
		++iter1;
	}
	return dd;
//	return m1*m2*adjoint(m1);
}

template<class ex1, class nt2, class st2>
inline MUL_STRUCTURE2(_matrix<nt2, st2>,_matrixExpr<ex1>)
	adjprop(const _matrixExpr<ex1> &m1, const _matrix<nt2, st2> &m2)
{
	typename _matrixExpr<ex1>::structure iter1(m1.rows(0), m1.cols(0));  //dd_ij
	//int rr=expr2_.rows(0);
	//int cc=expr2_.cols(0);
	MUL_STRUCTURE2(_matrix<nt2, st2>,_matrixExpr<ex1>) dd(m1.rows(0), m1.cols(0));
	typedef OutType(typename _matrixExpr<ex1>::numtype, nt2) numtype;
	numtype ddij;
	while(iter1){
		typename st2::structure::iterator iter2(m2.rows(0), m2.cols(0)); //in_lk
		ddij=ZeroType<numtype>::zero();
		while(iter2){
			ddij+=conj(m1(iter1.row(),iter2.row()))*  //conj(U_li)
						m2(iter2.row(), iter2.col())*				//<in>_lk
						(m1(iter1.col(),iter2.col()));				//U_kj
			++iter2;
		}
		dd.put(iter1.row(), iter1.col(), ddij);
		++iter1;
	}
	return dd;

	//return m1*m2*adjoint(m1);
}

template<class nt1, class st1, class ex2 >
inline MUL_STRUCTURE2(_matrix<nt1, st1>,_matrixExpr<ex2>)
	adjprop(const _matrix<nt1, st1> &m1, const _matrixExpr<ex2> &m2)
{
	typename st1::iterator iter1(m1.rows(0), m1.cols(0));  //dd_ij
	//int rr=expr2_.rows(0);
	//int cc=expr2_.cols(0);
	MUL_STRUCTURE2(_matrix<nt1, st1>,_matrixExpr<ex2>) dd(m1.rows(0), m1.cols(0));
	typedef OutType(nt1, typename _matrixExpr<ex2>::numtype) numtype;
	numtype ddij;
	while(iter1){
		typename _matrixExpr<ex2>::structure::iterator iter2(m2.rows(0), m2.cols(0)); //in_lk
		ddij=ZeroType<numtype>::zero();
		while(iter2){
			ddij+=conj(m1(iter1.row(),iter2.row()))*  //conj(U_li)
						m2(iter2.row(), iter2.col())*				//<in>_lk
						(m1(iter1.col(),iter2.col()));				//U_kj
			++iter2;
		}
		dd.put(iter1.row(), iter1.col(), ddij);
		++iter1;
	}
	return dd;
	//return m1*m2*adjoint(m1);
}

template<class ex1, class ex2 >
inline  MUL_STRUCTURE(_matrixExpr<ex1>,_matrixExpr<ex2>)
	adjprop(const _matrixExpr<ex1> &m1, const _matrixExpr<ex2> &m2)
{
	typename _matrixExpr<ex1>::iterator iter1(m1.rows(0), m1.cols(0));  //dd_ij
	//int rr=expr2_.rows(0);
	//int cc=expr2_.cols(0);
	MUL_STRUCTURE(_matrixExpr<ex1>,_matrixExpr<ex2>) dd(m1.rows(0), m1.cols(0));
	typedef OutType( typename _matrixExpr<ex1>::numtype, typename _matrixExpr<ex2>::numtype) numtype;
	 numtype ddij;
	while(iter1){
		typename _matrixExpr<ex2>::structure::iterator iter2(m2.rows(0), m2.cols(0)); //in_lk
		ddij=ZeroType<numtype>::zero();
		while(iter2){
			ddij+=conj(m1(iter1.row(),iter2.row()))*  //conj(U_li)
						m2(iter2.row(), iter2.col())*				//<in>_lk
						(m1(iter1.col(),iter2.col()));				//U_kj
			++iter2;
		}
		dd.put(iter1.row(), iter1.col(), ddij);
		++iter1;
	}
	return dd;

	//return m1*m2*adjoint(m1);
}

*/
/***************************** TRACE OPERATIONS **********************/

/* Single matrix trace

	trace(U)=Sum <i|U|i>

*/
//matrix
template<class nt1, class st1>
inline SumType(nt1)
trace(const _matrix<nt1, st1> &m1){
	if(m1.rows() != m1.cols()){
		BLEXCEPTION(" Matrix must be square...")
	}
	SumType(nt1) tm=ZeroType<SumType(nt1)>::zero();
	for(int i=0;i<m1.rows();++i){
		tm+=m1(i,i);
	}
	return tm;
}

//matrix Expr
template<class expr>
inline SumType(typename expr::numtype)
trace(const _matrixExpr<expr> &m1){
	if(m1.rows(0) != m1.cols(0)){
		BLEXCEPTION(" Matrix must be square...")
	}
	SumType(typename expr::numtype) tm=ZeroType<SumType(typename expr::numtype) >::zero();
	for(int i=0;i<m1.rows(0);++i){
		tm+=m1(i,i);
	}
	return tm;
}

/* Double matrix traces
	trace(U*F)=Sum_i Sum_j <i|U|j><j|F|i>
*/

//matrix matrix
template<class nt1, class st1, class nt2, class st2>
inline SumType(OutType(nt1, nt2))
trace(const _matrix<nt1, st1> &m1, const _matrix<nt2, st2> &m2){
	if(m1.rows() !=m2.cols() || m1.cols() != m2.rows()){
		BLEXCEPTION(" Matrix must be square...")
	}
	SumType(OutType(nt1, nt2)) tm=ZeroType<SumType(SumType(OutType(nt1, nt2)))>::zero();
	//typename MatOutType(st1, st2) structure;
	for(int i=0;i<m1.rows();++i){
		for(int j=0;j<m1.cols();++j){
			tm+=m1(i,j)*m2(j, i);
		}
	}
	return tm;
}

//matrix matrixExpr
template<class nt1, class st1, class expr>
inline SumType(OutType(nt1,typename expr::numtype))
trace(const _matrix<nt1, st1> &m1, const _matrixExpr<expr> &m2){
	if(m1.rows() !=m2.cols(0) || m1.cols() != m2.rows(0)){
		BLEXCEPTION(" Matrix must be square...")
	}
	SumType(OutType(nt1, typename expr::numtype)) tm=ZeroType<SumType(SumType(OutType(nt1, typename expr::numtype)))>::zero();
	//typename MatOutType(st1, st2) structure;
	for(int i=0;i<m1.rows();++i){
		for(int j=0;j<m1.cols();++j){
			tm+=m1(i,j)*m2(j, i);
		}
	}
	return tm;
}

//matrixExpr matrix
template<class nt1, class st1, class expr>
inline SumType(OutType(nt1,typename expr::numtype))
trace(const _matrixExpr<expr> &m1, const _matrix<nt1, st1> &m2){
	if(m1.rows(0) !=m2.cols() || m1.cols(0) != m2.rows()){
		BLEXCEPTION(" Matrix must be square...")
	}
	SumType(OutType(nt1, typename expr::numtype)) tm=ZeroType<SumType(OutType(nt1, typename expr::numtype))>::zero();
	//typename MatOutType(st1, st2) structure;
	for(int i=0;i<m2.cols();++i){
		for(int j=0;j<m2.rows();++j){
			tm+=m1(i,j)*m2(j, i);
		}

	}
	return tm;
}

//matrixExpr matrixExpr
template<class expr1, class expr2>
inline SumType(OutType(typename expr1::numtype,typename expr2::numtype))
trace(const _matrixExpr<expr1> &m1,const _matrixExpr<expr2> &m2){
	if(m1.rows(0) !=m2.cols(0) || m1.cols(0) != m2.rows(0)){
		BLEXCEPTION(" Matrix must be square...")
	}
	SumType(OutType(typename expr2::numtype, typename expr1::numtype)) tm=ZeroType<SumType(SumType(OutType(typename expr2::numtype, typename expr1::numtype)))>::zero();
	//typename MatOutType(st1, st2) structure;
	for(int i=0;i<m1.rows(0);++i){
		for(int j=0;j<m1.cols(0);++j){
			tm+=m1(i,j)*m2(j, i);
		}

	}
	return tm;
}



END_BL_NAMESPACE


#endif


