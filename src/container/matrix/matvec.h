/**************************************************
 matvec.h

 performs the matrix vector operations

 1) Matrix*Vector --> VExpr
 2) cross(Vector, Vector) --> _matrixExpr

 should be included from '_matrix.h'

 **************************************************/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-23-01
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

#ifndef _matvec_h_
#define _matvec_h_ 1


BEGIN_BL_NAMESPACE


 //************************************************************************
 //		!MATRIX VECTOR OPERATIONS! returns a VectorExpr
 // much like Mat*Mat we need a special class to handle the little
 // random sums that need to be done

 /************8 SUB MAtrix * VEctro Template functions */
template<class T, class INstructure, class T2>
SumType(OutType(T,T2)) _mul_mat_vec_(const _matrixReference<T, INstructure> &a,const Vector<T2> &v,  int row)
{
	SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	for(int k=0;k<a.cols(0);k++)	 tmp+=a(row,k) * v(k);
	return tmp;
}

template<class T, class INstructure, class expr1>
SumType(OutType(T,typename expr1::numtype))
  _mul_mat_vec_(const _matrixReference<T, INstructure> &a,const VExpr<expr1> &v,  int row)
{
	SumType(OutType(T,typename expr1::numtype))
	  tmp=ZeroType<SumType(OutType(T,typename expr1::numtype))>::zero();
	for(int k=0;k<a.cols(0);k++)	 tmp+=a(row,k) * v(k);
	return tmp;
}


template<class T, class T2>
SumType(OutType(T,T2))
  _mul_mat_vec_(const _matrixReference<T, DiagonalMatrix> &a,const Vector<T2> &v,  int row)
{
	//SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	//tmp+=a(i,k) * v(k);
	return a(row, row)*v(row);
}

template<class T, class expr1>
SumType(OutType(T,typename expr1::numtype))
  _mul_mat_vec_(const _matrixReference<T, DiagonalMatrix> &a,const VExpr<expr1> &v,  int row)
{
	//SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	//tmp+=a(i,k) * v(k);
	return a(row, row)*v(row);
}

template<class T, class T2>
SumType(OutType(T,T2)) _mul_mat_vec_(const _matrixReference<T, IdentityMatrix> &a,const Vector<T2> &v,  int row)
{
	//SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	//tmp+=a(i,k) * v(k);
	return v(row);
}

template<class T, class expr1>
SumType(OutType(T,typename expr1::numtype))
  _mul_mat_vec_(const _matrixReference<T, IdentityMatrix> &a,const VExpr<expr1> &v,  int row)
{
	//SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	//tmp+=a(i,k) * v(k);
	return v(row);
}

 /************8 SUB MAtrix * VEctro Template functions */
template<class expr, class T2>
SumType(OutType(typename expr::numtype,T2)) _mul_mat_vec_(const _matrixExpr<expr> &a,const Vector<T2> &v,  int row)
{
	SumType(OutType(typename expr::numtype,T2)) tmp=ZeroType<SumType(OutType(typename expr::numtype,T2))>::zero();
	for(int k=0;k<a.cols(0);k++)	 tmp+=a(row,k) * v(k);
	return tmp;
}

template<class expr, class expr2>
SumType(OutType(typename expr::numtype,typename expr2::numtype))
	_mul_mat_vec_(const _matrixExpr<expr> &a,const VExpr<expr2> &v,  int row)
{
	SumType(OutType(typename expr::numtype,typename expr2::numtype))
		tmp=ZeroType<SumType(OutType(typename expr::numtype,typename expr2::numtype))>::zero();
	for(int k=0;k<a.cols(0);k++)	 tmp+=a(row,k) * v(k);
	return tmp;
}


template<class T, class T2>
SumType(OutType(T,T2)) _mul_mat_vec_(const _matrixReference<T, TriDiagonalMatrix> &a,const Vector<T2> &v,  int row)
{
	//SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	//tmp+=a(i,k) * v(k);
	if(row==0 && a.rows(v.size())>1){	return a(0,0)*v(0)+a(0,1)*v(1);	}
	else if(row==0){	return a(0,0)*v(0);	}

	if(row==a.rows(v.size())-1 && a.rows(v.size())-1>=0){	return a(row,row)*v(row)+a(row,row-1)*v(row-1);	}
	else if(row==a.rows(v.size())-1){	return a(row,row)*v(row);	}

	return a(row, row-1)*v(row-1)+a(row,row)*v(row)+a(row, row+1)*v(row+1);
}


 template<class mat_, class vec_>
 class VExprBinOpMatMul{
 	private:
 		mat_ expr1_;
 		vec_ expr2_;

 		const void lenerr()const{
 			BLEXCEPTION(" matrix.cols() ==Must be== vector.length()")
 		}

 		VExprBinOpMatMul(){};

 	public:

 		typedef mat_ expr1;
 		typedef vec_ expr2;
 		typedef typename expr1::numtype numtype1;
 		typedef typename expr2::numtype numtype2;
 		typedef SumType(OutType(numtype1, numtype2)) numtype;

 		Vector<numtype> outish;

 		VExprBinOpMatMul(const VExprBinOpMatMul<mat_, vec_> &in):
 			expr1_(in.expr1_), expr2_(in.expr2_), outish(in.outish){

 			/*int i=0;
 			for(i=0;i<length(0);i++){
 				outish.put(i, mul(i));
 			}*/
 		}
 		VExprBinOpMatMul(const mat_ &A, const vec_ &B):
 			expr1_(A), expr2_(B), outish(B.length(0)){

 			int i=0;
 			for(i=0;i<length(0);i++){
 				outish.put(i, mul(i));
 			}
 		}

 		numtype operator()(int i) const{
 			return outish(i);
 		}

 		numtype operator[](int i) const{
 			return outish(i);
 		}

 		numtype mul(int i) const{
 			return _mul_mat_vec_(expr1_, expr2_, i);
 		}

		int begin() const {	return expr2_.begin();	}
		int end()	const	{	return expr2_.end();	}

		int begin(int glen) const {	return expr2_.begin(glen);	}
		int end(int gend)	const	{	return expr2_.end(gend);	}


 		int length(int guessCols)const{
 #ifdef VBoundCheck
 			if(expr1_.rows(guessCols)!=expr2_.length(guessCols)) lenerr();
 #endif
 			return expr1_.cols(guessCols);
 		}

 		int guessLen()const{
#ifdef VBoundCheck
  			if(expr1_.rows(0)!=expr2_.length(0)) lenerr();
#endif
  			return expr1_.cols(0);
 		}
 };


 //************************************
 //	the multiplication operation
 //  Matrix*vector

 template<class T, class INstructure, class T2>
 inline VExpr<VExprBinOpMatMul< _matrixReference<T, INstructure>, Vector<T2> > >
 operator*(const _matrix<T, INstructure>& m1,const Vector<T2>& v1)
 {
     typedef VExprBinOpMatMul< _matrixReference<T, INstructure>, Vector<T2> > expr;
     return VExpr<expr>(expr(m1,v1));
 }

 template<class T, class INstructure, class T2>
 inline VExpr<VExprBinOpMatMul< _matrixReference<T, INstructure>, VExpr<T2> > >
 operator*(const _matrix<T, INstructure>& m1,const VExpr<T2>& v1)
 {
     typedef VExprBinOpMatMul< _matrixReference<T, INstructure>, VExpr<T2> > expr;
     return VExpr<expr>(expr(m1,v1));
 }

 template<class expr1, class T2>
 inline VExpr<VExprBinOpMatMul< _matrixExpr<expr1>, Vector<T2> > >
 operator*(const _matrixExpr<expr1>& m1,const Vector<T2>& v1)
 {
     typedef VExprBinOpMatMul< _matrixExpr<expr1>, Vector<T2> > expr;
     return VExpr<expr>(expr(m1,v1));
 }

 template<class expr1,  class expr2>
 inline VExpr<VExprBinOpMatMul< _matrixExpr<expr1>, VExpr<expr2> > >
 operator*(const _matrixExpr<expr1>& m1,const VExpr<expr2>& v1)
 {
     typedef VExprBinOpMatMul< _matrixExpr<expr1>, VExpr<expr2> > expr;
     return VExpr<expr>(expr(m1,v1));
 }

 //************************************************************************
 //		!MATRIX VECTOR<coord<> > OPERATIONS! returns a VectorExpr
 // much like Mat*Mat we need a special class to handle the little
 // random sums that need to be done
 //
 //  NOTE:: this is a spcialization for matrix<T>*Vector<coord<T, N> >
 //         where matrix.cols()=N*length(Vector)

 /************8 SUB MAtrix * VEctro Template functions */
template<class T, class INstructure, class T2, int N>
 SumType(OutType(T,T2)) _mul_mat_vec_(const _matrixReference<T, INstructure> &a,const Vector<coord< T2, N > > &v,  int row)
{
	typedef SumType(OutType(T,T2)) numtype;
	numtype tmp=ZeroType<numtype >::zero();
	int ct=0, k, j;
	for(k=0;k<v.length(0);k++)
		for(j=0;j<N;j++)
			tmp+=a(row,ct++) * v[k][j];


	return tmp;
}


 template<class mat_, class T, int N>
 class VExprBinOpMatMulCoord{
 	private:
 		mat_ expr1_;
 		Vector<coord<T, N> >  expr2_;

 		const void lenerr()const{
 			BLEXCEPTION(" matrix.cols() ==Must be== N*vector.length()")
 		}

 		VExprBinOpMatMulCoord(){};

 	public:

 		typedef mat_ expr1;
 		typedef typename expr1::numtype numtype1;
 		typedef T numtype2;
 		typedef coord<SumType(OutType(numtype1,T)), N> numtype;
		typedef Vector<numtype > expr2;

 		Vector<numtype> outish;

 		VExprBinOpMatMulCoord(const VExprBinOpMatMulCoord<mat_, T, N > &in):
 			expr1_(in.expr1_), expr2_(in.expr2_), outish(in.outish){

 			/*int i=0;
 			for(i=0;i<length(0);i++){
 				outish.put(i, mul(i));
 			}*/
 		}

 		VExprBinOpMatMulCoord(const mat_ &A, const Vector<coord<T, N> > &B):
 			expr1_(A), expr2_(B), outish(A.rows(0)/N){
 			int i=0, j=0, ct=0;
 			for(i=0;i<outish.length(0);++i)
 				for(j=0;j<N;++j){
 					outish[i][j]=mul(ct++);
 				}

 		}

 		numtype operator()(int i) const{
 			return outish(i);
 		}

 		numtype operator[](int i) const{
 			return outish(i);
 		}

 		SumType(OutType(numtype1,T)) mul(int i) const{
 			return _mul_mat_vec_(expr1_, expr2_, i);
 		}

		int begin() const {	return expr2_.begin();	}
		int end()	const	{	return expr2_.end();	}

		int begin(int glen) const {	return expr2_.begin(glen);	}
		int end(int gend)	const	{	return expr2_.end(gend);	}


 		int length(int guessCols)const{
 #ifdef VBoundCheck
 			if(expr1_.cols(guessCols)!=N*expr2_.length(guessCols)) lenerr();
 #endif
 			return expr1_.rows(guessCols)/N;
 		}

 		int guessLen()const{
#ifdef VBoundCheck
  			if(expr1_.cols(0)!=N*expr2_.length(0)) lenerr();
#endif
  			return expr1_.rows(0)/N;
 		}
 };

 template<class T, class INstructure, class T2, int N>
 inline VExpr<VExprBinOpMatMulCoord< _matrixReference<T, INstructure>, T2, N  > >
 operator*(const _matrix<T, INstructure>& m1,const Vector<coord<T2, N> >& v1)
 {
     typedef VExprBinOpMatMulCoord< _matrixReference<T, INstructure>, T2, N > expr;
     return VExpr<expr>(expr(m1,v1));
 }


 template<class expr1, class T2, int N>
 inline VExpr<VExprBinOpMatMulCoord< _matrixExpr<expr1>, T2, N  > >
 operator*(const _matrixExpr<expr1>& m1,const Vector<coord<T2, N> >& v1)
 {
     typedef VExprBinOpMatMulCoord< _matrixExpr<expr1>,  T2, N  > expr;
     return VExpr<expr>(expr(m1,v1));
 }

template<class expr1,class expr2,  class T2, int N>
inline VExpr<VExprBinOpMatMulCoord< _matrixExpr<expr2>, T2 , N> >
operator*(const _matrixExpr<expr1>& m1,
		  const VExpr<VExprBinOpMatMulCoord< _matrixExpr<expr1>, T2, N  > > & v1)
{
 typedef VExprBinOpMatMulCoord< _matrixExpr<expr2>, T2,N > expr;
 return VExpr<expr>(expr(m1,Vector<coord<T2, N> >(v1)));
}

template<class T, class INstructure, class expr1,  class T2, int N>
inline VExpr<VExprBinOpMatMulCoord<  _matrixReference<T, INstructure>,T2, N > >
operator*(const _matrix<T,INstructure> & m1,
		  const VExpr<VExprBinOpMatMulCoord< _matrixExpr<expr1>, T2, N  >  >& v1)
{
 typedef VExprBinOpMatMulCoord< _matrixReference<T, INstructure>,  T2,N > expr;
 return VExpr<expr>(expr(m1,Vector<coord<T2, N> >(v1)));
}

template<class T, class INstructure, class Tm2, class INstructure2, class expr2,  class T2, int N>
inline VExpr<VExprBinOpMatMulCoord<  _matrixReference<T, INstructure>,T2, N > >
operator*(const _matrix<T,INstructure> & m1,
		  const VExpr<VExprBinOpMatMulCoord< _matrixReference<Tm2, INstructure2>, T2, N  >  >& v1)
{
 typedef VExprBinOpMatMulCoord< _matrixReference<T, INstructure>,  T2,N > expr;
 return VExpr<expr>(expr(m1,Vector<coord<T2, N> >(v1)));
}


/*
// MAtrix Coord functions
// Matrix*coord
 template<class T, class INstructure, class T2, int N2>
 inline VExpr<VExprBinOpMatMul< _matrixReference<T, INstructure>, VIterConst<T2> > >
 operator*(const _matrix<T, INstructure>& m1,const coord<T2,N>& v1)
 {
     typedef VExprBinOpMatMul< _matrixReference<T, INstructure>, VIterConst<T2> > expr;
     return VExpr<expr>(expr(m1,v1.begin()));
 }
*/

/***************** For Vector*Matrix...**/
 template<class vec_, class mat_>
 class VExprBinOpMatMul2{
 	private:
 		mat_ expr2_;
 		vec_ expr1_;

 		const void lenerr()const{
 			BLEXCEPTION(" matrix.rows() ==Must be== vector.length()")
 		}
 	public:

 		typedef mat_ expr2;
 		typedef vec_ expr1;
 		typedef typename expr1::numtype numtype1;
 		typedef typename expr2::numtype numtype2;
 		typedef SumType(OutType(numtype1, numtype2)) numtype;

 		Vector<numtype> outish;

 		VExprBinOpMatMul2(const VExprBinOpMatMul<vec_, mat_> &in):
 			expr1_(in.expr1_), expr2_(in.expr2_), outish(expr2_.cols(0)){

 			int i=0;
 			for(i=0;i<expr2_.cols(0);i++){
 				outish.put(i, mul(i));
 			}
 		}
 		VExprBinOpMatMul2(const vec_ &A, const mat_ &B):
 			expr1_(A), expr2_(B), outish(B.cols(0)){

 			int i=0;
 			for(i=0;i<B.cols(0);i++){
 				outish.put(i, mul(i));
 			}
 		}

 		numtype operator()(int i) const{
 			return outish(i);
 		}

 		numtype operator[](int i) const{
 			return outish(i);
 		}

 		numtype mul(int i) const{
 			numtype tmp=ZeroType<numtype>::zero();
 			for(int k=0;k<expr2_.rows(0);k++)	 tmp+=expr1_(k) * expr2_(k,i);

 			return tmp;
 		}

		int begin() const {	return expr1_.begin();	}
		int end()	const	{	return expr1_.end();	}

		int begin(int glen) const {	return expr1_.begin(glen);	}
		int end(int gend)	const	{	return expr1_.end(gend);	}


 		int length(int guessCols)const{
 #ifdef VBoundCheck
 			if(expr2_.cols(guessCols)!=outish.length(guessCols)) lenerr();
 #endif
 			return outish.length(expr2_.cols());
 		}

 		int guessLen()const{
#ifdef VBoundCheck
  			if(expr2_.cols(0)!=outish.length(0)) lenerr();
#endif
  			return outish.guessLen();
 		}
 };

 //************************************
 //	the multiplication operation
 //  Vector*matrix

 template<class T, class INstructure, class T2>
 inline VExpr<VExprBinOpMatMul2<Vector<T2>, _matrixReference<T, INstructure> > >
 operator*(const Vector<T2>& v1,const _matrix<T, INstructure>& m1)
 {
     typedef VExprBinOpMatMul2<Vector<T2> , _matrixReference<T, INstructure> > expr;
     return VExpr<expr>(expr(v1,m1));
 }

 template<class T, class INstructure, class T2>
 inline VExpr<VExprBinOpMatMul2< Vector<T2>,  _matrixReference<T, INstructure> > >
 operator*(const VExpr<T2>& v1,const _matrix<T, INstructure>& m2)
 {
     typedef VExprBinOpMatMul2<VExpr<T2>, _matrixReference<T, INstructure> > expr;
     return VExpr<expr>(expr(v1,m2));
 }


//**********************************************************
//  Vector Vector cross product to output a matrix


 template<class v1, class v2>
 class _matrixVecVecBinOp{
 	private:
 		v1 expr1_;
 		v2 expr2_;
	public:

 		typedef v1 expr1;
 		typedef v2 expr2;
 		typedef typename expr1::numtype numtype1;
 		typedef typename expr2::numtype numtype2;
 		typedef SumType(OutType(numtype1, numtype2)) numtype;
 		typedef FullMatrix structure;
		typedef FullMatrix::invers_structure invers_structure;

 		_matrixVecVecBinOp(const v1 &A, const v2 &B): expr1_(A), expr2_(B){}

 		numtype operator()(int i, int j) const{
			return expr1_(i)*expr2_(j);
 		}

 		int rows(int guessCols)const{
 			return expr1_.length(guessCols);
 		}

 		int cols(int guessCols)const{
 			return expr2_.length(guessCols);
 		}
 };

//****The Cross Operations between two vectors*****

template<class nt1, class nt2>
inline _matrixExpr< _matrixVecVecBinOp< Vector<nt1>, Vector<nt2> > >
cross(const Vector<nt1> &v1,const Vector<nt2> &v2){
	typedef _matrixVecVecBinOp< Vector<nt1>, Vector<nt2> > expr;
	return _matrixExpr<expr>(expr(v1, v2));
}

template<class nt1, class nt2>
inline _matrixExpr< _matrixVecVecBinOp< Vector<nt1>, VExpr<nt2> > >
cross(const Vector<nt1> &v1,const VExpr<nt2> &v2){
	typedef _matrixVecVecBinOp< Vector<nt1>, VExpr<nt2> > expr;
	return _matrixExpr<expr>(expr(v1, v2));
}

template<class nt1, class nt2>
inline _matrixExpr< _matrixVecVecBinOp< VExpr<nt1>, Vector<nt2> > >
cross(const VExpr<nt1> &v1,const Vector<nt2> &v2){
	typedef _matrixVecVecBinOp< VExpr<nt1>, Vector<nt2> > expr;
	return _matrixExpr<expr>(expr(v1, v2));
}

template<class nt1, class nt2>
inline _matrixExpr< _matrixVecVecBinOp< VExpr<nt1>, VExpr<nt2> > >
cross(const VExpr<nt1> &v1,const VExpr<nt2> &v2){
	typedef _matrixVecVecBinOp< VExpr<nt1>, VExpr<nt2> > expr;
	return _matrixExpr<expr>(expr(v1, v2));
}

END_BL_NAMESPACE


#endif

