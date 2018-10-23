




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

/****************************************************************************
 * _matrix.h
 	A really nice expression templated class for matricies
 	this is probably as fast as C++ can get for a general matrices of objects
 	uses template expressions to expand all operations, removes all
 	temporary objects, and excutes ONE loop, for each big long
 	line of operations

 	based roughly on the excelent teachings of blitz++
 	(http://oonumerics.org/blitz)

 	if you never seen expression templates, this code will
 	look like a pile of raging letters with little meaning
 	the coding is not intuitive.

 	**note1a:: becuase matrix multiplication inherently involves summation
 	over previous elements...this particular operation cannot be
 	expanded into the 'one loop' instead a temporary must be created.
 	it does, however, use the 'matrix form' to speed up all operations
 	so if our matrix was diagonal, then only the i==j elements would ever
 	be itterated.


 	**note2::there is NO automatic matrix type conversion (or i should say it is not
 	smart about type conversion).  whatever is on the right hand side of the expression
 	std::string is always the outtype..for instance
 		<full>=<symetric> op <diagonal>...
 	will result in a  FULL matrix not the correct 'symmetric' type...one can use this method
 	to manualy type convert any series of operations, but be carefull becuase
 		<diagonal>=<full> op <full>
 	will always give you a diagonal matrix back

 	**note3:: Matrix multiplication creats the CORRECT output matrix type(using the methods found
 	in 'genmatrix.h') but the final type is still at the whim of the right-hand side of the assignment
 	so
 		<diagonal>=<full>*<symmetrical>
 	would initially give
 		<diagonal>=<full>
 	but upon assignment it would take only its diagonal elements

 	**note4:: numerical types should be auto-up converted
 	so
 		<double>*<complex>*<int>
 	would give us a
 		<complex>
 	BUT, again if the rhs is NOT the correct out put, it will try to 'down-convert' it
 	this will work for double-->int, but complex-->double is NOT defined, and you will get
 	a compiler error


 */

/******NOTE NOTE NOTE NOTE NOTE *****
  While Writing code you should enable the "MatDeBug" flag
  it turns on the matrix length and structure matching errors
  It uses ALOT of CPU cycles to check every time
  SO turn if OFF when your program is finished

*/



#ifndef __matrix_h_
#define __matrix_h_ 1

#ifndef ON_WINDOWS
#include "blochconfig.h"
#endif

#include<math.h>
#include "container/rankType.h"
#include<string>
#include "container/complex.h"
#include "container/MemChunk.h"
#include<iostream>
#include<iomanip>
#include "container/Vector/Vector.h"
//#include "container/grids/coords.h"
#include "container/range.h"

#include "container/matrix/imatrix.h"
#include "container/matrix/dmatrix.h"
#include "container/matrix/fmatrix.h"
#include "container/matrix/smatrix.h"
#include "container/matrix/hmatrix.h"
#include "container/matrix/trimatrix.h"

#include "container/matrix/matrixconfig.h"

#include "container/operations.h"  		//contains the applicative template operations like Add, Mul, Div, etc


BEGIN_BL_NAMESPACE

//the 'master' _matrix class that deals with all the
// 'type structures' above

template<class T, int N>
class coord;

template<class T>	class _matrixExpr;
template<class T>	class MIter;
template<class T>	class MIterConst;
template<class T, class st> class _matrixReference;
template<class nt2> class _matrixExprConst;


template<class T, class INstructure>
class _matrix : protected MemChunkRef<T>{
	private:
		INstructure thetype;
		void structureError(std::string oo)  {
			std::string mess="Error in matrix assignment \n"+
			" The left hand side matrix is of type \""+MatName(type())+"\"\n"+
			" But the right is of 'greater' type \""+oo+"\"\n"+
			" cannot do <"+MatName(type())+"> = <"+oo+">\n"+
			" without loss of information\n";
			BLEXCEPTION(mess)
		}
		void ColErr() const {
			BLEXCEPTION(" column desired does not exsists...")
		}

		const void RowErr() const {
			BLEXCEPTION(" row desired does not exsists...")
		}

		const void AssErr() const {
			BLEXCEPTION(" Matrix assignment not possible...mis match in sizes")
		}


		//matrix_& operator=(_matrix<nt1, st1> mat){}
	public:
		const void MulErr() const {
			BLEXCEPTION(" size mismatch  m1.rows()==m2.cols()")
		}

		const void PropErr() const {
			BLEXCEPTION(" size mismatch  m1.rows()==m1.cols()")
		}


    	typedef T        numtype;
    	typedef typename INstructure::iterator        iterator;
    	typedef INstructure structure;
    	typedef typename INstructure::invers_structure invers_structure;
		typedef _matrix<numtype, structure>   matrix_;

		/***************************************/
		// Constructors


		_matrix():
			MemChunkRef<T>(), thetype(0,0){}

		_matrix(int rows, int cols):
			MemChunkRef<T>(), thetype(rows, cols)
		{
			MemChunkRef<T>::newBlock(thetype.numelements());
		}

		//copy
		_matrix(const _matrix &in) :
			MemChunkRef<T>(), thetype(in.rows(), in.cols())
		{
			MemChunkRef<T>::newBlock(thetype.numelements());
			typename structure::iterator iter(rows(), cols());
			while (iter){
				put(iter.row(), iter.col(), in(iter.row(), iter.col()));
				++iter;
			}
			//memcpy(data(), in.data(), thetype.numelements() * sizeof(T));

		}

		//converts one type to another
		template<class T1,class instruct>
		_matrix(const _matrix<T1, instruct> &in) :
			MemChunkRef<T>(),  thetype(in.rows(), in.cols())
		{
			MemChunkRef<T>::newBlock(thetype.numelements());
			typename structure::iterator iter(rows(), cols());
			while (iter){
				put(iter.row(), iter.col(), T(in(iter.row(), iter.col())));
				++iter;
			}
		}

		template<class T1,class instruct>
		_matrix(const _matrixReference<T1, instruct> &in) :
			MemChunkRef<T>(*in._matrix_),  thetype(in.rows(), in.cols())
		{
			//MemChunkRef<T>::newBlock(thetype.numelements());
			//typename structure::iterator iter(rows(), cols());
			//while (iter){
			//	put(iter.row(), iter.col(), T(in(iter.row(), iter.col())));
			//	++iter;
			//}
			//data_=in.data_;
		}

		template<class exp1>
		_matrix(const _matrixExpr<exp1> &in)
		{
			resize(in.rows(rows()), in.cols(cols()));
			(*this)=in;
		}

		//fills matrix with some default value
		_matrix(int rows, int cols, T in):thetype(rows, cols){
			MemChunkRef<T>::newBlock(thetype.numelements());
			fill(in);
		}
		//fills matrix with some default value
		template<class T1>
		_matrix(int rows, int cols, T1 in): thetype(rows, cols){
			MemChunkRef<T>::newBlock(thetype.numelements());
			fill(in);
		}

		//fills matrix from a 2D datyapointer
		template<class T1>
		_matrix(int rows, int cols, T1 **in): thetype(rows, cols){
			MemChunkRef<T>::newBlock(thetype.numelements());
			for(int i=0;i<rows;++i){
				for(int j=0;j<cols;++j){
					put(i,j,in[i][j]);
				}
			}
		}

		//fills matrix from a 1D datyapointer
		//meant really for coping the data array
		//from one array to the next (useful for MPI passing)
		// three policies are inplace,
		//-matrix::newcopy--means
		//it will copy te data and you will need to delete it
		//-matrix::point--means the data will point to the
		// input data, and the matrix subclass 'memchunk' will
		// delete it when it is finished
		// -matrix::pointnodelete--means the matrix will point
		//  to the data and NOT delete it
		enum matCopyOps{newCopy, point, pointNoDelete};

		_matrix(int rows, int cols, T *in, matCopyOps mops=newcopy):
		  thetype(rows, cols)
		{
			if(mops==point){
				MemChunkRef<T>::createChunk(thetype.numelements(), in, duplicateData);
			}else if(mops==pointNoDelete){
				MemChunkRef<T>::createChunk(thetype.numelements(), in, neverDeleteData);
			}else{
				MemChunkRef<T>::newBlock(thetype.numelements());
				typename INstructure::iterator iter(rows, cols);
				int ct=0;
				while(iter){
					put(iter.row(), iter.col(), in[iter.col()*rows+iter.row()]);
					++ct;
					++iter;
				}
			}
		}

	//copy
		void copy(const _matrix &in)
		{
			this->resize(in.rows(), in.cols());
			if(!empty())
				std::memcpy(&data()[0], &in.data()[0], thetype.numelements() * sizeof(T));

		}
//aplies a 'function' across an entire vector
//The input class MUST have the operator(T, int) defined

		template<class Func>
		void apply( Func &f)
		{
			typename INstructure::iterator iter(rows(), cols());
			while(iter){
				put(iter.row(), iter.col(), f(iter.row(), iter.col()));
				++iter;
			}
		}

//aplies a 'function' across an a Range in the vector
//The input class MUST have the operator(T, int) defined

		template<class Func>
		void apply( Func &f, const Range &r, const  Range &r2)
		{
			for(int i=r.first(0); i<r.last(rows);++i){
				for(int j=r2.first(0); j<r2.last(cols);++j){
					put(i,j, f(i,j));
				}
			}
		}



		//******************************************/
		//   getting the _matrix identity parts

		MatrixType type()const {	return thetype.type();	}
		bool issquare()const{ return thetype.issquare();	}

		int cols()const{ return thetype.cols();	}
		int cols(int guessCols)const{ return thetype.cols();	}
		int columns(){	return thetype.cols();	}
		int rows()const{	return thetype.rows();	}
		int rows(int guessRows)const{	return thetype.rows();	}
		int numElements(){	return thetype.numEle();	}
		int numEle(){	return thetype.numEle();	}
		int numElements()const{	return thetype.numEle();	}
		int numEle()const{	return thetype.numEle();	}
		int position(int r, int c) const { return thetype.position(r,c);	}
		bool empty() const { return (thetype.rows()==0&&thetype.cols()==0)?true:false;	}


		//************************
		///		getting and assining bits of the _matrix

		T operator()(int i, int j) const{	return thetype.get(i,j,data_);	}
		T& operator()(int i, int j){	return thetype.get(i,j,data_);	   	}
		T operator()(int i) const{	return thetype.get(i,data_);	}
		T& operator()(int i){	return thetype.get(i,data_);	   	}

//Grabs a Column Vector in the Range
//Range operators (Range, int)
		Vector<T> operator()(const Range &r, int c) const
		{
#ifdef MatDeBug
			if(r.first(0) > rows() || r.last(rows()-1)>rows()-1) RowErr();
			if(c > cols() || c<0) ColErr();
#endif

			Vector<T> out(r.length(rows()));
			for(int i=0;i<r.length(rows());i++)
			{
				out(i)=get(r(i), c);
			}
			return out;
		}

//Grabs a ROW Vector in the Range
//Range operators (int, Range)
		Vector<T> operator()(int r, const Range &c) const
		{
#ifdef MatDeBug
			if(c.first(0) > cols() || c.last(cols()-1)>cols()-1) ColErr();
			if(r > rows() || r<0) RowErr();
#endif
			Vector<T> out(c.length(cols()));
			for(int i=0;i<c.length(cols());i++)
			{
				out( i)=get(r, c(i));
			}
			return out;
		}

//Sub Matrix Getting
//Range operators (Range, Range)
		_matrix operator()(const Range &r, const Range &c)
		{
#ifdef MatDeBug
			if(c.last(cols()-1)>=cols()-1) ColErr();
			if(r.last(rows()-1)>=rows()-1) RowErr();

#endif
			if(r.first(0)>=rows() || c.first(0)>=cols()) return _matrix(0,0,0);
			_matrix out(r.length(),c.length());
			typename structure::iterator iter(out.rows(), out.cols());
			while(iter){
				out(iter.row(), iter.col())=get(r(iter.row()), c(iter.col()));
				++iter;
			}
			return out;
		}

//Sub Matrix Getting (beginRow, beingColumn, endRow, endColumn)
		_matrix operator()(int br, int bc, int er, int ec)
		{
#ifdef MatDeBug
			if(br>=cols() || ec>=cols()) ColErr();
			if(br>=rows() || er>=rows()) RowErr();
#endif

			return operator()(Range(br, er), Range(bc, ec));
		}

//SETTING elements in the ranges....

//sets a COLUMN of elements in the Range
//The Vector length should be the length of the range...
//the vector should be ordered like [0...length(range)]
//the range simply determins the posision in the matrix..
		template<class newT>
		void put(const Range &r, int c, const Vector<newT> &in)
		{
#ifdef MatDeBug
			if(r.length()>in.length() || in.length()>rows()) RowErr();
			if(c>=cols() || c<0) ColErr();
#endif
			for(int i=0;i<r.length(rows());++i)
			{
				put(r(i),c, in(i));
			}
		}

//Vector Expression version
		template<class expr>
		void put(const Range &r, int c, const VExpr<expr> &in)
		{
#ifdef MatDeBug
			if(r.length(rows())>in.length(rows()) || in.length(rows())>rows()) RowErr();
			if(c>=cols() || c<0) ColErr();
#endif
			for(int i=0;i<r.length(rows());++i)
			{
				put(r(i),c, in(i));
			}
		}

//sets a ROW of elements in the Range
//The Vector length should be the length of the range...
//the vector should be ordered like [0...length(range)]
//the range simply determins the posision in the matrix..
		template<class newT>
		void put(int r, const Range &c, const Vector<newT> &in)
		{
#ifdef MatDeBug
			if(c.length()>in.length() || in.length()>cols()) ColErr();
			if(r>=rows() || r<0) RowErr();
#endif
			for(int i=0;i<c.length(cols());++i)
			{
				put(r,c(i), in(i));
			}
		}

//Vector Expression version
		template<class expr>
		void put(int r, const Range &c, const VExpr<expr> &in)
		{
#ifdef MatDeBug
			if(c.length(cols())>in.length(cols()) || in.length(cols())>cols()) ColErr();
			if(r>=rows() || r<0) RowErr();
#endif
			for(int i=0;i<c.length(cols());++i)
			{
				put(r,c(i), in(i));
			}
		}

//sets a TOTAL ROW in the matrix
//The Vector length should be the length of the cols()...
		template<class newT>
		void putRow(int r, const Vector<newT> &in)
		{
#ifdef MatDeBug
			if(r>=rows() || r<0) RowErr();
			if(in.length() != cols()) ColErr();
#endif
			for(int i=0;i<cols();++i)
			{
				put(r,i, in(i));
			}
		}

//sets a TOTAL ROW in the matrix
//The Vector length should be the length of the cols()...
		template<class newT, int N>
		void putRow(int r, const coord<newT, N> &in)
		{
#ifdef MatDeBug
			if(r>=rows() || r<0) RowErr();
			if(in.length() != cols()) ColErr();
#endif
			for(int i=0;i<cols();++i)
			{
				put(r,i, in(i));
			}
		}

//Vector Expresion version
		template<class expr>
		void putRow(int r, const VExpr<expr> &in)
		{
#ifdef MatDeBug
			if(r>=rows() || r<0) RowErr();
			if(in.length(cols()) != cols()) ColErr();
#endif
			for(int i=0;i<cols();++i)
			{
				put(r,i, in(i));
			}
		}

//sets a TOTAL COLUMN in the matrix
//The Vector length should be the length of the rows()...
		template<class newT>
		void putCol(int c, const Vector<newT> &in)
		{
#ifdef MatDeBug
			if(c>=cols() || c<0) ColErr();
			if(in.length() != rows()) RowErr();
#endif
			for(int i=0;i<rows();++i)
			{
				put(i,c, in(i));
				//cout<<i<<" "<<in(i)<<endl;
			}
		}

//sets a TOTAL COLUMN in the matrix
//The Vector length should be the length of the rows()...
		template<class newT,int N>
		void putCol(int c, const coord<newT,N> &in)
		{
#ifdef MatDeBug
			if(c>=cols() || c<0) ColErr();
			if(in.length() != rows()) RowErr();
#endif
			for(int i=0;i<rows();++i)
			{
				put(i,c, in(i));
				//cout<<i<<" "<<in(i)<<endl;
			}
		}

//Vector Expresion version
		template<class expr>
		void putCol(int c, const VExpr<expr> &in)
		{
#ifdef MatDeBug
			if(c>=cols() || c<0) ColErr();
			if(in.length(rows()) != rows()) RowErr();
#endif
			for(int i=0;i<rows();++i)
			{
				put(i,c, in(i));
			}
		}

//sets a TOTAL SUB MATRIX in the matrix
//THe first Range determins WHERE to put it in ROWS
//THe second Range determins WHERE to put it in COLUMNS
//The matrix  should be as big or smaller then the original matrix (this)
		template<class newT, class struc>
		void put(const Range &r, const Range &c, const  _matrix<newT, struc> &in)
		{
#ifdef MatDeBug
			if(c.length(cols())>=cols() || c.first(0)<0 || c.last(cols())>=cols()) ColErr();
			if(r.length(rows())>=rows() || r.first(0)<0 || r.last(rows())>=rows()) RowErr();
			if(in.rows()>rows() || in.cols() >cols()) AssErr();
#endif
			for(int i=0;i<r.length(rows());++i)
			{
				for(int j=0;j<c.length(cols());++j){
					put(r(i),c(j), in(i,j));
				}
			}
		}

//Matrix Expresion version
		template<class expr>
		void put(const Range &r, const Range &c, const  _matrixExpr<expr> &in)
		{
#ifdef MatDeBug
			if(c.length(cols())>=cols() || c.first(0)<0 || c.last(cols())>=cols()) ColErr();
			if(r.length(rows())>=rows() || r.first(0)<0 || r.last(rows())>=rows()) RowErr();
			if(in.rows(rows())>rows() || in.cols(cols()) >cols()) AssErr();
#endif
			for(int i=0;i<r.length(rows());++i)
			{
				for(int j=0;j<c.length(cols());++j){
					put(r(i),c(j), in(i,j));
				}
			}
		}

//the 'single' number version
		template<class newT>
		void put(const Range &r, const Range &c, const  newT &in)
		{
#ifdef MatDeBug
			if(c.length(cols())>=cols() || c.first(0)<0 || c.last(cols())>=cols()) ColErr();
			if(r.length(rows())>=rows() || r.first(0)<0 || r.last(rows())>=rows()) RowErr();
#endif
			for(int i=0;i<r.length(rows());++i)
			{
				for(int j=0;j<c.length(cols());++j){
					put(r(i),c(j), in);
				}
			}
		}

		T get(int i, int j) const { return thetype.get(i,j, data_);	}
		T& get(int i, int j){ return thetype.get(i,j, data_);	}


		Vector<T> getRow(int therow){
#ifdef MatDeBug
			if(therow<0 || therow>rows()) RowErr();
#endif
			Vector<T> out(cols());
			for(int i=0;i<cols();++i){
				out(i)=data_[position(therow, i)];
			}
			return out;
		}

		Vector<T> getCol(int thecol){
#ifdef MatDeBug
			if(thecol<0 || thecol>cols()) ColErr();
#endif
			Vector<T> out(rows());
			for(int i=0;i<rows();++i){
				out(i)=(*this)(i, thecol);
			}
			return out;
		}



		void put(int i, int j, T in)const{ return thetype.put(i,j, data_, in);	}

		template<class T1>
		void put(int i, int j, T1 in)const{ return thetype.put(i,j, data_, T(in));	}

//***************SPECIFIC TO TRIDIAGONALS **********************
//Diagonal pUtter
		void putD(int i, T in)const{ return thetype.putD(i, data_, in);	}

		template<class T1>
		void putD(int i, T1 in)const{ return thetype.putD(i, data_, T(in));	}

//Left Column pUtter
		void putL(int i, T in)const{ return thetype.putL(i, data_, in);	}

		template<class T1>
		void putL(int i, T1 in)const{ return thetype.putL(i, data_, T(in));	}

//Right Column pUtter
		void putR(int i, T in)const{ return thetype.putR(i, data_, in);	}

		template<class T1>
		void putR(int i, T1 in)const{ return thetype.putR(i, data_, T(in));	}

//sets The DIAGONAL Range from a Number
//THe first Range determins WHERE to put it in ROWS
		void putD(const Range &r, const  numtype &in)
		{
#ifdef MatDeBug
			if(r.length(rows())>=rows() || r.first(0)<0 || r.last(rows())>=rows()) RowErr();
#endif
			for(int i=0;i<r.length(rows());++i)	putD(r(i), in);
		}

//sets The RIGHT COLUMN Range from a Number
//THe first Range determins WHERE to put it in ROWS
		void putR(const Range &r, const  numtype &in)
		{
#ifdef MatDeBug
			if(r.length(rows())>=rows() || r.first(0)<0 || r.last(rows())>=rows()-1) RowErr();
#endif
			for(int i=0;i<r.length(rows());++i)	putR(r(i), in);

		}

//sets The LEFT COLUMN Range from a Number
//THe first Range determins WHERE to put it in ROWS
		void putL(const Range &r, const  numtype &in)
		{
#ifdef MatDeBug
			if(r.length(rows())>=rows() || r.first(0)<1 || r.last(rows())>=rows()) RowErr();
#endif
			for(int i=0;i<r.length(rows());++i)	putL(r(i), in);

		}

//sets The DIAGONAL Range from a Vector
//THe first Range determins WHERE to put it in ROWS
		template<class newT>
		void putD(const Range &r, const  Vector<newT> &in)
		{
#ifdef MatDeBug
			if(r.length(rows())>=rows() || r.first(0)<0 || r.last(rows())>=rows()) RowErr();
			if(in.size()>rows() || in.size() >cols()) AssErr();
#endif
			for(int i=0;i<r.length(rows());++i)	putD(r(i), in(i));
		}

//sets The RIGHT COLUMN Range from a Vector
//THe first Range determins WHERE to put it in ROWS
		template<class newT>
		void putR(const Range &r, const  Vector<newT> &in)
		{
#ifdef MatDeBug
			if(r.length(rows())>=rows() || r.first(0)<0 || r.last(rows())>=rows()) RowErr();
			if(in.size()>rows() || in.size() >cols()) AssErr();
#endif
			for(int i=0;i<r.length(rows());++i)	putR(r(i), in(i));

		}

//sets The LEFT COLUMN Range from a Vector
//THe first Range determins WHERE to put it in ROWS
		template<class newT>
		void putL(const Range &r, const  Vector<newT> &in)
		{
#ifdef MatDeBug
			if(r.length(rows())>=rows() || r.first(0)<0 || r.last(rows())>=rows()) RowErr();
			if(in.size()>rows() || in.size() >cols()) AssErr();
#endif
			for(int i=0;i<r.length(rows());++i)	putL(r(i), in(i));

		}
/**************************END TRIDIAG SPECIFIC BITS *********************/


		T *data(){ return data_;	}
		const T *data()const{ return data_;	}

		template<class T1>
		void fill(T1 in){
			typename structure::iterator iter(rows(), cols());
			while(iter){
				put(iter.row(), iter.col(),T(in));
				++iter;
			}
		}

		void identity(){
			typename structure::iterator iter(rows(), cols());
			while(iter){
				if(iter.row()==iter.col()){	put(iter.row(), iter.col(),OneType<T>::one());}
				else{ put(iter.row(), iter.col(),ZeroType<T>::zero()); }
				++iter;
			}
		}

		void identity(int n){
			resize(n,n);
			typename structure::iterator iter(rows(), cols());
			while(iter){
				if(iter.row()==iter.col()){	put(iter.row(), iter.col(),OneType<T>::one());}
				else{ put(iter.row(), iter.col(),ZeroType<T>::zero()); }
				++iter;
			}
		}

		void identity(int n, int c){
			resize(r,c);
			typename structure::iterator iter(rows(), cols());
			while(iter){
				if(iter.row()==iter.col()){	put(iter.row(), iter.col(),OneType<T>::one());}
				else{ put(iter.row(), iter.col(),ZeroType<T>::zero()); }
				++iter;
			}
		}

		//reshape the entity...gives new cols and row defs
		// the rows*cols must be the same for the new and old ones
		void reshape(int nr, int nc)
		{
			if(nr*nc != rows()*cols()){
				BLEXCEPTION(" number of elements MUST be the same for the new rows and cols")
			}
			thetype.resize(nr,nc);
		}

		void resize(int r, int c){
			if(rows()!=r || cols()!=c){
				thetype.resize(r,c);
				MemChunkRef<T>::newBlock(thetype.numEle());
			}
	    }

		void resizeAndPreserve(int r, int c){
			if(rows()!=r || cols()!=c){
				_matrix tmm((*this));
				resize(r,c);
				int i,j, mr=(rows()>tmm.rows())?tmm.rows():rows(), mc=(cols()>tmm.cols())?tmm.cols():cols();
				for(i=0;i<mr;++i)
					for(j=0;j<mc;++j) put(i,j,tmm(i,j));

				//MemChunkRef<T>::changeBlock(tmm, thetype.numEle());
			}
	    }

	    template<class newT>
		void resizeAndPreserve(int r, int c, const newT &in){
			if(rows()!=r || cols()!=c){
				_matrix tmm((*this));
				resize(r,c, in);
				int i,j, mr=(rows()>tmm.rows())?tmm.rows():rows(), mc=(cols()>tmm.cols())?tmm.cols():cols();
				for(i=0;i<mr;++i)
					for(j=0;j<mc;++j) put(i,j,tmm(i,j));

				//MemChunkRef<T>::changeBlock(tmm, thetype.numEle());
			}
	    }

	    void resize(int r, int c, T in){
			if(rows()!=r || cols()!=c){
				thetype.resize(r,c);
				MemChunkRef<T>::newBlock(thetype.numEle());
			}
			fill(in);
	    }

		void resizeNew(int r, int c){
			thetype.resize(r,c);
			MemChunkRef<T>::newBlock(thetype.numEle());
	    }

		//used in the template expansions of operations
		// passes a pointer around rather then the big
		// _matrix parts...speeds up the works a bit
		_matrixReference<numtype, structure> getRef() const{
    		return _matrixReference<numtype, structure>(*this);
		}

		//************assignment operations  See "matassign.h"
		void operator=(const _matrix &mat);

		template<class theexpr>
    	void operator=(const _matrixExpr<theexpr> &INexpr);

		template<class nt2, class st2>
		void operator=(const _matrix<nt2, st2> &mat);

		template<class nt1>
		void operator=(const nt1 &num){		fill(num);		}
/*
		template<class nt1, class st1>
		bool operator==(const _matrix<nt1, st1> &rhs)
		{
			if(rows()!=rhs.rows()) return false;
			if(cols()!=rhs.cols()) return false;

			iterator i(rows(), cols());
			while(i){
				if(this->get(i.row(), i.col()) != rhs(i.row(),i.col()) )	return false;
				++i;
			}
			return true;
		}

		template<class expr>
		bool operator==(const _matrixExpr<expr> &rhs)
		{
			if(rows()!=rhs.rows(0)) return false;
			if(cols()!=rhs.cols(0)) return false;

			iterator i(rows(), cols());
			while(i){
				if(this->get(i.row(), i.col()) != rhs(i.row(),i.col()) )	return false;
				++i;
			}
			return true;
		}

		template<class nt1, class st1>
		bool operator!=(const _matrix<nt1, st1> &rhs)
		{
			if(rows()!=rhs.rows()) return true;
			if(cols()!=rhs.cols()) return true;

			iterator i(rows(), cols());
			while(i){
				if(this->get(i.row(), i.col()) != rhs(i.row(),i.col()))	return true;
				++i;
			}
			return false;
		}

		template<class expr>
		bool operator!=(const _matrixExpr<expr> &rhs)
		{
			if(rows()!=rhs.rows(0)) return true;
			if(cols()!=rhs.cols(0)) return true;

			iterator i(rows(), cols());
			while(i){
				if(this->get(i.row(), i.col()) != rhs(i.row(),i.col()))	return true;
				++i;
			}
			return false;
		}
*/
	//*************************row, column extraction...returns a vector of type T

		Vector<T> row(int therow) const
		{
#ifdef MatDeBug
			if(therow<0 || therow>rows()) RowErr();
#endif
			Vector<T> out(cols());
			for(int i=0;i<cols();++i){
				out(i)=(*this)(therow, i);
			}
			return out;
		}

		Vector<T> col(int thecol) const
		{
#ifdef MatDeBug
			if(thecol<0 || thecol>cols()) ColErr();
#endif
			Vector<T> out(rows());
			for(int i=0;i<rows();++i){
				out(i)=(*this)(i, thecol);
			}
			return out;
		}

		T *row_p(int therow){
#ifdef MatDeBug
			if(therow<0 || therow>rows()) RowErr();
#endif
			T *out;
			if(type()==Mfull || type()==Mgeneral){
				out=new T[cols()];
				for(int i=0;i<rows();++i){
					out[i]=&data_[position(i, thecol)];
				}
				return out;
			}else{
				BLEXCEPTION(" Can only use this function for 'FullMatrix' types")
			}
			return out;
		}

		T *col_p(int therow){
#ifdef MatDeBug
			if(therow<0 || therow>cols()) ColErr();
#endif
			T *out;
			if(type()==Mfull || type()==Mgeneral){
				 out=new T[rows()];
				for(int i=0;i<rows();++i){
					*out[i]=&data_[position(i, therow)];
				}
				return out;
			}else{
				 BLEXCEPTION(" Can only use this function for 'FullMatrix' types")
			}
			return out;
		}

		//************************************
		//	self operations...see "matassign.h" for these functions..
		template<class nt1, class st1> inline void operator+=(const _matrix<nt1,st1 > &);
		template<class nt1, class st1> inline void operator-=(const _matrix<nt1,st1 > &);
		template<class nt1, class st1> inline void operator*=(const _matrix<nt1,st1 > &);
		//template<class nt1, class st1> inline matrix_& operator*=( _matrix<nt1,st1 > &);

		template<class nt1, class st1>
		inline void operator/=(const _matrix<nt1,st1 > &);

		template<class expr> inline void operator+=(const _matrixExpr<expr> & );
		template<class expr> inline void operator-=(const _matrixExpr<expr> & );
		//template<class expr> inline matrix_& operator/=(const _matrixExpr<expr> & );
		template<class expr> inline void operator*=(const _matrixExpr<expr> & );
		template<class nt2> inline void operator *= (const _matrixExpr<_matrixExprConst<nt2> > &);

		template<class nt2> inline void operator /= (const _matrixExpr<_matrixExprConst<nt2> > &);


		inline void operator+=(const complex &);
		inline void operator-=(const complex &);
		inline void operator*=(const complex &);
		inline void operator/=(const complex &);

		inline void operator+=(const double &);
		inline void operator-=(const double &);
		inline void operator*=(const double &);
		inline void operator/=(const double &);

		inline void operator+=(const float &);
		inline void operator-=(const float &);
		inline void operator*=(const float &);
		inline void operator/=(const float &);

		inline void operator+=(const int &);
		inline void operator-=(const int &);
		inline void operator*=(const int &);
		inline void operator/=(const int &);

		inline void operator+=(const short &);
		inline void operator-=(const short & );
		inline void operator*=(const short & );
		inline void operator/=(const short &);

		inline void operator+=(const long &);
		inline void operator-=(const long &);
		inline void operator*=(const long &);
		inline void operator/=(const long &);

		inline void operator+=(const char &);
		inline void operator-=(const char &);
		inline void operator*=(const char &);
		inline void operator/=(const char &);

		inline void operator+=(const bool&);
		inline void operator-=(const bool &);
		inline void operator*=(const bool &);
		inline void operator/=(const bool & );

//internla multiply for the same type..should be faster
//This does This=rhs*this
//NOTE:: THE NMR PROPOGATOR CASE!!! NOT THE NORMAL ONE!!!
		_matrix<complex,FullMatrix> &
		  operator*=(const _matrix<Complex<double>, FullMatrix>  &rhs);

		_matrix<Complex<float>,FullMatrix> &
		  operator*=(const _matrix<Complex<float>, FullMatrix>  &rhs);

		_matrix<double,FullMatrix> &
		  operator*=(const _matrix<double, FullMatrix> &rhs);

/**** 'INternal Adjoint and transpose ***/
//this solves the A=adjoint(A) memory problem
//
		void adjoint();
		void transpose();

/**** 'Internal Multiplications' ****/
// things like B=A*B*adjoint(A), B=A*B;
// where we can avoid memory copy
// note-->B=B*A is already optimized

		//this is B=A*B*adjoint(A)-->similarity tranformations
		template<class nt1, class st1>
		void simTrans(const _matrix<nt1, st1> &A);
		template<class nt1, class st1>
		void prop(const _matrix<nt1, st1> &A); //same as above
		template<class nt1, class st1>
		void adjprop(const _matrix<nt1, st1> &A); //B=adjoint(A)*B*A

		template<class Ctype>
		void simTrans(const _matrix<Complex<Ctype>, FullMatrix> &A);
		void simTrans(const _matrix<double, FullMatrix> &A);
		void simTrans(const _matrix<float, FullMatrix> &A);

		template<class nt1>
		void prop(const _matrix<nt1, FullMatrix> &A); //same as above
		template<class nt1>
		void adjprop(const _matrix<nt1, FullMatrix> &A); //B=adjoint(A)*B*A

		//multiply and assign B=A*B
		template<class nt1, class st1>
		void mulAss(const _matrix<nt1, st1> &inmat);

/*** Interanl Diag routines ***/
		template<class nt1>
		void diag( _matrix<nt1, DiagonalMatrix> &dia,  _matrix<nt1, FullMatrix> &U, const bool destroy=false);

		template<class nt1>
		void diagonalize( _matrix<nt1, DiagonalMatrix> &dia,  _matrix<nt1, FullMatrix> &U, const bool destroy=false);

		//returns ONLY the eignvalues inside dis (does not calc evects)
		template<class nt1>
		void eigenvalues(_matrix<nt1, DiagonalMatrix> &dia);

//**** LU decomp functions

	//these 3 functions technically should not be here, but it is the only way to
	//get the header files to be accepted
		template<class T1>
		int LUsolve_(Vector<T1> &b);
		int LUdecomp_( Vector<int> &rowperm, float &d);

	//NON-destructive LU decomp into L and U into the same matrix
		int LUdecomp_( Vector<int> &rowperm, float &d, _matrix<T, FullMatrix> &LU);
	//NON-destructive LU decomp into L and U
		int LUdecomp_( Vector<int> &rowperm, float &d, _matrix<T, FullMatrix> &L, _matrix<T, FullMatrix> &U);

		template<class T1>
		int LUbackSub_( Vector<int> &rowperm, Vector<T1> &b);

	//spcial case for Vector<coord<> >
		template<int N>
		int LUbackSub_( Vector<int> &rowperm, Vector<coord<T,N> > &b);

	//returns -1 if signluar matrix
		Vector<int> LUdecomp();
		int LUdecomp(Vector<int> &rowperm);
		int LUdecomp(Vector<int> &rowperm, float &d);
		int LU(); //destructive single shot LU

	//'Normal' LUs
	//NON destructive LUs ...(slower, but potential more usful)
		int LU(_matrix<T, FullMatrix> &LU);
		int LU(_matrix<T, FullMatrix> &L,_matrix<T, FullMatrix> &U);

		template<class T1>
		void LUbackSub(Vector<int> &rowperm, Vector<T1> &b);

	//spcial case for Vector<coord<> >
		template<int N>
		void LUbackSub(Vector<int> &rowperm, Vector<coord<T, N> > &b);

		template<class T1>
		void LUsolve(Vector<T1> &b);

///matrix inverse...uses LUdecomp
		_matrix<T, invers_structure> inv() const;

//the linear solve function (returns the solutions)
		Vector<T> solve(const Vector<T> &r);

//QR and Gram Schmit orthogonalization
//see 'matqr.h' for these

	//returns the Q and R matrix from the QR decomposistion
	// this=QR...Q is orthonormal, and R is upper-trianular (with a diagonal)
	//NOT destructive to this matrix
		void QR(_matrix<T, FullMatrix> &Q, _matrix<T, FullMatrix> &R);

	//returns the W and R matrix from the QR decomposistion
	//The R amtrix is in (THIS) i.e. this is destructive
	// this=QR...Q is orthonormal, and R is upper-trianular (with a diagonal)
	//destructive to this matrix
	// W--> is the 'packed' version of 'Q' (use UnPackQR(W) to unpack it)
	// R-->is theupper tri matrix (THIS)
		void PackedQR(_matrix<T, FullMatrix> &W);

	//unpacks the 'packed' W from the 'PackedQR' function above
	//Q--> becomes a full orthonormal matrix
		void UnPackQR(_matrix<T, FullMatrix> &W, _matrix<T, FullMatrix> &Q);
		template<class T2>
		void chopper(T2 cutoff);
		template<class T2>
		void chop(T2 cutoff);
		void ReNorm(double cutoff);

};


//**************************************
//	_matrix EXPRESION
//		the main widget where all _matrix aoperation evenutally find their way to
//		it is the drive that does the recursive expansion of the expressions

template<class INexpr>
class _matrixExpr {
	private:
		INexpr expr_;

	public:
		const void MulErr() const {
			BLEXCEPTION(" size mismatch  m1.rows()==m2.cols()")
		}
		typedef INexpr expr;
		typedef typename expr::numtype numtype;
		typedef typename expr::structure structure;
		typedef typename expr::invers_structure invers_structure;


		_matrixExpr(const expr &in): expr_(in) {}
		numtype operator()(int i, int j)const {
			return expr_(i,j);
		}

		int cols(int guessCols)const{ return expr_.cols(guessCols);	}
		int rows(int guessRows)const{ return expr_.rows(guessRows);	}

};

//**************************************
//	_matrix Doubel Operation EXPRESION
//		holds the Expresion types for double operations (i.e. A+B)

template<class INexpr1, class INexpr2, class INop>
class _matrixExprBinOp{
	private:
		INexpr1 expr1_;
		INexpr2 expr2_;
		void lenerr()const{
			 BLEXCEPTION(" _matrix sizes MUST agree")
		}
	public:
		typedef INexpr1 expr1;
		typedef INexpr2 expr2;
		typedef typename expr1::numtype numtype1;
		typedef typename expr2::numtype numtype2;
		typedef OutType(numtype1, numtype2) numtype;
		typedef typename MatOutType(expr1, expr2) structure;
		typedef typename MatOutTypeS(typename INexpr1::invers_structure, typename INexpr2::invers_structure) invers_structure;

		typedef INop op;

		//_matrixExprBinOp(expr1 A, expr2 B): expr1_(A), expr2(B) {}
		_matrixExprBinOp(const expr1 &A, const expr2 &B): expr1_(A), expr2_(B) {
			}

		numtype operator()(int i, int j) const{
			return op::apply(expr1_(i,j), expr2_(i,j));
		}

		int cols(int guessCols)const{
			//if(expr1_.cols(guessCols)==-1) return expr2_.cols(guessCols);
			//if(expr2_.cols(guessCols)==-1) return expr1_.cols(guessCols);
			//if(expr1_.cols(guessCols)!=expr2_.cols(guessCols)){
			//	std::cout<<"BINOPcc: "<<expr1_.cols(guessCols)<<" "<<expr2_.cols(guessCols)<<std::endl;
			//	lenerr();
			//}
			if(expr1_.cols(guessCols)>expr2_.cols(guessCols)) return expr1_.cols(guessCols);
			return expr2_.cols(guessCols);
		}

		int rows(int guessRows)const{
			//if(expr1_.rows(guessRows)==-1) return expr2_.rows(guessRows);
			//if(expr2_.rows(guessRows)==-1) return expr1_.rows(guessRows);
			//if(expr1_.rows(guessRows)!=expr2_.rows(guessRows)) {
			//std::cout<<"BINOPrr: "<<expr1_.rows(guessRows)<<" "<<expr2_.rows(guessRows)<<std::endl;
			//	lenerr();
			//}
			if( expr1_.rows(guessRows)>expr2_.rows(guessRows)) return expr1_.rows(guessRows);
			return expr2_.rows(guessRows);
		}
};

//**************************************
//	_matrix Doubel Operation EXPRESION FOR !!MULTIPLICATION!!
//		holds the Expresion types for double operations (i.e. A*B)
//		need a special one becuase the (i,j) element is acctually
//		a sum -->  c(i,j)= Sum_k [ a(i,k)*b(k,j) ]
//		and so the number of cols in a must be the number of rows in b

template<class INexpr1, class INexpr2>
class _matrixExprBinOpMUL{
	private:
		INexpr1 expr1_;
		INexpr2 expr2_;

		const void lenerr()const{
			BLEXCEPTION(" m1 cols() ==Must be== m2 rows()")
		}
	public:

		typedef INexpr1 expr1;
		typedef INexpr2 expr2;
		typedef typename expr1::numtype numtype1;
		typedef typename expr2::numtype numtype2;
		typedef SumType(OutType(numtype1, numtype2)) numtype;
		typedef typename MatOutType(expr1, expr2) structure;
		typedef typename MatOutTypeS(typename expr1::invers_structure, typename expr2::invers_structure) invers_structure;

		_matrix<numtype, structure> dd;

		_matrixExprBinOpMUL(const expr1 &A, const expr2 &B): expr1_(A), expr2_(B),
			dd(rows(expr1_.rows(1)), cols(expr2_.cols(1))){

			typename structure::iterator iter1(dd.rows(), dd.cols());
			while(iter1){
				dd.put(iter1.row(), iter1.col(), mul(iter1.row(), iter1.col()));
				++iter1;
			}
		}

		numtype operator()(int i, int j) const{
			return dd(i,j);
		}

		numtype mul(int i, int j) const{
			numtype tmp=ZeroType<numtype>::zero();
			for(int k=0;k<cols(0);++k)	 tmp+=expr1_(i,k) * expr2_(k,j);

			return tmp;
		}


		int cols(int guessCols)const{
#ifdef MatDeBug
			if(expr1_.cols(guessCols)!=expr2_.cols(guessCols)) lenerr();
#endif
			return expr1_.cols(guessCols);
		}

		int rows(int guessRows)const{
#ifdef MatDeBug
			if(expr1_.rows(guessRows)!=expr2_.rows(guessRows)) lenerr();
#endif
			return expr1_.rows(guessRows);
		}
};

//**************************************
//	_matrix Single Operation EXPRESION
//		holds the Expresion types for double operations (i.e. abs(A))

template<class INexpr, class INop>
class _matrixExprUnrOp{
	private:
		INexpr expr_;

	public:
		typedef INexpr expr;
		typedef typename INop::numtype numtype;
		typedef INop op;
		typedef typename expr::structure structure;
		typedef typename expr::invers_structure invers_structure;

		//_matrixExprUnrOp(expr A): expr_(A) {}
		_matrixExprUnrOp(const expr &A): expr_(A) {}

		numtype operator()(int i, int j) const{
			return op::apply(expr_(i,j));
		}

		int cols(int guessCols)const{	return expr_.cols(guessCols);	}
		int rows(int guessRows)const{	return expr_.rows(guessRows);	}
};


//**************************************
//	_matrix Constant Operation EXPRESION
//		holds the Expresion types for double operations (i.e. A+3)
//		basically converts the '3' into somthing the _matrixExpr can handle

template<class INtype>
class _matrixExprConst{
	private:
		INtype expr_;

	public:
		typedef INtype numtype;
		typedef INtype structure;
		typedef INtype invers_structure;

		//_matrixExprConst(INtype A): expr_(A) {}
		_matrixExprConst(const INtype &A): expr_(A) {}

		numtype apply(int i, int j)const {
			return expr_;
		}

		numtype operator()(int i, int j) const{
			return expr_;
		}

		int cols(int guessCols) const {	return guessCols;	}
		int rows(int guessRows) const {	return guessRows;	}
};


//**************************************
//	_matrix Reference
//		takes care of multiple _matrix references...instead of passing the
//		entire _matrix data, pass the pointer to it instead


template<class INnumtype, class INstructure>
class _matrixReference {
	private:
	    _matrixReference() { } //kill the null constructor from use

  public:
		typedef INnumtype numtype;
		typedef INstructure structure;
		typedef typename INstructure::invers_structure invers_structure;
		const _matrix<INnumtype, INstructure>* _matrix_;

		_matrixReference(const _matrix<numtype, INstructure>& m): _matrix_(&m){ }

		numtype operator()(int i, int j) const{ return (*_matrix_)(i,j); }
		numtype get(int i, int j) const{ return (*_matrix_)(i,j); }

		int rows(int GuessR) const  { return _matrix_->rows(); }
		int cols(int GuessC) const  { return _matrix_->cols(); }
		int rows() const  { return _matrix_->rows(); }
		int cols() const  { return _matrix_->cols(); }

};

//********************************************************************************
//*****************OPERATIONS!! OPERATIONS!! OPERATIONS!!*************************
//********************************************************************************
//
//  For these, to save a bit of typing some Marco's will be used, but in their definitions the
// following nameing convention is used
//
//		nt1-->the number type for the first _matrix (float, bool, int, complex, etc)
//		nt2--> "    "      "   "  "   "  second "
//		st1-->the structure type for the first _matrix (Full, Diagonal, etc)
//		st2--> "       "      "   "    " second   "
//		m1 --> the 'lefthand side' _matrix
//		m2 --> the 'righthand side' _matrix
//		ex1 --> left _matrix expresion
//		ex2 --> right _matrix expr

template<class nt1, class st1, class nt2, class st2>
bool operator==(const _matrix<nt1, st1> &lhs, const _matrix<nt2, st2> &rhs)
{
	if(lhs.rows()!=rhs.rows()) return false;
	if(lhs.cols()!=rhs.cols()) return false;

	typename _matrix<nt1, st1>::iterator i(lhs.rows(), lhs.cols());
	while(i){
		if(lhs(i.row(), i.col()) != rhs(i.row(),i.col()) )	return false;
		++i;
	}
	return true;
}

template<class nt2, class st2, class expr>
bool operator==(const _matrixExpr<expr> &lhs, const _matrix<nt2, st2> &rhs)
{
	if(lhs.rows(0)!=rhs.rows()) return false;
	if(lhs.cols(0)!=rhs.cols()) return false;

	typename _matrixExpr<expr>::structure::iterator i(rows(0), cols(0));
	while(i){
		if(lhs(i.row(), i.col()) != rhs(i.row(),i.col()) )	return false;
		++i;
	}
	return true;
}

template<class nt1, class st1, class expr>
bool operator==(const _matrix<nt1, st1> &lhs, const _matrixExpr<expr>  &rhs)
{
	if(lhs.rows()!=rhs.rows(0)) return false;
	if(lhs.cols()!=rhs.cols(0)) return false;

	typename _matrix<nt1, st1>::iterator i(lhs.rows(), lhs.cols());
	while(i){
		if(lhs(i.row(), i.col()) != rhs(i.row(),i.col()) )	return false;
		++i;
	}
	return true;
}

template<class expr1, class expr2>
bool operator==(const _matrixExpr<expr1> &lhs, const _matrixExpr<expr2> &rhs)
{
	if(lhs.rows(0)!=rhs.rows(0)) return false;
	if(lhs.cols(0)!=rhs.cols(0)) return false;

	typename _matrixExpr<expr1>::structure::iterator i(lhs.rows(0), lhs.cols(0));
	while(i){
		if(lhs(i.row(), i.col()) != rhs(i.row(),i.col()) )	return false;
		++i;
	}
	return true;
}

template<class nt1, class st1, class nt2, class st2>
bool operator!=(const _matrix<nt1, st1> &lhs, const _matrix<nt2, st2> &rhs)
{
	if(lhs.rows()!=rhs.rows()) return true;
	if(lhs.cols()!=rhs.cols()) return true;

	typename _matrix<nt1, st1>::iterator i(lhs.rows(), lhs.cols());
	while(i){
		if(lhs(i.row(), i.col()) != rhs(i.row(),i.col()) )	return true;
		++i;
	}
	return false;
}

template<class nt2, class st2, class expr>
bool operator!=(const _matrixExpr<expr> &lhs, const _matrix<nt2, st2> &rhs)
{
	if(lhs.rows()!=rhs.rows()) return true;
	if(lhs.cols()!=rhs.cols()) return true;

	typename _matrixExpr<expr>::structure::iterator i(lhs.rows(0), lhs.cols(0));
	while(i){
		if(lhs(i.row(), i.col()) != rhs(i.row(),i.col()) )	return true;
		++i;
	}
	return false;
}

template<class nt1, class st1, class expr>
bool operator!=(const _matrix<nt1, st1> &lhs, const _matrixExpr<expr>  &rhs)
{
	if(lhs.rows()!=rhs.rows(0)) return true;
	if(lhs.cols()!=rhs.cols(0)) return true;

	typename _matrix<nt1, st1>::iterator i(lhs.rows(), lhs.cols());
	while(i){
		if(lhs(i.row(), i.col()) != rhs(i.row(),i.col()) )	return true;
		++i;
	}
	return false;
}

template<class expr1, class expr2>
bool operator!=(const _matrixExpr<expr1> &lhs, const _matrixExpr<expr2> &rhs)
{
	if(lhs.rows(0)!=rhs.rows(0)) return true;
	if(lhs.cols(0)!=rhs.cols(0)) return true;

	typename _matrixExpr<expr1>::structure::iterator i(lhs.rows(0), lhs.cols(0));
	while(i){
		if(lhs(i.row(), i.col()) != rhs(i.row(),i.col()) )	return true;
		++i;
	}
	return false;
}

//*******************************
//	SELF operations (*=, -=, etc)
//*******************************

//******Mat OP NUM
#define MatFromConst(op, TYPE)	\
	template<class T, class  INstructure>\
	inline void _matrix<T, INstructure>::operator op (const TYPE &a){	\
		(*this) op _matrixExpr<_matrixExprConst<T> >(_matrixExprConst<T>(a));		\
	} \

MatFromConst(+=, complex)
MatFromConst(-=, complex)
MatFromConst(/=, complex)
MatFromConst(*=, complex)

MatFromConst(+=, double)
MatFromConst(-=, double)
MatFromConst(/=, double)
MatFromConst(*=, double)

MatFromConst(+=, float)
MatFromConst(-=, float)
MatFromConst(/=, float)
MatFromConst(*=, float)

MatFromConst(+=, int)
MatFromConst(-=, int)
MatFromConst(/=, int)
MatFromConst(*=, int)

MatFromConst(+=, short)
MatFromConst(-=, short)
MatFromConst(/=, short)
MatFromConst(*=, short)

MatFromConst(+=, long)
MatFromConst(-=, long)
MatFromConst(/=, long)
MatFromConst(*=, long)

MatFromConst(+=, char)
MatFromConst(-=, char)
MatFromConst(/=, char)
MatFromConst(*=, char)

MatFromConst(+=, bool)
MatFromConst(-=, bool)
MatFromConst(/=, bool)
MatFromConst(*=, bool)

template<class T, class  INstructure>
 template<class nt2>
inline void _matrix<T, INstructure>::operator *= (const _matrixExpr<_matrixExprConst<nt2> > &m2){
	typename structure::iterator iter(rows(), cols());
	while (iter){
		put(iter.row(), iter.col(), this->get(iter.row(), iter.col())*m2(iter.row(), iter.col()));
		++iter;
	}
}

template<class T, class  INstructure>
 template<class nt2>
inline void _matrix<T, INstructure>::operator /= (const _matrixExpr<_matrixExprConst<nt2> > &m2){
	typename structure::iterator iter(rows(), cols());
	while (iter){
		put(iter.row(), iter.col(), this->get(iter.row(), iter.col())/m2(iter.row(), iter.col()));
		++iter;
	}

}


//***********Matrix op Matrix


//MatFromMat(+=)  See "matassign.h"
//MatFromMat(-=)  See "matassign.h"
// ***
// ******SPECIAL FOR operator '_matrix*=_matrix'
// ***
template<class T, class  INstructure>
 template<class nt2, class st2>
inline void _matrix<T, INstructure>::operator *= ( const _matrix<nt2, st2> &m2){

#ifdef MatDeBug
	if(rows() != m2.cols()){
		BLEXCEPTION(" m1 rows() ==Must be== m2 cols()")
	}
#endif
	(*this)=(*this) * m2;
}

//****************************************************************************
//		See "matassign.h" allow for a _matrix<nt1, st1> from _matrixExpr<ex1>
//		 for use with unaray ops i.e (+=, -=)
//****************************************************************************

//***
//******SPECIAL FOR operator '_matrix*=_matrixExpr'
//***
template<class T, class INstructure> template<class ex2>
inline void _matrix<T, INstructure>::operator *= (const _matrixExpr<ex2> &m2){
#ifdef MatDeBug
	if(rows() != m2.cols(cols())){
		BLEXCEPTION(" m1 rows() ==Must be== m2 cols()")
	}
#endif
	(*this)=(*this)*m2;
}


//*******************************
//		+, -
//*******************************

//**the order of the operators given in the macro

// _matrix OP _matrix
// Matrix OP _matrixExpr
// MatriExpr OP _matrix
// MatriExpr OP _matrixExpr


#define MatMakeAS(OP, OPNAME) \
template<class nt1, class st1, class nt2, class st2> \
inline	\
_matrixExpr< _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2>, \
	OPNAME<nt1, nt2> > >	\
operator OP (const _matrix<nt1, st1>& m1, const _matrix<nt2, st2>& m2){	\
    typedef _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixReference<nt2, st2>, OPNAME<nt1, nt2> > expr;	\
    return _matrixExpr<expr>(expr(m1.getRef(),m2.getRef()));	\
}	\
	\
	\
template<class nt1, class st1, class ex2>	\
inline	\
_matrixExpr< _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixExpr<ex2>, \
	OPNAME<nt1,typename _matrixExpr<ex2>::numtype> > >	 \
operator OP (const _matrix<nt1, st1>& m1, const _matrixExpr<ex2>& m2){	\
    typedef _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixExpr<ex2>, OPNAME<nt1, typename _matrixExpr<ex2>::numtype> > expr;	\
    return _matrixExpr<expr>(expr(m1.getRef(),m2));	\
}	\
	\
	\
template<class ex1, class nt2, class st2>	\
inline	\
_matrixExpr< _matrixExprBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2>, \
	OPNAME<typename _matrixExpr<ex1>::numtype,nt2> > >	 \
operator OP (const _matrixExpr<ex1>& m1, const _matrix<nt2,st2>& m2){	\
    typedef _matrixExprBinOp< _matrixExpr<ex1>, _matrixReference<nt2, st2>, OPNAME<typename _matrixExpr<ex1>::numtype, nt2> > expr;	\
    return _matrixExpr<expr>(expr(m1,m2.getRef()));	\
}																								\
																								\
																								\
template<class ex1, class ex2>																	\
inline																							\
_matrixExpr< _matrixExprBinOp< _matrixExpr<ex1>, _matrixExpr<ex2>,								\
	OPNAME<typename _matrixExpr<ex1>::numtype,typename _matrixExpr<ex2>::numtype> > >			\
operator OP (const _matrixExpr<ex1>& m1, const _matrixExpr<ex2>& m2){							\
    typedef _matrixExprBinOp< _matrixExpr<ex1>, _matrixExpr<ex2>, 								\
    	OPNAME<typename _matrixExpr<ex1>::numtype, typename _matrixExpr<ex2>::numtype> > expr;	\
    return _matrixExpr<expr>(expr(m1,m2));									\
}																			\
																			\
																			\

MatMakeAS(+, ApAdd)
MatMakeAS(-, ApSubtract)


//the below '/' operators are correct (but the above macro redefines them)
//these are expressed BUT NOT VALID matrix divisions/inversions...need to get the LU decomp working first

//*******************************
//		*, /
//*******************************

//**the order of the operators given in the macro
//**NOTE:: the Matrix*Matrix types MUST be handled sparately
//			these handle the num*/Matrix types
//
//**NOTE: divide '/' uses ApMultiply and just inverts the number first

// _matrixExpr OP num
// num OP _matrixExpr

//**
//*********Division**
//**

// matrixExpr/num
template<class ex1, class nt2>
inline
_matrixExpr< _matrixExprBinOp< _matrixExpr<ex1>, _matrixExprConst<FloatType(nt2)>,
	ApMultiply<typename _matrixExpr<ex1>::numtype, FloatType(nt2)> > >
operator / (const _matrixExpr<ex1>& m1, nt2 m2){
    typedef _matrixExprBinOp< _matrixExpr<ex1>, _matrixExprConst<nt2>,
    	ApMultiply<typename _matrixExpr<ex1>::numtype, FloatType(nt2)> > expr;
    return _matrixExpr<expr>(expr(m1,FloatType(nt2)(1./m2)));
}

// matrix/num
template<class nt1, class st1, class nt2>
inline
_matrixExpr< _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixExprConst<FloatType(nt2)>,
	ApMultiply<nt1, FloatType(nt2)> > >
operator / (const _matrix<nt1, st1>& m1, nt2 m2){
    typedef _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixExprConst<FloatType(nt2)>,
    	ApMultiply<nt1, FloatType(nt2)> > expr;
    return _matrixExpr<expr>(expr(m1.getRef(),FloatType(nt2)(1./m2)));
}




/*****************************
	Need a special class for division and inversion things
	1) inversion of matricies (i.e. 1/m)
	2) Division of matrices by other matrices...uses Inverse
******************************/

//for this class the entity on the RIGHT SIDE is the one that needs the invers
template<class m1, class m2>
class divMatMatBinOp{
	private:
		m1 expr1_;
		m2 expr2_;


	public:

		typedef m1 expr1;
		typedef m2 expr2;
		typedef typename expr1::numtype numtype1;
		typedef typename expr2::numtype numtype2;
		typedef SumType(OutType(numtype1, numtype2)) numtype;
		typedef FullMatrix structure;
		_matrix<numtype, FullMatrix> dd;

		divMatMatBinOp(const m1 &A, const m2 &B): expr1_(A), expr2_(B),
			dd(expr2_.rows(0), expr2_.cols(0))
		{
			dd.fill(ZeroType<numtype>::zero());
			dd=A*B.inv();
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



// num/matrix
template<class nt1, class nt2, class st2>
inline
_matrix< OutType(nt1, nt2),  typename st2::invers_structure>
operator / (const nt1 &m1, const _matrix<nt2, st2>& m2){
    return inv(m2)*m1;
}

// num/matrixExpr
template<class nt1, class ex2>
inline
_matrix< OutType(nt1, typename ex2::numtype),   typename ex2::invers_structure>
operator / (const nt1 &m1, const _matrixExpr<ex2>& m2){
    return inv(m2)*m1;
}

// matrix/matrix
template<class nt1, class st1, class nt2, class st2>
inline
_matrix<typename _matrixExprBinOpMUL<_matrix<nt1, st1>, _matrix<nt2, typename st2::invers_structure> >::numtype,
        typename _matrixExprBinOpMUL<_matrix<nt1, st1>, _matrix<nt2, typename st2::invers_structure> >::structure>
operator / (const _matrix<nt1, st1>& m1, const _matrix<nt2, st2>& m2){
   return m1*m2.inv();
}

// matrix/=matrix
template<class T, class INstructure>
template<class nt1, class st1>
inline
void _matrix<T, INstructure>::operator
/= (const _matrix<nt1, st1>& m1){
   *this=(*this)*m1.inv();
//   return *this;
}

// matrix/matrixExpr
template<class nt1, class st1, class ex2>
inline
_matrix<typename _matrixExprBinOpMUL<_matrix<nt1, st1>, ex2 >::numtype,
        typename MatOutTypeS( st1,typename ex2::invers_structure) >
operator / (const _matrix<nt1, st1>& m1, const _matrixExpr<ex2>& m2){
    return m1*inv(m2);
}

// matrixExpr/matrix
template<class expr1, class nt2, class st2>
inline
_matrix<typename _matrixExprBinOpMUL<expr1, _matrix<nt2, st2> >::numtype,
        typename  MatOutTypeS( typename expr1::structure, typename st2::invers_structure)  >
operator / (const _matrixExpr<expr1>& m1, const _matrix<nt2, st2>& m2){
    return m1*m2.inv();
}

// matrixExpr/matrixExpr
template<class ex1, class ex2>
inline
_matrix<typename _matrixExprBinOpMUL<ex1, ex2 >::numtype,
       typename  MatOutTypeS(typename ex1::structure, typename ex2::invers_structure)>
operator / (const _matrixExpr<ex1>& m1, const _matrixExpr<ex2>& m2){
   return m1*inv(m2);
}


//**
//*********multiplication**...
// these are expressed only for the
// 'builin' types (double, float...)...this allows OTHER classes
// to provide there own definitions of the operators with
// a matrix without an 'ambigous overload' problem
//**

#define MatBUILTINop(OP, APFUNC, T1) \
template<class ex1>	\
inline	\
_matrixExpr< _matrixExprBinOp< _matrixExpr<ex1>, _matrixExprConst< T1 >,	\
	APFUNC<typename _matrixExpr<ex1>::numtype, T1> > >	\
operator OP (const _matrixExpr<ex1>& m1, T1 m2){	\
    typedef _matrixExprBinOp< _matrixExpr<ex1>, _matrixExprConst<T1>,	\
    	APFUNC<typename _matrixExpr<ex1>::numtype, T1 > > expr;	\
    return _matrixExpr<expr>(expr(m1,m2));	\
}	\
	\
template< class ex2>	\
inline		\
_matrixExpr< _matrixExprBinOp< _matrixExprConst<T1>, _matrixExpr<ex2>,	\
	APFUNC<T1, typename _matrixExpr<ex2>::numtype> > >	\
operator OP ( T1 m1, const _matrixExpr<ex2>& m2){	\
    typedef _matrixExprBinOp< _matrixExprConst< T1 >, _matrixExpr<ex2>,	\
    	ApMultiply<T1, typename _matrixExpr<ex2>::numtype> > expr;	\
    return _matrixExpr<expr>(expr(m1,m2));	\
}	\
	\
template<class nt1, class st1>	\
inline	\
_matrixExpr< _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixExprConst<T1>,	\
	APFUNC<nt1, T1> > >	\
operator OP (const  _matrix<nt1, st1>& m1,  T1 m2){	\
    typedef _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixExprConst<T1>,	\
    	APFUNC<nt1, T1> > expr;	\
    return _matrixExpr<expr>(expr(m1.getRef(),m2));\
}	\
	\
template<class nt2, class st2>	\
inline	\
_matrixExpr< _matrixExprBinOp< _matrixExprConst<T1>, _matrixReference<nt2, st2>,	\
	APFUNC<T1, nt2> > >	\
operator OP(  T1 m1, const _matrix<nt2, st2>& m2){	\
    typedef _matrixExprBinOp< _matrixExprConst< T1 >, _matrixReference<nt2, st2>,	\
    	APFUNC<T1, nt2> > expr;	\
    return _matrixExpr<expr>(expr(m1,m2.getRef()));	\
}	\
	\


//specail case for '/' as num/matrix reuqires an inverse....and is defined above

#define MatBUILTINopDiv(OP, APFUNC, T1) \
template<class ex1>	\
inline	\
_matrixExpr< _matrixExprBinOp< _matrixExpr<ex1>, _matrixExprConst< T1 >,	\
	APFUNC<typename _matrixExpr<ex1>::numtype, T1> > >	\
operator OP (const _matrixExpr<ex1>& m1, T1 m2){	\
    typedef _matrixExprBinOp< _matrixExpr<ex1>, _matrixExprConst<T1>,	\
    	APFUNC<typename _matrixExpr<ex1>::numtype, T1 > > expr;	\
    return _matrixExpr<expr>(expr(m1,m2));	\
}	\
	\
	\
template<class nt1, class st1>	\
inline	\
_matrixExpr< _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixExprConst<T1>,	\
	APFUNC<nt1, T1> > >	\
operator OP (const  _matrix<nt1, st1>& m1,  T1 m2){	\
    typedef _matrixExprBinOp< _matrixReference<nt1, st1>, _matrixExprConst<T1>,	\
    	APFUNC<nt1, T1> > expr;	\
    return _matrixExpr<expr>(expr(m1.getRef(),m2));\
}	\
	\


MatBUILTINop(*, ApMultiply, complex);
MatBUILTINop(+, ApAdd, complex);
MatBUILTINop(-, ApSubtract, complex);
MatBUILTINopDiv(/, ApDivide, complex);

MatBUILTINop(*, ApMultiply, double);
MatBUILTINop(+, ApAdd, double);
MatBUILTINop(-, ApSubtract, double);
MatBUILTINopDiv(/, ApDivide, double);

MatBUILTINop(*, ApMultiply, float);
MatBUILTINop(+, ApAdd, float);
MatBUILTINop(-, ApSubtract, float);
MatBUILTINopDiv(/, ApDivide, float);

MatBUILTINop(*, ApMultiply, int);
MatBUILTINop(+, ApAdd, int);
MatBUILTINop(-, ApSubtract, int);
MatBUILTINopDiv(/, ApDivide, int);

MatBUILTINop(*, ApMultiply, long);
MatBUILTINop(+, ApAdd, long);
MatBUILTINop(-, ApSubtract, long);
MatBUILTINopDiv(/, ApDivide, long);

MatBUILTINop(*, ApMultiply, short);
MatBUILTINop(+, ApAdd, short);
MatBUILTINop(-, ApSubtract, short);
MatBUILTINopDiv(/, ApDivide, short);

MatBUILTINop(*, ApMultiply, char);
MatBUILTINop(+, ApAdd, char);
MatBUILTINop(-, ApSubtract, char);
MatBUILTINopDiv(/, ApDivide, char);

MatBUILTINop(*, ApMultiply, bool);
MatBUILTINop(+, ApAdd, bool);
MatBUILTINop(-, ApSubtract, bool);
MatBUILTINopDiv(/, ApDivide, bool);




//*******************************
//		unuary ops:: (i.e. abs)
//*******************************

//**order of the macro below
// _matrix
// _matrixExpr

#define MatMakeUnary(OP, OPNAME)														\
	template<class nt1, class st1>														\
	inline _matrixExpr<_matrixExprUnrOp<_matrixReference<nt1, st1>, OPNAME<nt1> > > 		\
	OP(const _matrix<nt1, st1>& m1)																\
	{																					\
		typedef _matrixExprUnrOp<_matrixReference<nt1, st1>, OPNAME<nt1>  > expr;							\
		return _matrixExpr<expr>(expr(m1.getRef()));											\
	}																					\
																						\
	template<class ex1>																	\
	inline _matrixExpr<_matrixExprUnrOp<_matrixExpr<ex1>, OPNAME<typename _matrixExpr<ex1>::numtype> > > 	\
	OP(const _matrixExpr<ex1>& m1)	\
	{																					\
		typedef _matrixExprUnrOp<_matrixExpr<ex1>,OPNAME<typename _matrixExpr<ex1>::numtype> > expr;				\
		return _matrixExpr<expr>(expr(m1));													\
	}																					\

MatMakeUnary(abs, ApAbs);
MatMakeUnary(sqrt, ApSqrt);
MatMakeUnary(exp, ApExp);
MatMakeUnary(log, ApLog);
MatMakeUnary(log10, ApLog10);
MatMakeUnary(neg, ApNeg);
MatMakeUnary(Re, ApRe);
MatMakeUnary(Im, ApIm);
MatMakeUnary(real, ApRe);
MatMakeUnary(imag, ApIm);
MatMakeUnary(conj, ApConj);
MatMakeUnary(chop, ApChop);
MatMakeUnary(sin, ApSin);
MatMakeUnary(cos, ApCos);
MatMakeUnary(tan, ApTan);
MatMakeUnary(sinh, ApSin);
MatMakeUnary(cosh, ApCos);
MatMakeUnary(tanh, ApTan);
MatMakeUnary(atan, ApAtan);
MatMakeUnary(acos, ApAcos);
MatMakeUnary(asin, ApAsin);
MatMakeUnary(asinh, ApAsinh);
MatMakeUnary(acosh, ApAcosh);
MatMakeUnary(atanh, ApAtanh);
MatMakeUnary(floor, ApFloor);
MatMakeUnary(ceil, ApCiel);
#ifdef HAVE_FINITE
MatMakeUnary(isnan, ApNan);
#elif HAVE_ISNAN
MatMakeUnary(isnan, ApNan);
#endif


//*****************************
//	handles the 'negation' operator
//	a weird on that does not fit the other cases well
//	as it is a unary o, but is also concidereed an 'operator'
template<class nt1, class st1>
inline _matrixExpr<_matrixExprUnrOp<_matrixReference<nt1, st1>, ApNeg<nt1> > >
operator -(const _matrix<nt1, st1>& m1)
{
	typedef _matrixExprUnrOp<_matrixReference<nt1, st1>, ApNeg<nt1>  > expr;
	return _matrixExpr<expr>(expr(m1.getRef()));
}

template<class ex1>
inline _matrixExpr<_matrixExprUnrOp<_matrixExpr<ex1>, ApNeg<typename _matrixExpr<ex1>::numtype> > >
operator -(const _matrixExpr<ex1>& m1)
{
	typedef _matrixExprUnrOp<_matrixExpr<ex1>,ApNeg<typename _matrixExpr<ex1>::numtype> > expr;
	return _matrixExpr<expr>(expr(m1));
}

//**********************std::max***************
// finds the std::max value inside the matrix
template<class nt1, class st1>
nt1 max(const _matrix<nt1, st1> &in){
	typename st1::iterator iter(in.rows(0), in.cols(0));
	if(in.empty()) return;
	double curv=in(0,0);
	while(iter){
		curv=std::max(curv,in(iter.row(), iter.col()));
		++iter;
	}
	return curv;
}

template<class ex1>
typename _matrixExpr<ex1>::numtype max(const _matrixExpr<ex1>& in)
{
	typename _matrixExpr<ex1>::structure::iterator iter(in.rows(0), in.cols(0));
	typename _matrixExpr<ex1>::numtype curv=in(0,0);
	while(iter){
		curv=std::max(curv,in(iter.row(), iter.col()));
		++iter;
	}
	return curv;
}


//**********************std::max***************
// returns a matrix with the std::max values of the 2 matrices
template<class nt1, class st1, class nt2, class st2>
_matrix<typename _matrixExprBinOpMUL<_matrix<nt1, st1>, _matrix<nt2, st2> >::numtype,
        typename _matrixExprBinOpMUL<_matrix<nt1, st1>, _matrix<nt2, st2> >::structure>
 max( _matrix<nt1, st1> &in,  _matrix<nt2, st2> &m2){

	typedef OutType(nt1, nt2) numtype;
	typedef typename _matrixExprBinOpMUL<_matrix<nt1, st1>, _matrix<nt2, st2> >::structure structure;

   	if(in.rows() != m2.rows() || in.cols() != m2.cols()){
		BLEXCEPTION(" matrix must be the same size ")
	}
    _matrix<numtype, structure> dd(in);

	typename structure::iterator iter(in.rows(0), in.cols(0));
	if(in.empty()) return dd;

	while(iter){
		dd(iter.row(), iter.col())=std::max(m2(iter.row(), iter.col()),in(iter.row(), iter.col()));
		++iter;
	}
	return dd;
}

//**********************std::min***************
// finds the std::min value inside the matrix
template<class nt1, class st1>
nt1 min(const _matrix<nt1, st1> &in){
	typename st1::iterator iter(in.rows(0), in.cols(0));
	if(in.empty()) return ZeroType<nt1>::zero();
	double curv=in(0,0);
	while(iter){
		curv=std::min(curv,in(iter.row(), iter.col()));
		++iter;
	}
	return curv;
}

template<class ex1>
typename _matrixExpr<ex1>::numtype min(const _matrixExpr<ex1>& in)
{
	typename _matrixExpr<ex1>::structure::iterator iter(in.rows(0), in.cols(0));
	double curv=in(0,0);
	while(iter){
		curv=std::min(curv,in(iter.row(), iter.col()));
		++iter;
	}
	return curv;
}

//**********************std::min***************
// returns a matrix with the MINIMAL values of the 2 matrices
template<class nt1, class st1, class nt2, class st2>
_matrix<typename _matrixExprBinOpMUL<_matrix<nt1, st1>, _matrix<nt2, st2> >::numtype,
        typename _matrixExprBinOpMUL<_matrix<nt1, st1>, _matrix<nt2, st2> >::structure>
 min( _matrix<nt1, st1> &in,  _matrix<nt2, st2> &m2){

	typedef OutType(nt1, nt2) numtype;
	typedef typename _matrixExprBinOpMUL<_matrix<nt1, st1>, _matrix<nt2, st2> >::structure structure;

   	if(in.rows() != m2.rows() || in.cols() != m2.cols()){
		BLEXCEPTION(" matrix must be the same size ")
	}
    _matrix<numtype, structure> dd(in);

	typename structure::iterator iter(in.rows(0), in.cols(0));
	if(in.empty()) return dd;

	while(iter){
		dd(iter.row(), iter.col())=std::min(m2(iter.row(), iter.col()),in(iter.row(), iter.col()));
		++iter;
	}
	return dd;
}


//**********************CHOPPER***************
//this keeps litte deviation in numerical error under control
//cuts out any number whose abslute value is below 'cutoff'
template<class nt1, class st1>
template<class T2>
void _matrix<nt1, st1>::chopper(T2 cutoff){

	typename st1::iterator iter(rows(), cols());
	while(iter){
		put(iter.row(), iter.col(),BlochLib::chop( get(iter.row(), iter.col()), cutoff ));
		++iter;
	}
}

template<class nt1, class st1>
template<class T2>
void _matrix<nt1, st1>::chop(T2 cutoff){
	chopper(cutoff);
}

//**********************ReNorm***************
//this keeps litte deviation in numerical error under control
//for matrices that should have value between 0 and 1 (i.e. Unitary)
template<class nt1, class st1>
void _matrix<nt1, st1>::ReNorm(double cutoff){
	if(BiggerType(nt1, complex)){
		double tmpr=0., tmpi=0.;
		typename st1::iterator iter(rows(), cols());
		while(iter){
			tmpr=Re(get(iter.row(),iter.col()));
			tmpi=Im(get(iter.row(),iter.col()));
			if(abs(tmpr)<=cutoff)	{	tmpr=0.;	}
			else if(abs(tmpr)>1.)	{	tmpr=sign(tmpr);	}
			if(abs(tmpi)<=cutoff)	{	tmpi=0.;	}
			else if(abs(tmpi)>1.)	{	tmpi=sign(tmpi);	}
			put(iter.row(), iter.col(),complex(tmpr, tmpi));
			++iter;
		}
	}else{
		nt1 tmpr=0;
		typename st1::iterator iter(rows(), cols());
		while(iter){
			tmpr=get(iter.row(),iter.col());
			if(abs(tmpr)<=cutoff)	{	tmpr=0.;	}
			else if(abs(tmpr)>1.)	{	tmpr=sign(tmpr);	}
			put(iter.row(), iter.col(),tmpr);
			++iter;
		}
	}
}


//*******************_matrix Output****************

template<class T, class structure>
std::ostream &operator<<(std::ostream &otr, _matrix<T, structure> oo){
	int i, j;
	otr<<MatName(oo.type());

	otr<<" matrix: "<<oo.rows()<<"x"<<oo.cols()<<std::endl;
	otr<<"[ ";
	for(i=0;i<oo.rows();i++){
		otr<<"[ ";
		for(j=0;j<oo.cols();j++){
			otr<<setw(13)<<oo(i,j)<<" ";
		}
		otr<<" ]"<<std::endl;
	}
	otr<<" ]"<<std::endl;
	return otr;
}

template<class structure>
std::ostream &operator<<(std::ostream &otr, _matrix<complex, structure> oo){
	int i, j;
	otr<<MatName(oo.type());

	otr<<" matrix: "<<oo.rows()<<"x"<<oo.cols()<<std::endl;
	otr<<"[ ";
	for(i=0;i<oo.rows();i++){
		otr<<"[ ";
		for(j=0;j<oo.cols();j++){
			otr<<"complex"<<oo(i,j)<<" ";
		}
		otr<<" ]"<<std::endl;
	}
	otr<<" ];"<<std::endl;
	return otr;
}

template<class T>
std::ostream &operator<<(std::ostream &otr, _matrixExpr<T> oo){
	int i, j;
	typename _matrixExpr<T>::structure s;
	otr<<MatName(s.type());
	otr<<" matrix: "<<oo.rows(1)<<"x"<<oo.cols(1)<<std::endl;
	otr<<"[ ";
	for(i=0;i<oo.rows(1);i++){
		otr<<"[ ";
		for(j=0;j<oo.cols(1);j++){
			otr<<setw(13)<<oo(i,j)<<" ";
		}
		otr<<" ]"<<std::endl;
	}
	otr<<" ]"<<std::endl;
	return otr;
}

END_BL_NAMESPACE


#include "container/matrix/matmatmul.h" 	//matrix*matrix functions
#include "container/matrix/matvec.h"		//matrix vector operations
#include "container/matrix/matmat.h" 	//holds specialized matrix matrix type operations
#include "container/matrix/matfft.h" 	//matrix FFT operations
#include "container/matrix/matdiagonalize.h" //diagonalization routines
#include "container/matrix/matexp.h" //diagonalization routines
#include "container/matrix/LU.h"
#include "container/matrix/matassign.h"		//assignment routines
#include "container/matrix/matqr.h"  	//QR decomposotion and gram schmit orthonorm

#endif

