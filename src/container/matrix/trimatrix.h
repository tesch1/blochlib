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

/***************************************************/
//
//	TriDiagonal MATRIX::a 'trait' type class for Matrix
//		has elements on the diagonl and +/- a column
//
/***************************************************/

#ifndef _trimatrix_h_
#define _trimatric_h_ 1

#include "container/matrix/genmatrix.h"
#include "container/matrix/matrixconfig.h"

BEGIN_BL_NAMESPACE


//the iterator class for the FullMatrix
class TriDiagonalMatrixItter : public GeneralMatrixItter {
	private:
		void iterr() const {
			BLEXCEPTION(" iterator bounds already reached limit")
		}
	public:
		TriDiagonalMatrixItter() : GeneralMatrixItter() {}
		TriDiagonalMatrixItter(int ro, int co): GeneralMatrixItter(ro,co) {}

		void operator++(){
			if(!underlimit_) iterr();
			++element_;
			++j;
			if (j > i+1){
				++i;
				j = i-1;
				if (i == rows_ )	underlimit_ = false;
       		}
       //		cout<<i<<" "<<j<<" "<<element_<<endl;
		}
};



//the full matrix structure
class TriDiagonalMatrix : public GeneralMatrix {
	private:
		void testsquare() const {
			if(!issquare()) {
				BLEXCEPTION(" TriDiagonalMatrix Must be square...")
			}
		}
	public:
		typedef TriDiagonalMatrix structure;
		typedef FullMatrix invers_structure;
		typedef TriDiagonalMatrixItter iterator;

		MatrixType type()const{ return type_;	}

		TriDiagonalMatrix(): GeneralMatrix(){
			type_=Mtridiagonal;
			numelements_=0;
		}

		TriDiagonalMatrix(int ro): GeneralMatrix(ro,ro){
			type_=Mtridiagonal;
			numelements_=ro*3;
		}

		TriDiagonalMatrix(int ro, int co) : GeneralMatrix(ro, co){
			testsquare();
			type_=Mtridiagonal;
			numelements_=ro*3;
		}

		int numElements(int ro, int col){	return ro*3;	}
		int end(){	return numelements_;	}

		void resize(int nr){
			GeneralMatrix::resize(nr, nr);
			testsquare();
			numelements_=nr*3;
		}

		void resize(int nr, int nc){
			GeneralMatrix::resize(nr, nc);
			testsquare();
			numelements_=nr*3;
		}

		//class should return a '0' of the type
		template<class T>
		inline T get(int r, int c, const T *data)const {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			if(r==c || r==c+1 || r==c-1) return data[position(r,c)];
			return ZeroType<T>::zero();
		}

		template<class T>
		inline T &get(int r, int c, T *data) {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			if(r==c || r==c+1 || r==c-1) return data[position(r,c)];
			return ZeroType<T>::zero();
		}
		//gets the diagonal element....
		template<class T>
		inline T getD(int r, const T *data)const {
#ifdef MatDeBug
			lenerr(r,r);
#endif
			return data[position(r,r)];

		}

		//gets the diagonal element....
		template<class T>
		inline T &getD(int r, const T *data) {
#ifdef MatDeBug
			lenerr(r,r);
#endif
			return data[position(r,r)];
		}

		//gets the Right Column element from the diagonal....
		template<class T>
		inline T getR(int r, const T *data)const {
#ifdef MatDeBug
			lenerr(r,r+1);
#endif
			return data[position(r,r+1)];
		}

		//gets the Right Column element from the diagonal....
		template<class T>
		inline T &getR(int r, const T *data) {
#ifdef MatDeBug
			lenerr(r,r+1);
#endif
			return data[position(r,r+1)];
		}

		//gets the LEFT Column element from the diagonal....
		template<class T>
		inline T getL(int r, const T *data)const {
#ifdef MatDeBug
			lenerr(r,r-1);
#endif
			return data[position(r,r-1)];
		}

		//gets the LEFT Column element from the diagonal....
		template<class T>
		inline T &getL(int r, const T *data) {
#ifdef MatDeBug
			lenerr(r,r-1);
#endif
			return data[position(r,r-1)];
		}


		//returns the posision in a ficticios data vector
		inline int position(int i, int j)const{
#ifdef MatDeBug
			lenerr(i,j);
#endif
			if (i == j)   return i;
			else if(j==i-1) return i+rows();
			else if(j==i+1) return i+2*rows();
			else          return 0;
		}

		template<class T>
		inline void put(int i, int j, T *data, const T &thing) const{
#ifdef MatDeBug
			lenerr(i,j);
#endif
			data[position(i,j)]=thing;
		}

	//put the Diagonal Element
		template<class T>
		inline void putD(int i, T *data, const T &thing) const{
#ifdef MatDeBug
			lenerr(i,i);
#endif
			data[position(i,i)]=thing;
		}

		//put the LEFT COLUMN Element
		template<class T>
		inline void putL(int i, T *data, const T &thing) const{
#ifdef MatDeBug
			lenerr(i,i-1);
#endif
			data[position(i,i-1)]=thing;
		}

		//put the RIGHT COLUMN Element
		template<class T>
		inline void putR(int i, T *data, const T &thing) const{
#ifdef MatDeBug
			lenerr(i,i+1);
#endif
			data[position(i,i+1)]=thing;
		}
};


END_BL_NAMESPACE


#endif


