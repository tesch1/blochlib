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
//	DIAGONAL MATRIX::a 'trait' type class for Matrix
//
/***************************************************/

#ifndef _dmatrix_h_
#define _dmatrix_h_ 1

#include "container/rankType.h"
#include "container/matrix/matrixconfig.h"
#include "container/matrix/genmatrix.h"

BEGIN_BL_NAMESPACE


//the iterator class for the FullMatrix
class DiagonalMatrixItter : public GeneralMatrixItter {
	private:
		void iterr() const {
			BLEXCEPTION(" iterator bounds already reached limit")
		}
	public:
		DiagonalMatrixItter() : GeneralMatrixItter() {}
		DiagonalMatrixItter(int ro, int co): GeneralMatrixItter(ro,co) {}

		void operator++(){
			if(!underlimit_) iterr();
			++element_;
			++i;
			++j;
			if(j==cols_ || i==rows_)	underlimit_=false;
		}
};


//the diagonal matrix structure
class DiagonalMatrix : public GeneralMatrix {
	private:
		void testsqare() const {
			if(!issquare()){
				BLEXCEPTION(" Diagonal Matrix Must be square...")
			}
		}

	public:
		typedef DiagonalMatrixItter iterator;
		typedef DiagonalMatrix structure;
		typedef DiagonalMatrix invers_structure;


		MatrixType type()const{	return type_;	}

		DiagonalMatrix(): GeneralMatrix(){
			type_=Mdiagonal;
			numelements_=0;
		}

		DiagonalMatrix(int ro): GeneralMatrix(ro,ro){
			type_=Mdiagonal;
			numelements_=1;
		}

		DiagonalMatrix(int ro, int co): GeneralMatrix(ro,co){
			testsqare();
			type_=Mdiagonal;
			numelements_=ro;
		}

		int numElements(int ro, int col){	return ro;	}
		int end(){	return numelements_;	}

		void resize(int nr, int nc){
			GeneralMatrix::resize(nr, nc);
			testsqare();
			numelements_=nr;
		}

		//class should return a '0' of the type
		template<class T>
		inline T get(int r, int c, const T *data)const {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			if(r==c) return data[r];
			else return ZeroType<T>::zero();
		}

		template<class T>
		inline T &get(int i, int j, T *data) {
#ifdef MatDeBug
			lenerr(i,j);
#endif
			if(i==j) return data[i];
			else  return ZeroType<T>::zero();
		}

		template<class T>
		inline T &get(int i, T *data) {
#ifdef MatDeBug
			//lenerr(i,j);
#endif
			return data[i];
		}

		//returns the posision in a ficticios data vector
		inline int position(int i, int j){
#ifdef MatDeBug
			lenerr(i,j);
#endif
			if(i==j && i<numelements_){ return i;	}
			return 0;
		}

		template<class T>
		inline void put(int i, int j, T *data, const T &thing) const{
#ifdef MatDeBug
			lenerr(i,j);
#endif
			if(i==j){ data[i]=thing;	}
		}
};


END_BL_NAMESPACE


#endif


