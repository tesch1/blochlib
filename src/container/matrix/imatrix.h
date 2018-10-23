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
//	IDENTITY MATRIX::a 'trait' type class for Matrix
//
/***************************************************/

#ifndef _imatrix_h_
#define _imatrix_h_ 1

#include "container/rankType.h"
#include "container/matrix/genmatrix.h"
#include "container/matrix/matrixconfig.h"

BEGIN_BL_NAMESPACE


//the iterator class for the IdentityMatrix
class IdentityMatrixItter : public GeneralMatrixItter {
	private:
		void iterr() const {
			BLEXCEPTION(" iterator bounds already reached limit")
		}
	public:
		IdentityMatrixItter() : GeneralMatrixItter() {}
		IdentityMatrixItter(int ro, int co): GeneralMatrixItter(ro,co) {}

		void operator++(){
			if(!underlimit_) iterr();
			++i;
			++j;
			if(j==cols_ || i==rows_)	underlimit_=false;
		}
};


//the Identity matrix structure

class IdentityMatrix : public GeneralMatrix{
	private:
		void testsqare() const {
			if(!issquare()){
				BLEXCEPTION(" Identity Matrix Must be square...")
			}
		}

	public:
		typedef IdentityMatrixItter iterator;
		typedef IdentityMatrix structure;
		typedef IdentityMatrix invers_structure;

		MatrixType type()const{	return type_;	}

		IdentityMatrix(): GeneralMatrix(){
			type_=Midentity;
			numelements_=0;
		}

		IdentityMatrix(int ro): GeneralMatrix(ro,ro){
			type_=Midentity;
			numelements_=1;
		}

		IdentityMatrix(int ro, int co): GeneralMatrix(ro,co){
			testsqare();
			type_=Midentity;
			numelements_=1;
		}

		int numElements(int ro, int col){	return 1;	}
		int end(){	return numelements_;	}

		void resize(int nr, int nc){
			GeneralMatrix::resize(nr, nc);
			testsqare();
			numelements_=1;
		}
		//class should return a '0' of the type
		template<class T>
		inline T get(int r, int c, const T *data)const {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			if(r==c) return OneType<T>::one();
			else return ZeroType<T>::zero();
		}

		template<class T>
		inline T &get(int r, int c, T *data) {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			if(r==c) return OneType<T>::one();
			else return ZeroType<T>::zero();
		}

		//returns the posision in a ficticios data vector
		inline int position(int i, int j){
#ifdef MatDeBug
			lenerr(i,j);
#endif
			if(i==j && i<numelements_){ return 0;	}
			return 0;
		}

		template<class T>
		inline void put(int i, int j, T *data, T thing)const{
#ifdef MatDeBug
			lenerr(i,j);
#endif
			if(i==j){ data[0]=thing;	}
		}
};


END_BL_NAMESPACE



#endif

