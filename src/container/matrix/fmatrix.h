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
//	FULL MATRIX::a 'trait' type class for Matrix
//
/***************************************************/

#ifndef _fmatrix_h_
#define _fmatric_h_ 1

#include "container/matrix/genmatrix.h"
#include "container/matrix/matrixconfig.h"

BEGIN_BL_NAMESPACE


//the iterator class for the FullMatrix
class FullMatrixItter : public GeneralMatrixItter {
	private:
		void iterr() const {
			BLEXCEPTION(" iterator bounds already reached limit")
		}
	public:
		FullMatrixItter() : GeneralMatrixItter() {}
		FullMatrixItter(int ro, int co): GeneralMatrixItter(ro,co) {}

		void operator++(){
			if(!underlimit_) iterr();
			++element_;
			++j;
			if(j==cols_){
				j=0;
				++i;
				if(i==rows_) 	underlimit_=false;
			}
		}
};



//the full matrix structure
class FullMatrix : public GeneralMatrix {
	private:

	public:
		typedef FullMatrixItter iterator;
		typedef FullMatrix structure;
		typedef FullMatrix invers_structure;


		MatrixType type()const{	return type_;	}

		FullMatrix(): GeneralMatrix(){
			type_=Mfull;
			numelements_=0;
		}

		FullMatrix(int ro, int co) : GeneralMatrix(ro, co){
			type_=Mfull;
			numelements_=ro*co;
		}

		int numElements(int ro, int col){	return ro*col;	}
		int end(){	return numelements_;	}

		void resize(int nr, int nc){
			GeneralMatrix::resize(nr, nc);
			numelements_=nr*nc;
		}

		template<class T>
		inline T get(int r, int c, const T *data)const {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			return data[c*rows_+r];
		}

		template<class T>
		inline T &get(int r, int c, T *data) {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			return data[c*rows_+r];
		}

		//returns the posision in a ficticios data vector
		inline int position(int i, int j)const{
#ifdef MatDeBug
			lenerr(i,j);
#endif
			{ return j*rows_+i;	}
		}

		template<class T>
		inline void put(int i, int j, T *data, T thing)const{
#ifdef MatDeBug
			lenerr(i,j);
#endif
			data[j*rows_+i]=thing;
		}
};

END_BL_NAMESPACE



#endif


