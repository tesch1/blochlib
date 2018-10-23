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
//	Hermietian MATRIX::a 'trait' type class for Matrix
//		(like Symmetric, w/ the conjugate rule)
//
/***************************************************/

#ifndef _smatrix_h_
#define _smatric_h_ 1

#include "container/matrix/genmatrix.h"
#include "container/matrix/matrixconfig.h"

BEGIN_BL_NAMESPACE


//the iterator class for the FullMatrix
class HermitianMatrixItter : public GeneralMatrixItter {
	private:
		void iterr() const {
			BLEXCEPTION(" iterator bounds already reached limit")
		}
	public:
		HermitianMatrixItter() : GeneralMatrixItter() {}
		HermitianMatrixItter(int ro, int co): GeneralMatrixItter(ro,co) {}

		void operator++(){
			if(!underlimit_) iterr();
			++element_;
			++j;
			if (j > i){
				j = 0;
				++i;
				if (i == rows_)	underlimit_ = false;
       		}
		}
};



//the full matrix structure
class HermitianMatrix : public GeneralMatrix {
	private:
		void testsquare() const {
			if(!issquare()){
				BLEXCEPTION(" Hermitian Matrix Must be square...")
			}
		}

	public:
		typedef HermitianMatrixItter iterator;
		typedef HermitianMatrix structure;
		typedef FullMatrix invers_structure;

		const void HermErr() const
		{
			BLEXCEPTION(std::string(" Hermitian Matrix Diagonal Elements contains IMAGINARY...") +
					std::string(" Components...this is strictly unallowed..change matrix type to 'Full'..."))
		}

		MatrixType type()const{ return type_;	}

		HermitianMatrix(): GeneralMatrix(){
			type_=Mhermitian;
			numelements_=0;
		}

		HermitianMatrix(int ro): GeneralMatrix(ro,ro){
			type_=Mhermitian;
			numelements_=ro*(ro+1)/2;
		}

		HermitianMatrix(int ro, int co) : GeneralMatrix(ro, co){
			testsquare();
			type_=Mhermitian;
			numelements_=ro*(ro+1)/2;
		}

		int numElements(int ro, int col){	return ro*(ro+1)/2;	}
		int end(){	return numelements_;	}

		void resize(int nr){
			GeneralMatrix::resize(nr, nr);
			testsquare();
			numelements_=nr*(nr+1)/2;
		}

		void resize(int nr, int nc){
			GeneralMatrix::resize(nr, nc);
			testsquare();
			numelements_=nr*(nr+1)/2;
		}

		//class should return a '0' of the type
		template<class T>
		inline T get(int r, int c, const T *data)const {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			if(r>=c)	return data[position(r,c)];
			//static T value_;	//used becuase we need to take a 'conj' making refernece passing impossible
			//value_=conj(data[position(r,c)]);
			return conj(data[position(r,c)]);
		}

		template<class T>
		inline T &get(int r, int c, T *data) {
#ifdef MatDeBug
			lenerr(r,c);
#endif
			if(r>=c)	return data[position(r,c)];
			static T value_;	//used becuase we need to take a 'conj' making refernece passing impossible
			value_= conj(data[position(r,c)]);
			return value_;
		}


		//returns the posision in a ficticios data vector
		inline int position(int i, int j)const{
#ifdef MatDeBug
			lenerr(i,j);
#endif
			if (i >= j)   return i*(i+1)/2 + j;
			else          return j*(j+1)/2 + i;
			/*if (i == j)		return i*cols()-(i*(i-1))/2;
			else if(i<j)	return i*cols()-(i*(i-1))/2+j-i;
			else			return j*cols()-(j*(j-1))/2+i-j;*/
		}

		template<class T>
		inline void put(int i, int j, T *data, T thing)const{
#ifdef MatDeBug
			lenerr(i,j);
#endif
			//if (i == j)  data[position(i,j)]=Re(thing);
			if(i>=j) data[position(i,j)]=thing;
			else         data[position(i,j)]=conj(thing);

		}
};

END_BL_NAMESPACE



#endif


