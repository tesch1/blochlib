/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-20-01
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
//	GENERAL BASE MATRIX::a 'trait' base class for the Matrix traits
//
/***************************************************/

#ifndef _genmatrix_h_
#define _genmatrix_h_ 1

#include "container/rankType.h"


BEGIN_BL_NAMESPACE



//*********************************
// global enumerator that tells me what type of matrix
// we are ....
enum MatrixType{Mgeneral,Mfull, Mdiagonal, Midentity, Msymmetric, Mhermitian, Mtridiagonal};

class FullMatrix;
class DiagonalMatrix;
class SymmetricMatrix;
class IdentityMatrix;
class GeneralMatrix;
class HermitianMatrix;
class TriDiagonalMatrix;

//**********************************
// spits out the 'std::string version of the emun type

inline static std::string MatName(MatrixType in) {
	if(in==Mfull) 		return "Full";
	else if(in==Msymmetric)	return "Symmetric";
	else if(in==Mdiagonal)	return "Diagonal";
	else if(in==Midentity)	return "Identity";
	else if(in==Mhermitian)	return "Hermitian";
	else if(in==Mtridiagonal) return "Tri-Diagonal";
	else					return "General";
}

//***********************************
// used to see if the 'left hand side' is of 'greater, more general type the the rhs
// used for assignement checking...i.e. <hermitian>=<full> would be bad
// but <full>=<diag> is okay

#ifdef MatAssignCheck

static bool CanAssign(MatrixType lhs, MatrixType rhs) {
	if(lhs==rhs){						return true;	}
	else if(lhs==Mfull || lhs==Mgeneral){	return true;	}
	else if(lhs==Mhermitian){
		if(rhs==Mdiagonal || rhs==Midentity || rhs==Mtridiagonal){		return true;	}
		else{ 							return false;	}
	}else if(lhs==Msymmetric){
		if(rhs==Mdiagonal|| rhs==Midentity ||  rhs==Mtridiagonal){				return true;	}
		else{							return false;	}
	}else if(lhs==Mtridiagonal){
		if(rhs==Mdiagonal|| rhs==Midentity ){				return true;	}
		else{							return false;	}
	}else if(lhs==Mdiagonal){
		if(rhs==Midentity){				return true;	}
		else{							return false;	}
	}else{
		return false;
	}
}

#endif

//*****************************************
//  given two matrix types and we do an operation
//	we may need to promote the matrix type based on
//  an operation performed on it, thus these little
//  guys below do just that
//  use "MatOutType(mat1, mat2)" to get the correct outtype

template<class OBJ>
struct MatPromoteObj{
	static const  int rank;
	typedef OBJ structure;
};


#define MatObj(OBJ, NUM) \
template<>		\
struct MatPromoteObj<class OBJ>{	\
	static const int rank=NUM;	\
	typedef OBJ structure;	\
}; \

#define MatObjNum(OBJ, NUM) \
template<>		\
struct MatPromoteObj<OBJ>{	\
	static const int rank=NUM;	\
	typedef OBJ structure;	\
}; \


template<class Num_t, class Struct_T>
class _matrix;

template<class INnumtype, class INstructure>
class _matrixReference ;

template<class INtype>
class _matrixExprConst;


MatObj(GeneralMatrix, 1000)
MatObj(FullMatrix, 800)
MatObj(HermitianMatrix, 700)
MatObj(SymmetricMatrix, 600)
MatObj(TriDiagonalMatrix, 500)
MatObj(DiagonalMatrix, 400)
MatObj(IdentityMatrix, 100)
MatObjNum(bool, 0)
MatObjNum(char, 1)
MatObjNum(short, 2)
MatObjNum(unsigned int, 3)
MatObjNum(int, 4)
MatObjNum(long int, 5)
MatObjNum(float, 6)
MatObjNum(double, 7)
//MatObjNum(long double, 8);
MatObjNum(Complex<float>, 9)
MatObjNum(Complex<double>, 10)

template<class T, class INstructure>
struct MatPromoteObj<_matrix<T, INstructure> >{
	static const int rank=MatPromoteObj<INstructure>::rank;
	typedef INstructure structure;
};

template<class INnumtype, class INstructure>
struct MatPromoteObj<_matrixReference<INnumtype, INstructure> >{
	static const int rank=MatPromoteObj<INstructure>::rank;
	typedef INstructure structure;
};

template<class INtype>
struct MatPromoteObj<_matrixExprConst<INtype> >{
	static const int rank=MatPromoteObj<INtype>::rank;
	typedef INtype structure;
};

template<class T1, class T2, int T1_wins> //t1 wins
struct MatUp2{
	typedef T1 outtype;
};

template<class T1, class T2>	//t2 wins
struct MatUp2<T1, T2, 0>{
	typedef T2 outtype;
};

template<class T1, class T2>
struct MatPromote{
	typedef MatPromoteObj<typename T1::structure> struct1;
	typedef MatPromoteObj<typename T2::structure> struct2;

	enum{
		aresame = int(struct1::rank) == int(struct2::rank),

		T1_wins = int(struct1::rank) > int(struct2::rank),

		T2_wins = int(struct1::rank) < int(struct2::rank)
	};

	enum{
		uptoBool= (aresame) ? 1 : ((T1_wins) ? 1: 0)
	};
	typedef typename MatUp2<typename struct1::structure,typename  struct2::structure, uptoBool>::outtype outtype;
};

template<class T1, class T2>
struct MatPromoteS{
	typedef MatPromoteObj<T1> struct1;
	typedef MatPromoteObj<T2> struct2;

	enum{
		aresame = int(struct1::rank) == int(struct2::rank),

		T1_wins = int(struct1::rank) > int(struct2::rank),

		T2_wins = int(struct1::rank) < int(struct2::rank)
	};

	enum{
		uptoBool= (aresame) ? 1 : ((T1_wins) ? 1: 0)
	};
	typedef typename MatUp2<typename struct1::structure,typename  struct2::structure, uptoBool>::outtype outtype;
};

#define MatOutType(OBJ1, OBJ2) MatPromote<OBJ1, OBJ2>::outtype
#define MatOutTypeS(OBJ1, OBJ2) MatPromoteS<OBJ1, OBJ2>::outtype



//the General itteratorsimply defines '++'  on the
// row/column steps...the ++ is simply advance 'j' until we hit a
// column size limit, then advance the 'i' until we run out of
// rows..the inhereted classes contain the '++' itterator operator

class GeneralMatrixItter {
	protected:
		int rows_;
		int cols_;
		bool underlimit_;
		int i;
		int j;
		int element_;

	public:

		GeneralMatrixItter():
			rows_(0), cols_(0), underlimit_(true), i(0), j(0), element_(0){}
		GeneralMatrixItter(int ro, int co){
			rows_=ro; cols_=co; i=0; j=0;
			if(ro+co==0){
				underlimit_=false;
			}else{
				underlimit_=true;
			} element_=0;
		}

		//specifies the actual spot in the matrix
		int row(){ return i;	}
		int col(){	return j;	}
		int column(){	return j;	}
		bool underlimit(){	return underlimit_;	}
		int element()const{	return element_;	}
		int Element()const{	return element_;	}
		int ele()const{		return element_;	}
		operator bool() const{ return underlimit_;	}
};


class GeneralMatrix {
	friend class FullMatrix;
	friend class HermitianMatrix;
	friend class DiagonalMatrix;
	friend class IdentityMatrix;
	friend class SymmetricMatrix;
	friend class TriDiagonalMatrix;
	private:
		int rows_;
		int cols_;
		int numelements_;

	public:
		MatrixType type_;
		typedef GeneralMatrixItter iterator;
		typedef GeneralMatrix invers_structure;

		GeneralMatrix(){
			rows_=0; cols_=0;
			type_=Mgeneral;
		}

		GeneralMatrix(int ro, int co){
			rows_=ro;
			cols_=co;
			type_=Mgeneral;
		}

		int cols()const{ return cols_;	}
		int columns()const{	return cols_;	}
		int rows()const{	return rows_;	}
		int cols(int guessC){ return cols_;	}
		int columns(int guessC){	return cols_;	}
		int rows(int guessR){	return rows_;	}

		int numElements()const{	return numelements_;	}
		int numelements()const{	return numelements_;	}
		int numEle()const{	return numelements_;	}
		int numElements(){	return numelements_;	}
		int numelements(){	return numelements_;	}
		int numEle(){	return numelements_;	}

		void resize(int nr, int nc){
			rows_=nr;
			cols_=nc;
		}
		bool issquare()const{ return rows_==cols_;	}

		void lenerr(int r, int c)const{
			if(r<0){
				std::cerr<<std::endl<<"error: Matrix[]"<<std::endl;
				std::cerr<<" row element accsess below 0"<<std::endl;
				    exit(1);
			}else if(r>=rows_){
				std::cerr<<std::endl<<"error: Matrix[]"<<std::endl;
				std::cerr<<" row element accsess above row count"<<std::endl;
				    exit(1);
			}else if(c<0){
				std::cerr<<std::endl<<"error: Matrix[]"<<std::endl;
				std::cerr<<" column element accsess below 0"<<std::endl;
				    exit(1);
			}else if(c>=cols_){
				std::cerr<<std::endl<<"error: Matrix[]"<<std::endl;
				std::cerr<<" column element accsess above column count"<<std::endl;
				    exit(1);
			}
		}
};


END_BL_NAMESPACE



#endif

