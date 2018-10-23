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


/* The assignment operations...this acts as a little compiler program
	(i.e. at compile time this is acctually unrolled completely)
	becuase each of the classes below must express itself at compilation
	the compile acctually unrolls the entire list inside the matrix
	upon assignment....*/


#ifndef _matassign_h_
#define _matassign_h_1


BEGIN_BL_NAMESPACE


template<class T, class INstructure>
void _matrix<T, INstructure>::operator=(const _matrix<T, INstructure> &INexpr)
{
//	this->resize(INexpr.rows(), INexpr.cols());					//resize 'this' to the expr size
//	iterator iter(rows(), cols());
//	while (iter){
//		put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
//		++iter;
//	}
//	_matrix<T, INstructure> &eee=const_cast< _matrix<T, INstructure>& >(INexpr);
//	changeBlock(eee,0);
//	newBlock(eee, 0);
//	std::memcpy(data(),INexpr.data(),rows()*cols()*sizeof(T));
	copy(INexpr);
//	return *this;
}

template<class T, class INstructure>
template<class nt2, class st2>
void _matrix<T, INstructure>::operator=(const _matrix<nt2, st2> &INexpr)
{
	//typedef structure outstruct ;
	//typedef INexpr::structure instruct;

	T *s1=NULL;
	st2 *intype=NULL;

#ifdef MatAssignCheck
	if(!CanAssign(type(),INexpr.type())) structureError(MatName(INexpr.type()));	//check to make sure structures are compatable
#endif
	this->resize(INexpr.rows(), INexpr.cols());					//resize 'this' to the expr size

	AssignMat(*this, INexpr, s1, intype);
//	return *this;
}

//the inclass assignment operator

template<class T, class INstructure>
template<class theexpr>
void _matrix<T, INstructure>::operator=(const _matrixExpr<theexpr> &INexpr)
{
	typedef structure outstruct ;

	typedef typename _matrixExpr<theexpr>::structure instruct;
	//typedef _matrix<nt1, st1> mymat;
	typedef  _matrixExpr<theexpr> myexpr;

//	typename  MatOutType(instruct, outstruct) *s1= new typename  MatOutType(instruct, outstruct);
	typename  MatOutType(instruct, outstruct) *s1=NULL;
	//typename _matrixExpr<theexpr>::structure *intype= new typename _matrixExpr<theexpr>::structure;
	typename _matrixExpr<theexpr>::structure *intype=NULL;

#ifdef MatAssignCheck
	if(!CanAssign(type(),s1.type())) structureError(MatName(s1.type()));	//check to make sure structures are compatable
#endif

	this->resize(INexpr.rows(rows()), INexpr.cols(cols()));					//resize 'this' to the expr size

	AssignMat(*this, INexpr, s1, intype);
	//delete s1;
	//delete intype;
//	return *this;
}

//the templated version of each type comparison

//the default one
//Full=Full uses default
//Hermitian=Hermiaitan uses default
//Symmetric=Symmetric uses default
//Diagonal=Diagonal uses default
//Identity=Identity uses default
//Diagonal=Identity uses default

template<class outStructure, class inStructure, class mat, class expr>
inline void AssignMat( mat &in,const expr &INexpr, const outStructure *outst, const inStructure *inst)
{
	typename mat::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		++iter;
	}
	/*static int i,j;
	for(i=0;i<in.rows();++i){
		for(j=0;j<in.cols();++j){
			in(i,j)=INexpr(i,j);
		}
	}*/
}



//Full=Hermitian
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const HermitianMatrix *ins)
{
	typename HermitianMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()!=iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
			in.put(iter.col(), iter.row(), conj(INexpr(iter.row(), iter.col())));
		}else{
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}

//Full=Symmetric
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const SymmetricMatrix *ins)
{
	typename SymmetricMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()!=iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
			in.put(iter.col(), iter.row(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}


//Full=Identity
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const IdentityMatrix *ins)
{
	typename FullMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(),ZeroType<typename mat::numtype>::zero());
		}
		++iter;
	}
}

//Full=Diagonal
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const DiagonalMatrix *ins)
{
	typename FullMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(),ZeroType<typename mat::numtype>::zero());
		}
		++iter;
	}
}

//Hermitian=Diagonal
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr, const HermitianMatrix *outs, const DiagonalMatrix *ins)
{
	typename HermitianMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(),ZeroType<typename mat::numtype>::zero());
		}
		++iter;
	}
}

//Hermitian=Identity
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr,const HermitianMatrix *outs, const IdentityMatrix *ins)
{
	typename HermitianMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(),ZeroType<typename mat::numtype>::zero());
		}
		++iter;
	}
}

//Hermitian=Symmetric
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr,const HermitianMatrix *outs, const SymmetricMatrix *ins)
{
	typename HermitianMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}

//Hermitian=Hermitian
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr,const HermitianMatrix *outs, const HermitianMatrix *ins)
{
	typename HermitianMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			//if(Im(INexpr(iter.row(), iter.col()))!=0) ins->HermErr(); //check for imaginary diagonal...
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}


//Symmetric=Identity
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr,const SymmetricMatrix *outs, const IdentityMatrix *ins)
{
	typename SymmetricMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(),ZeroType<typename mat::numtype>::zero());
		}
		++iter;
	}
}

//Symmetric=Diagonal
template<class mat, class expr>
inline void AssignMat(const mat &in,const expr &INexpr,const SymmetricMatrix *outs, const DiagonalMatrix *ins)
{
	typename SymmetricMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(),ZeroType<typename mat::numtype>::zero());
		}
		++iter;
	}
}



template<class T, class INstructure>
template<class nt1, class st1>
inline void _matrix<T, INstructure>::operator += (const _matrix<nt1, st1> &m2)
{
	typedef structure outstruct ;
	typedef typename _matrix<nt1, st1>::structure instruct;

	static typename MatOutType(instruct, outstruct) s1;
	static instruct intype;

#ifdef MatAssignCheck
	if(!CanAssign(type(),s1.type())) structureError(MatName(s1.type()));	//check to make sure structures are compatable
#endif

	AddAssignMat(*this, m2, &s1, &intype);
//	return *this;
}

template<class T, class INstructure>
template<class expr>
inline void _matrix<T, INstructure>::operator += (const _matrixExpr<expr> &m2)
{
	typedef structure outstruct ;
	typedef typename _matrixExpr<expr>::structure instruct;
	static typename MatOutType(instruct, outstruct) s1;
	static typename _matrixExpr<expr>::structure intype;

#ifdef MatAssignCheck
	if(!CanAssign(type(),s1.type())) structureError(MatName(s1.type()));	//check to make sure structures are compatable
#endif
	AddAssignMat(*this, m2, &s1, &intype);
//	return *this;
}



//for "+=" operations

//the default one
//Full+=Full uses default
//Hermitian+=Hermiaitan uses default
//Symmetric+=Symmetric uses default
//Diagonal+=Diagonal uses default
//Identity+=Identity uses default (not really legal)
//Diagonal+=Identity uses default
//Hermitian+=Symmetric uses default

template<class outStructure, class inStructure, class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr, const outStructure *outst, const inStructure *inst)
{
	typename outStructure::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(),in.get(iter.row(),iter.col())+INexpr(iter.row(), iter.col()));
		++iter;
	}
}



//Full=Hermitian
template<class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const HermitianMatrix *ins)
{
	typename HermitianMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()!=iter.col()){
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
			in.put(iter.col(), iter.row(), in.get(iter.col(), iter.row())+conj(INexpr(iter.row(), iter.col())));
		}else{
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}

//Full=Symmetric
template<class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const SymmetricMatrix *ins)
{
	typename SymmetricMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()!=iter.col()){
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
			in.put(iter.col(), iter.row(), in.get(iter.col(), iter.row())+INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}


//Full=Identity
template<class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const IdentityMatrix *ins)
{
	typename IdentityMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
		++iter;
	}
}

//Full=Diagonal
template<class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const DiagonalMatrix *ins)
{
	typename DiagonalMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}

//Hermitian=Diagonal
template<class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr, const HermitianMatrix *outs, const DiagonalMatrix *ins)
{
	typename DiagonalMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
		++iter;
	}
}

//Hermitian=Identity
template<class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr,const HermitianMatrix *outs, const IdentityMatrix *ins)
{
	typename IdentityMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
		++iter;
	}
}


//Symmetric=Identity
template<class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr,const SymmetricMatrix *outs, const IdentityMatrix *ins)
{
	typename IdentityMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
		++iter;
	}
}

//Symmetric=Diagonal
template<class mat, class expr>
inline void AddAssignMat(const mat &in,const expr &INexpr,const SymmetricMatrix *outs, const DiagonalMatrix *ins)
{
	typename DiagonalMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())+INexpr(iter.row(), iter.col()));
		++iter;
	}
}



template<class T, class  INstructure>
template<class nt1, class st1>
inline void _matrix<T, INstructure>::operator -= (const _matrix<nt1, st1> &m2)
{
	typedef structure outstruct ;
	typedef typename _matrix<nt1, st1>::structure instruct;

	static typename MatOutType(instruct, outstruct) s1;
	static instruct intype;

#ifdef MatAssignCheck
	if(!CanAssign(type(),s1.type())) structureError(MatName(s1.type()));	//check to make sure structures are compatable
#endif

	SubAssignMat(*this, m2, &s1, &intype);
//	return *this;
}

template<class T, class  INstructure>
template<class expr>
inline void _matrix<T, INstructure>::operator -= (const _matrixExpr<expr> &m2)
{
	typedef structure outstruct ;
	typedef typename _matrixExpr<expr>::structure instruct;
	static typename MatOutType(instruct, outstruct) s1;
	static typename _matrixExpr<expr>::structure intype;
#ifdef MatAssignCheck
	if(!CanAssign(type(),s1.type())) structureError(MatName(s1.type()));	//check to make sure structures are compatable
#endif

	SubAssignMat(*this, m2, &s1, &intype);
//	return *this;
}



//for "+=" operations

//the default one
//Full-=Full uses default
//Hermitian-=Hermiaitan uses default
//Symmetric-=Symmetric uses default
//Diagonal-=Diagonal uses default
//Identity-=Identity uses default (not really legal)
//Diagonal-=Identity uses default
//Hermitian-=Symmetric uses default

template<class outStructure, class inStructure, class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr, const outStructure *outst, const inStructure *inst)
{
	typename outStructure::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(),in.get(iter.row(),iter.col())-INexpr(iter.row(), iter.col()));
		++iter;
	}
}



//Full=Hermitian
template<class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const HermitianMatrix *ins)
{
	typename HermitianMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()!=iter.col()){
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
			in.put(iter.col(), iter.row(), in.get(iter.row(), iter.col())-conj(INexpr(iter.row(), iter.col())));
		}else{
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}

//Full=Symmetric
template<class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const SymmetricMatrix *ins)
{
	typename SymmetricMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()!=iter.col()){
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
			in.put(iter.col(), iter.row(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		}else{
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}


//Full=Identity
template<class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const IdentityMatrix *ins)
{
	typename IdentityMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		++iter;
	}
}

//Full=Diagonal
template<class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr, const FullMatrix *outs, const DiagonalMatrix *ins)
{
	typename DiagonalMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		if(iter.row()==iter.col()){
			in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		}
		++iter;
	}
}

//Hermitian=Diagonal
template<class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr, const HermitianMatrix *outs, const DiagonalMatrix *ins)
{
	typename DiagonalMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		++iter;
	}
}

//Hermitian=Identity
template<class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr,const HermitianMatrix *outs, const IdentityMatrix *ins)
{
	typename IdentityMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		++iter;
	}
}


//Symmetric=Identity
template<class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr,const SymmetricMatrix *outs, const IdentityMatrix *ins)
{
	typename IdentityMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		++iter;
	}
}

//Symmetric=Diagonal
template<class mat, class expr>
inline void SubAssignMat(const mat &in,const expr &INexpr,const SymmetricMatrix *outs, const DiagonalMatrix *ins)
{
	typename DiagonalMatrix::iterator iter(in.rows(), in.cols());
	while (iter){
		in.put(iter.row(), iter.col(), in.get(iter.row(), iter.col())-INexpr(iter.row(), iter.col()));
		++iter;
	}
}



/*
template<int rows, int columns, int I, int J>
class mat_Assign2 {
public:
    enum { go = (J < columns - 1) ? 1 : 0 };

    template<class matin, class expr, class upd>
    static inline void assign(matin& mat, expr expr, upd u)
    {
        u.update(mat(I,J), expr(I,J));
        mat_Assign2<rows * go, columns * go, I * go, (J+1) * go>
            ::assign(mat, expr, u);
    }
};

template<>
class mat_Assign2<0,0,0,0> {
public:
    template<class matin, class expr, class updater>
    static inline void assign(matin& mat, expr expr_, upd u)
    { }
};

template<int rows, int columns, int I>
class mat_Assign {
public:
    enum { go = (I < N_rows-1) ? 1 : 0 };

    template<class matin, class expr, class upd>
    static inline void assign(matin& mat, expr expr_, upd u)
    {
        mat_Assign2<rows, columns, I, 0>::assign(mat, expr_, u);
        mat_Assign<rows * go, columns * go, (I+1) * go>
            ::assign(mat, expr, u);
    }
};

template<>
class mat_Assign<0,0,0> {
public:
    template<class matin, class expr, class upd>
    static inline void assign(matin& mat, expr expr_, upd u)
    { }
};
*/

END_BL_NAMESPACE


#endif
