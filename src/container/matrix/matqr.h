
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

/**************************************************
 matqr.h


 does all the QR realted transforms

 should be included from '_matrix.h'

 the 'external' functions can be used, but it is probably
 better to use the 'Internal' function as all the 'external' functions
 as you may need to use the 'std::' namespace call for them in your program

 **************************************************/


#ifndef _matQR_h_
#define _matQR_h_



BEGIN_BL_NAMESPACE



#define DoubleSign(a,b) ((b) >= 0.0 ? abs(a) : -abs(a))

/*************************BEGIN QR EXTERNAL TO MATRIX*******************/

//function External to Matrix class
// QR decomposition... The PACKED version
//THIS IS DETRUCTIVE to A
// it returns '1' is sigularity is found
// '0' if not
//
// THis is the 'packed' version of the QR...the
// A packed HouseHolder/R matrix where the upper-tri (and diag) is the 'R' matrix
// and the lower-tri is a packed Householder form of Q...
//
// UpperTri+diag of A will have the 'R' matrix (rest will be 0)
// lower diag of W return the Q householder non zero elements on the off diagonal and the
// 'scaling' coiefs on the diagonal

template<class T, class INstructure>
int PackedQR(_matrix<T, INstructure> &a, _matrix<T, FullMatrix> &W)
{

	/* %Matlab Equivilent
	%except here the W(j,j)=coief(j)
	%and W(i,j)=W(i,j) (for the lower diag)
	%and R(j,i)=a(j,i) (for upper diag and diag)

	[m,n] = size(A);
	R = A;
	W = zeros(m,n);
	for j = 1:n
	  W(j:m,j) = R(j:m,j);
	  if W(j,j) ~= 0
	    z = -sign(W(j,j));
	  else
	    z = 1;
	  end
	  R(j,j) = z*norm(W(j:m,j));
	  R(j+1:m,j) = 0;
	  if R(j,j) ~= 0
	     W(j,j) = W(j,j) - R(j,j);
	     W(j:m,j) = W(j:m,j)*(1/sqrt(2*abs(R(j,j)*W(j,j))));
	     R(j:m,j+1:n) = R(j:m,j+1:n) - 2*W(j:m,j)*(W(j:m,j)'*R(j:m,j+1:n));
  	end
	*/

	int n=a.rows();
	int m=a.cols();
	int i,j,k, st=0;
	T z,sum;
	int sing=0;

	W.resize(n,m,0);

	for(j=0; j<n;j++)
	{
		sum=ZeroType<T>::zero();
		for(i=j;i<m;i++){
			W(i,j)=a(i,j);
			sum+=W(i,j)*W(i,j);
		}

		a(j,j)=-DoubleSign(sqrt(sum),W(j,j));
		for(i=j+1;i<m;i++) a(i,j)=0;
		if(a(j,j)!=ZeroType<T>::zero())
		{
			W(j,j)-=a(j,j);
			T scale=1./(sqrt(2.0*abs(a(j,j)*W(j,j))));
			for(i=j;i<m;i++)	W(i,j)*=scale;

			if(j<n-1)
			{
				Range rr(j,m-1), cc(j+1, n-1);
				//this is R(j:m,j+1:n) = R(j:m,j+1:n) - 2*W(j:m,j)*(W(j:m,j)'*R(j:m,j+1:n));
				a.put(rr,cc, a(rr,cc)-2.0*cross(W(rr, j),W(rr, j)*a(rr, cc)));
			}
		}
	}
	return 1;
}

//simple case for the diagonal and identity types
template<class T>
int PackedQR(_matrix<T, DiagonalMatrix> &a, _matrix<T, FullMatrix> &W)
{
	coief.resize(a.rows());
	coief.fill(a.rows());
	return 0;
}


//simple case for the diagonal and identity types
template<class T>
int PackedQR(_matrix<T, IdentityMatrix> &a, _matrix<T, FullMatrix> &W)
{
	coief.resize(a.rows());
	coief.fill(a.rows());
	return 0;
}

//the symmetric and Hermitian matrix types cannot be
//performed yet...they MUST be full matrices to use the 'PackedQR'
template<class T>
int PackedQR(_matrix<T, HermitianMatrix> &a, _matrix<T, FullMatrix> &W)
{
	std::cerr<<"\n**Hermitian Matrix Types CANNOT have PackedQR() performed\n  ..convert to full matrix first"<<std::endl;
	return 0;
}

//the symmetric and Hermitian matrix types cannot be
//performed yet...they MUST be full matrices to use the 'PackedQR'
template<class T>
int PackedQR(_matrix<T, SymmetricMatrix> &a, _matrix<T, FullMatrix> &W)
{
	std::cerr<<"\n**Symmetric Matrix Types CANNOT have PackedQR() performed\n  ..convert to full matrix first"<<std::endl;
	return 0;
}


// UNpacks the Packed form of the QR algo above...
// W--> the matrix that contains the 'packed' form of the Q matrix (from PackedQR)
// Q--> the output Q matrix

/*
function Q = formQ(W);
%FORMQ    Compute Q after house(A)
%  Q = FORMQ(W) returns the matrix Q from the QR decomposition
%  of a matrix A computed by [W,R] = HOUSE(A).

[m,n] = size(W);
Q = eye(m,m);

for j = n:-1:1
  if W(j,j) ~= 0
     Q(j:m,j:m) = Q(j:m,j:m) - 2*W(j:m,j)*(W(j:m,j)'*Q(j:m,j:m));
  end
end
*/

template<class T>
void UnPackQR(_matrix<T, FullMatrix> &W,  _matrix<T, FullMatrix> &Q)
{
	int n=W.rows(), i;
	Q.identity(n);
	for(i=n;i--;)
	{
		if(W(i,i)!=0)
		{
			Range R(i, n-1), C(i, n-1);
			Q.put(R,C, Q(R,C)-2.0*cross(W(R,i),W(R,i)*Q(R,C)));
		}
	}
}

//EXTERNAL to Matrix class...NOT sdestructive to the input matrix
// QR decomposition...
// The R matrix is stored upper-trangular (with components on the diagonal)
// the Q matrix is stored in 'Q'--> is orthronormal matrix
// (*this)=QR

//default
template<class T, class INstructure>
void QR(_matrix<T, INstructure> &A, _matrix<T, FullMatrix> &Q, _matrix<T, FullMatrix> &R)
{
	//our workspace
	R=A;
	_matrix<T, FullMatrix> W;
	std::PackedQR(R,W);
	std::UnPackQR(W,Q);
}

/******************************END EXTERNAL QR'S************************/

/*******************************BEGIN INTERANAL MATRIX QR FUNCTIONS ***/

//most of these simply call the 'external' functions...as they are typically
//how people call these functions themselves

//internal to Matrix class...DESTRCUTIVE to the input matrix
// PackedQR decomposition...
// The R matrix is stored upper-trangular (with components on the diagonal)
// The W matrix conatins the 'packed' Housholder Q matrix

//default
template<class T, class INstructure>
void _matrix<T, INstructure>::PackedQR(_matrix<T, FullMatrix> &W)
{
	std::PackedQR((*this),W);
}

// UNpacks the Packed form of the QR algo above...
// W--> the matrix that contains the 'packed' form of the Q matrix (from PackedQR)
// Q--> the output Q matrix

template<class T, class INstructure>
void _matrix<T, INstructure>::UnPackQR(_matrix<T, FullMatrix> &W,_matrix<T, FullMatrix> &Q)
{
	std::UnPackQR(W,Q);
}

//internal to Matrix class...NOT sdestructive to the input matrix
// QR decomposition...
// The R matrix is stored upper-trangular (with components on the diagonal)
// the Q matrix is stored in 'Q'--> is orthronormal matrix
// (*this)=QR

//default
template<class T, class INstructure>
void _matrix<T, INstructure>::QR(_matrix<T, FullMatrix> &Q, _matrix<T, FullMatrix> &R)
{
	//our workspace
	R=(*this);
	_matrix<T, FullMatrix> W;
	std::PackedQR(R,W);
	std::UnPackQR(W,Q);
}

/***********************END INTERNAL QR MATRIX FUNCTIONS***************/


/************************BEING EXTERNAL GRAM-SCHMIT ORHTNORMALIZATION *********/


//this is the main recursive function for the ortho process
template<class T, class INstructure>
void GramSchmidtRecurNorm( _matrix<T, INstructure> &y,
					   _matrix<T, INstructure> &z,
					   Vector<double> &norms,
					  int curpos)
{
	if(curpos==y.cols()) return;
	if(curpos==0)
	{
		norms(0)=norm(z.col(0));
		if(!norms(0))
		{
			std::cerr<<" Warning: GramSchmidtRecurNorm"<<std::endl;
			std::cerr<<" Signular matrix...cannot continue..."<<std::endl;
			return;
		}
		y.putCol(0, z.col(0)/norms(0));
	//	cout<<"Norms: "<<norms(0)<<endl<<"col: "<<z.col(0)<<" ycol: "<<y.col(0)<<endl;
		GramSchmidtRecurNorm(y,z,norms, curpos+1);
	}else{
		int i=0;
		Vector<T> tz=z.col(curpos);
		static Vector<T> ty(y.rows());
		while(i<curpos)
		{
			ty=y.col(i);
			if(!norms(i))
			{
				std::cerr<<" Warning: GramSchmidtRecurNorm"<<std::endl;
				std::cerr<<" Signular matrix...cannot continue..."<<std::endl;
				return;
			}
			tz-=dot(ty,tz)*ty;
			i++;

		}
		norms(curpos)=norm(tz);
		y.putCol(curpos,tz/norms(curpos));
	//	norms(curpos)=norm(tz);
	//	cout<<" curpose: "<<curpos<<" Norms: "<<norms(curpos)<<endl<<"Y: "<<y.col(curpos)<<endl<<" tz: "<<tz<<endl;

		GramSchmidtRecurNorm(y,z,norms, curpos+1);
	}
}

template<class T, class INstructure>
_matrix<T, INstructure> GramSchmidt( _matrix<T, INstructure> &z)
{
	if(!z.issquare())
	{
		std::cerr<<std::endl<<" Error: GramSchmidt()"<<std::endl;
		std::cerr<<" to orthnormalize a basis, the matrix must be square..."<<std::endl;
		std::cerr<<" returning original matrix...."<<std::endl;
		return z;
	}
	_matrix<T, INstructure> y(z.rows(), z.cols(),0);
	Vector<double> norms(z.rows(),0);
	GramSchmidtRecurNorm(y,z,norms,0);
	return y;
}

template<class T, class INstructure>
_matrix<T, INstructure> GramSchmidt( _matrix<T, INstructure> &z, Vector<double> &norms)
{
	if(!z.issquare())
	{
		std::cerr<<std::endl<<" Error: GramSchmidt()"<<std::endl;
		std::cerr<<" to orthnormalize a basis, the matrix must be square..."<<std::endl;
		std::cerr<<" returning original matrix...."<<std::endl;
		return z;
	}
	_matrix<T, INstructure> y(z.rows(), z.cols(),0);
	norms.resize(z.rows(),0);
	GramSchmidtRecurNorm(y,z,norms,0);
	return y;
}

template<class T, class INstructure>
_matrix<T, INstructure> OrthoNorm( _matrix<T, INstructure> &z)
{
	if(!z.issquare())
	{
		std::cerr<<std::endl<<" Error: GramSchmidt()"<<std::endl;
		std::cerr<<" to orthnormalize a basis, the matrix must be square..."<<std::endl;
		std::cerr<<" returning original matrix...."<<std::endl;
		return z;
	}
	_matrix<T, INstructure> y(z.rows(), z.cols(),0);
	Vector<double> norms(z.rows(),0);
	GramSchmidtRecurNorm(y,z,norms,0);
	return y;
}

END_BL_NAMESPACE


#endif
