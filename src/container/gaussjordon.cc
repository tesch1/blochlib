/* guass-Jordon elimination

	performs a gauss-jordan elimination

	solves this set of linear equations

	A*x_i*Y=b_i  where A and Y is a matrix, x_i and b_i are vectors

	It also solves for the inverse of A-->A^-1 ==Y

	the below function simply takes in the A matrix (NxN) and
	a matrix consisting of all the x_i (NxM) vectors to be solved for
	A is then replaced with  A^-1 and x_i is replaced with b_i
*/

#ifndef _gauss_jordon_cc_
#define _gauss_jordon_cc_

#include<mat.h>
#include<iostream>
#include<string>

BEGIN_BL_NAMESPACE


//first overload...only 1 b to solve for and thus is a vector

template< class Any > template< class T>
int mat<Any>::GaussJordon(mat<Any> &a, vector<T> &b){
	int i, j, len=b.size();
	mat<Any> tm(len,1);
	for(i=0;i<len;i++)	tm[i][0]=b[i];
	if(GaussJordon(a, tm)==-1){
		return -1;
	}

	for(i=0;i<len;i++)	b[i]=tm[i][0];
	return 1;
}


template< class Any > template< class T >
int mat<Any>::GaussJordon(mat<Any> &a, mat<T> &b){
	//test squareness
	if(!a.issquare()){
		cerr<<"GaussJordon error::"<<endl;
		cerr<<" 'a' matrix MUST be square"<<endl;
		return -1;
	}

	//test that b.rows()==a.rows()
	int n=a.rows(), i,j,k;
	int m=b.cols();
	int irow=0, icol=0;

	if(b.rows()!=n){
		cerr<<"GaussJordon error::"<<endl;
		cerr<<" dimensions must be A(NxN) b(NxM)"<<endl;
		return -1;
	}
	vector<int> ixc(n,0);
	vector<int> ixr(n,0);
	vector<int> piv(n,0);

	Any bb, pinv, tt;

	for(i=0;i<n;i++){
		bb=T(0.);
		for(j=0;j<n;j++){
			if(piv[j]!=1){
				for(k=0;k<n;k++){
					if(piv[k]==0){
						if(abs(a[j][k])>=bb){
							bb=abs(a[j][k]);
							irow=j;
							icol=k;
						}
					}else if(piv[k]>1){
						cerr<<"Error: GaussJordon::"<<endl;
						cerr<<" singular matrix 'a'..."<<endl;
						return -1;
					}
				}
			}
		}
		++piv[icol];

		if(irow!=icol){
			for(j=0;j<n;j++) swap(a[irow][j], a[icol][j]);
			for(j=0;j<m;j++) swap(b[irow][j], b[icol][j]);
		}
		ixr[i]=irow;
		ixc[i]=icol;

		if(a[icol][icol]==0.){
			cerr<<"Error:: GaussJordon"<<endl;
			cerr<<" singular matrix 'a'..."<<endl;
			return -1;
		}

		pinv=T(1.)/(a[icol][icol]);
		a[icol][icol]=T(1.);
		for(j=0;j<n;j++) a[icol][j] *= pinv;
		for(j=0;j<m;j++) b[icol][j] *= pinv;

		for(j=0;j<n;j++){
			if(j!=icol){
				tt=a[j][icol];
				a[j][icol]=T(0.);

				for(k=0;k<n;k++) a[j][k] -= a[icol][k]*tt;
				for(k=0;k<m;k++) b[j][k] -= b[icol][k]*tt;
			}
		}
	}

	for(i=n-1;i>=0;i--){
		if(ixr[i]!=ixc[i]){
			for(k=0;k<n;k++){
				swap(a[k][ixr[i]],a[k][ixc[i]]);
			}
		}
	}
	return 1;
}


template<class Any> template<class T>
void mat<Any>::GaussJordon(vector<T> &b){
	if(GaussJordon(*this, b)==-1){
		cerr<<endl<<"Error: mat<Any>::GaussJordon(vector<T>)"<<endl;
		cerr<<" ooch, the slate feels cold to the unwanted...."<<endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}
}



template<class Any> template<class T>
void mat<Any>::GaussJordon(mat<T> &b){
	if(GaussJordon(*this, b)==-1){
		cerr<<endl<<"Error: mat<Any>::GaussJordon(vector<T>)"<<endl;
		cerr<<" ooch, the slate feels cold to the unwanted...."<<endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}
}

END_BL_NAMESPACE



#endif



