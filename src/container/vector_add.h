#ifndef _vector_add_h_
#define _vector_add_h_

#include<vector.h>
#include<complex.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<string>
#include<fstream>

BEGIN_BL_NAMESPACE


/* WARNING?NOTE?LOOKOUT?RUN?HIDE*****************
this little file ONLY adds a few functions to the exsisting STL vector container
it is NOT effecient by any means...but for quick and dirty things, it should be fine..

The super-fast vector class can be found in 'Vector.h' (capital 'V')
where the object has been expression-templated
*/
/*
 *here is set up a simple vector
 *it simply adds on to the exsisting vector routines...
 *i.e. it allows one to creat a vector from a ptr array
 *of number types (the size of the array must be known!!)
 *
 */

template<class Comparable>
vector<Comparable> put_array(const Comparable* in, const int &start,const int &stop){
	vector<Comparable> tmp;
	for(int i=start;i<stop;i++){tmp.push_back(in[i]);}
	return tmp;
}

template<class Comparable>
vector<Comparable> put_array(const Comparable* in, const int &stop){
	return put_array(in, 0, stop);
}

//VECTOR MULTIPLICATION**************

//return[i]=l[i]*r[i]
template<class T1, class T2>
void mul(vector<T1> &eq,const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: mul(vector,vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	eq=l;
	for(i=0;i<r1;i++){
		eq[i]*=T1(r[i]);
	}
}


//return[i]=l[i]*r[i]
template<class T1, class T2>
vector<T1> mul(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: mul(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]*=T1(r[i]);
	}
	return tmp;
}

//return[i]=l[i]*r[i]
template<class T1, class T2>
vector<T1> operator*(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: *(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]*=T1(r[i]);
	}
	return tmp;
}

//return[i]*=r[i]
template<class T1, class T2>
vector<T1> operator*=(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: *=(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]*=T1(r[i]);
	}
	return tmp;
}


//VECTOR ADDTION*******************

//eq[i]=l[i]+r[i]
template<class T1, class T2>
void add(vector<T1> &eq, const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: add(vector,vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	eq=l;
	for(i=0;i<r1;i++){
		eq[i]+=T1(r[i]);
	}
}

//return[i]=l[i]+r[i]
template<class T1, class T2>
vector<T1> add(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: add(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]+=T1(r[i]);
	}
	return tmp;
}

//return[i]=l[i]+r[i]
template<class T1, class T2>
vector<T1> operator+(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: +(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]+=r[i];
	}
	return tmp;
}

//return[i]+=r[i]
template<class T1, class T2>
vector<T1> operator+=(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: +=(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]+=r[i];
	}
	return tmp;
}


//VECTOR SUBTRACTION***************************

//eq[i]=l[i]-r[i]
template<class T1, class T2>
void sub(vector<T1> &eq, const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: sub(eq, vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	eq=l;
	for(i=0;i<r1;i++){
		eq[i]-=T1(r[i]);
	}
}

//return[i]=l[i]-r[i]
template<class T1, class T2>
vector<T1> sub(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: sub(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]-=T1(r[i]);
	}
	return tmp;
}

//return[i]=l[i]-r[i]
template<class T1, class T2>
vector<T1> operator-(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: -(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]-=T1(r[i]);
	}
	return tmp;
}

//return[i]-=r[i]
template<class T1, class T2>
vector<T1> operator-=(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: -=(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]-=T1(r[i]);
	}
	return tmp;
}

//VECTOR DIVISION**********************************
//eq[i]=l[i]/r[i]
template<class T1, class T2>
void div(vector<T1> &eq, const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: div(vector, vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	eq=l;
	for(i=0;i<r1;i++){
		eq[i]/=T1(r[i]);
	}
}

//return[i]=l[i]/r[i]
template<class T1, class T2>
vector<T1> div(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: div(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]/=T1(r[i]);
	}
	return tmp;
}

//return[i]=l[i]/r[i]
template<class T1, class T2>
vector<T1> operator/(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: /(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]/=T1(r[i]);
	}
	return tmp;
}

//return[i]/=r[i]
template<class T1, class T2>
vector<T1> operator/=(const vector<T1> &l, const vector<T2> &r){
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: /=(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	vector<T1> tmp=l;
	for(i=0;i<r1;i++){
		tmp[i]/=T1(r[i]);
	}
	return tmp;
}


//VECTOR "OTHER" FUNCTIONS*****************************
template<class T1>
T1 sum(vector<T1> &moo){		//SUM
	int i, len=moo.size();
	T1 tmp=T1(0.);
	for(i=0;i<len;i++){
		tmp+=moo[i];
	}
	return tmp;
}

template<class T1>
T1 prod(vector<T1> &moo){		//Product
	int i, len=moo.size();
	T1 tmp=T1(1);
	for(i=0;i<len;i++){
		tmp*=moo[i];
	}
	return tmp;
}

template<class T1>
T1 product(vector<T1> &moo){		//Product
	int i, len=moo.size();
	T1 tmp=T1(1);
	for(i=0;i<len;i++){
		tmp*=moo[i];
	}
	return tmp;
}

template<class T1, class T2>
T1 dot(vector<T1> &l, vector<T2> &r){			//DOT
	int r1=l.size(), r2=r2.size();
	if(r1!=r2){
		cerr<<"Error: dot(vector, vector)"<<endl;
		cerr<<" vectors must be the same length"<<endl;
		exit(0);
	}
	int i;
	T1 tmp=T1(0.);
	for(i=0;i<r1;i++){
		tmp+=l[i]*r[i];
	}
}

template<class T1, class T2>
vector<T1> ExpMul(vector<T1> &a, T2 &fact){		//does exp(-i fact)*a[i]
	int i, len=a.size();
	if(len==0) return a;
	vector<T1> tmp=a;
	for(i=0;i<len;i++){
		tmp[i]=exp(-float(i)*T2(fact))*a[i];
	}
	return tmp;
}


template<class T1>
vector<complex> FFT(vector<T1> &in){				//FFT OUTPUT WILL BE COMPEX!!!
	vector<complex> data;
	int mm;
	int size=in.size();
	if(size==0) return data;
	data.resize(size, complex(0.,0.));
	for(mm=0;mm<size;mm++){
		data[mm]=complex(in[mm]);
	}
	complex w, wp, temp;
	double theta, sto2;
	static const double Pi=3.14159265358979323846;
	int i,j,m,mmax,istep,ii;

	for(j=0,i=0; i<size; i++){
		if(j>i){
			temp = data[j];
			data[j] = data[i];
			data[i] = temp;
		}
	    m = size/2;
	    while((m>1) && (j>=m)){
			j = j-m;
			m = m/2;
		}
	    j=j+m;
	}
	mmax=1;
	while (size>mmax){
		istep = 2*mmax;
		theta = Pi/(mmax);
		sto2 = sin(theta/2.);
		wp = complex(-2*sto2*sto2, sin(theta));
		w=1;
		for(ii=0; ii<mmax; ii++){
			for(i=ii; i<size; i+=istep){
				j = i+mmax;			// Danielson Lanczos
				temp = w*data[j];
				data[j] = data[i] - temp;
				data[i] += temp;
			}
			w+=w*wp;
		}
		mmax *= 2;
	}
	for(i=0,j=size/2; j<size; i++,j++){
		temp = data[j];
		data[j] = data[i];
		data[i] = temp;
	}
	return data;
}


template< class T >
vector<complex> IFFT(vector<T> &in){  //inverse FFT OUTPUT IS COMPLEX VECTOR!!
	vector<complex> data;
	int mm;
	int size=in.size();
	if(size==0) return data;
	data.resize(size, complex(0.,0.));
	for(mm=0;mm<size;mm++){
		data[mm]=complex(in[mm]);
	}
	complex w, wp, temp;
	double theta, sto2;
	static const double Pi=3.14159265358979323846;
	int i,j,m,mmax,istep,ii;
	for(i=0,j=size/2; j<size; i++,j++){
		temp = data[j];
		data[j] = data[i];
		data[i] = temp;
	}
	for(j=0,i=0; i<size; i++){
		if(j>i){
			temp = data[j];
			data[j] = data[i];
			data[i] = temp;
		}
		m = size/2;
		while((m>1) && (j>=m)){
			j = j-m;
			m = m/2;
		}
		j=j+m;
	}
	mmax=1;
	while (size>mmax){
		istep = 2*mmax;
		theta = Pi/(-mmax);
		sto2 = sin(theta/2.);
		wp = complex(-2*sto2*sto2, sin(theta));
		w=1;
		for(ii=0; ii<mmax; ii++){
			for(i=ii; i<size; i+=istep){
				j = i+mmax;			//Danielson Lan
				temp = w*data[j];
				data[j] = data[i] - temp;
				data[i] += temp;
			}
			w+=w*wp;
		}
		mmax *= 2;
	}

	return data;
}

//VECTOR OUTPUT************************************
template<class T>
ostream& operator<<(ostream &otr,const vector<T> &outty){
	copy(outty.begin(), outty.end(),  ostream_iterator<T>(otr, " "));
	return otr;
}


//FORCED VECTOR TEMPLATE LENGTH******************
template<class T, int N>
class NVector{
	private:
		vector<T> data;
	public:
		NVector(){	data.resize(N); }
		NVector(T in){	data.resize(N, in);	}
		NVector(const NVector &in){ data=in.data;	}

		const T &operator[](int r)const { return data[r];	}
		T operator[](int r)	{ return data[r];	}
		const T &operator()(int r)const{	return data[r];	}
		T operator()(int r){	return data[r];	}
		NVector &operator=(const NVector &r){
			if(&r!=this) {
				data=r.data;
				return *this;
			}
		}

		NVector &operator=(const T r){
			data.resize(N);
			for(int i=0;i<N;i++){
				data[i]=r;
			}
		}

		int size()const { return N;	}

		typedef vector<T>::iterator iterator;
		typedef vector<T>::const_iterator const_iterator;

		iterator begin(){ return data.begin();	}
		const_iterator begin()const{ return data.begin();	}
		iterator end(){ return data.end();	}
		const_iterator end()const{ return data.begin();	}
};

template<class T, int N>
ostream& operator<<(ostream &otr,const NVector<T,N> &outty){
	for(int i=0;i<N;i++)
		otr<<outty[i]<<" ";
	return otr;
}

END_BL_NAMESPACE


#endif

