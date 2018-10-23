
/* vector_speed.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-26-01
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

/*
	vector_speed.h--> performs speed tests for 'Vector<>' class

	to compile you must set "-ftemplate-depth-NN " where 'NN' must be larger then the
	Coord length...'N'

*/


#include "blochlib.h"

#include <valarray>



using namespace BlochLib;
using namespace std;

#if 1
 #define fdaxpy   fdaxpy_
 #define fidaxpy  fidaxpy_
 #define fidaxpyo fidaxpyo_
#endif

#if 0
 #define fdaxpy   FDAXPY
 #define fidaxpy  FIDAXPY
 #define fidaxpyo FIDAXPYO
#endif

extern "C" {
  void fdaxpy(const int& N, const double& da, double* x,
    const int& xstride, const double* y, const int& ystride);

  void fidaxpy(const double& a, double* x, const double* y,
    const int& length, const int& iters);

  void fidaxpyo(const double& a, double* x, const double* y,
    const int& length, const int& iters);
}


timer stopwatch;
void printTime(int nrounds, std::string mess=""){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds)*1.e6 << " microseconds for "<<mess<<"\n";
}

void Info(std::string mess)
{
	cout<<mess;
	cout.flush();
}

class myV{
	public:
		int size;
		double *d;

		myV(int i):
		 size(i)
		{
			d= new double[i];
		}

		~myV(){ delete [] d;	}

		myV &operator+=(const myV &in)
		{
			for(int i=0;i<size;++i){ d[i]+=in.d[i];	}
			return *this;
		}

		myV &operator=(const myV &in)
		{
			for(int i=0;i<size;++i){ d[i]=in.d[i];	}
			return *this;
		}
		double &operator()(int i)
		{
			return d[i];
		}
		double &operator[](int i)
		{
			return d[i];
		}
};

myV operator*(myV &b, double in)
{
	myV a(b.size);
	for(int i=0;i<a.size;++i){ a.d[i]=b.d[i]*in;	}
	return a;
}


template<class T>
void dummyf(T &oo){ }

void dummyf(double *oo){ }

void Cdummy(Vector<double> &a, Vector<double> &b, double con, int N)
{
	dummyf(a);
}

void Cdummy(myV &a, myV &b, double con, int N)
{
	dummyf(a);
}
void Cdummy(valarray<double> &a, valarray<double>& b, double con, int N)
{
		dummyf(a);
}

void Cdummy(double *a, double* b, double con, int N)
{
	dummyf(a);
}

void cdax(double *a, double *b, double con, int N){
	double cp=-con;
	for(int i=0;i<N;++i){
		a[i]+=b[i]*con;
	}
	for(int i=0;i<N;++i){
		a[i]+=b[i]*cp;
	}
}

void myVdaxpy(myV &a, myV &b, double con, int N=0){
	double cp=-con;
	a+=b*con;
	a+=b*cp;
}

template<class T>
void vecdax(Vector<T> &a, Vector<T> &b, T con, int N=0){
	double cp=-con;
	a+=b*con;
	a+=b*cp;
}

template<class T>
void valdax(valarray<T> &a, valarray<T> &b, T con, int N=0){
	double cp=-con;
	a+=b*con;
	a+=b*cp;
}

void blasdax(double *a, double *b, double con, int N){
	double cp=-con;
	const int st=1;
	fdaxpy(N,con, a,st,b, st);
	fdaxpy(N,cp, a,st,b, st);
}

#define LOOPER(doto, dummy, NOPs, mess, VECT) \
	stopwatch.reset();	\
	for(int i=0;i<loops;++i){ \
		doto \
	}	\
	calctime=stopwatch(); \
 	stopwatch.reset(); 	\
 	for(int i=0;i<loops;++i){	\
   dummy \
  } \
  ltime=stopwatch();	\
	tottime=calctime-ltime;	\
	cout<<"** "<<mess<<" **"<<std::endl;	\
	cout<<"time taken: "<<tottime*1.e6/double(loops)<<" microseconds"<<std::endl;	\
	cout<<"MFLOPS: "<<NOPs/(tottime/double(loops))/1.e6<<std::endl;	\
	VECT.push_back(NOPs/(tottime/double(loops))/1.e6); \

#define F77LOOPER(doto, dummy, NOPs, mess, VECT) \
	stopwatch.reset();	\
	doto \
	calctime=stopwatch(); \
 	stopwatch.reset(); 	\
 dummy \
  ltime=stopwatch();	\
	tottime=calctime-ltime;	\
	cout<<"** "<<mess<<" **"<<std::endl;	\
	cout<<"time taken: "<<tottime*1.e6/double(loops)<<" microseconds"<<std::endl;	\
	cout<<"MFLOPS: "<<NOPs/(tottime/double(loops))/1.e6<<std::endl;	\
	VECT.push_back(NOPs/(tottime/double(loops))/1.e6); \

int main(int argc, char **argv){

	int loops=0;
	int q=1;
	int Nmax, whx;
	query_parameter(argc, argv, q++, "Enter number of loops for each test: ", loops);
	query_parameter(argc, argv, q++, "Enter max exponent for vector length (10^N): ", Nmax);
	query_parameter(argc, argv, q++, "F77 Blas[0], F77 'optim'[1], valarray[2], 'C' array[3], Vector[4], all[5]: ", whx);
	int oldloops=loops;
		Vector<double> cdaxpyF, daxpyF,f77daxpyF, valdaxpyF, blasdaxpyF,myVdaxpyF;


	double calctime, ltime, tottime;
	Vector<int> logiter(Nmax*3,0);
	int ct=1;
	for(int i=0;i<logiter.size();i+=3){
		 logiter[i]=ct;
		 logiter[i+1]=ct*2;
		 logiter[i+2]=ct*5;
		 ct*=10;
	}


if(whx==1 || whx==5){
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);
		std::string presstr= " precision: double";
		double* x = new double[N];
  	double* y = new double[N];
		double con=Rand();
		for(int i=0;i<N;++i){	x[i]=Rand();y[i]=Rand();	}
		loops=oldloops/N;
		if(loops<2){	loops=2;	}
		F77LOOPER(fidaxpy(con, x, y, N, loops);,fidaxpyo(con, x, y, N, loops);, 2*N*4,  "a=a+const*b..."+lenstr+presstr, f77daxpyF)
		delete [] x;
		delete [] y;
	}
}
if(whx==0 || whx==5){
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);
		std::string presstr= " precision: double";
		double* x = new double[N];
  	double* y = new double[N];
		double con=Rand();
		for(int i=0;i<N;++i){	x[i]=Rand();y[i]=Rand();	}
	  int st=1;
		loops=oldloops/N;
		if(loops<2){	loops=2;	}
		LOOPER(blasdax(x,y,con,N);, fidaxpyo(con, x, y, N, loops);,4*N*2,  "a=a+const*b..."+lenstr+presstr, blasdaxpyF)
		delete [] x;
		delete [] y;
	}
}
if(whx==2 || whx==5){
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);
		std::string presstr= " precision: double";
		valarray<double> x(N), y(N);
		double con=Rand();
		for(int i=0;i<N;++i){	x[i]=Rand();y[i]=Rand();	}
	  int st=1;
	  loops=oldloops/N;
		if(loops<2){	loops=2;	}
		LOOPER(valdax(x,y,con);, Cdummy(x,y,con,N); , 4*N*2,  "a+=const*b..."+lenstr+presstr, valdaxpyF)
	}
}
if(whx==3 || whx==5){
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);
		std::string presstr= " precision: double";
		Vector<double> a(N), b(N);
		double con=Rand();
		loops=oldloops/N;
		if(loops<2){	loops=2;	}
		for(int i=0;i<N;++i){	a(i)=Rand();b(i)=Rand();	}
		LOOPER(vecdax(a,b,con);, Cdummy(a,b,con,N); , 4*N*2,  "a+=const*b..."+lenstr+presstr, daxpyF)
	}
}

if(whx==4 || whx==5){
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);
		std::string presstr= " precision: double";
		double *a=new double[N];
		double *b=new double[N];
		double con=Rand();
		loops=oldloops/N;
		if(loops<2){	loops=2;	}
		for(int i=0;i<N;++i){	a[i]=Rand();b[i]=Rand();	}
		LOOPER(cdax(a,b,con,N) ;, Cdummy(a,b,con,N); , 4*N*2,  "a+=const*b..."+lenstr+presstr, cdaxpyF)
	}
}

if(whx==5){
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);
		std::string presstr= " precision: double";
		myV a(N);
		myV b(N);
		double con=Rand();
		loops=oldloops/N;
		if(loops<2){	loops=2;	}
		for(int i=0;i<N;++i){	a[i]=Rand();b[i]=Rand();	}
		LOOPER(myVdaxpy(a,b,con,N) ;, Cdummy(a,b,con,N); , 4*N*2,  "a+=const*b..."+lenstr+presstr, myVdaxpyF)
	}
}


	ofstream oo("vectorspeed.m");
	oo<<"R=["<<logiter<<"];"<<endl;
	oo<<"daxpyF=["<<daxpyF<<"];"<<endl;
	oo<<"f77daxpyF=["<<f77daxpyF<<"];"<<endl;
	oo<<"blasdaxpyF=["<<blasdaxpyF<<"];"<<endl;
	oo<<"valdaxpyF=["<<valdaxpyF<<"];"<<endl;
	oo<<"cdaxpyF=["<<cdaxpyF<<"];"<<endl;
	oo<<"myVdaxpyF=["<<myVdaxpyF<<"];"<<endl;
	oo<<"semilogx(R, daxpyF, 'k', R, f77daxpyF, 'b', R, blasdaxpyF, 'r', R, valdaxpyF, 'g',R, cdaxpyF, 'k.', R, myVdaxpyF);"<<endl;
	oo<<"title('DAXPY(a+=const*b) 700 MHz Pentium III Xeon');"<<endl;
	oo<<"xlabel('log(N)');"<<endl;
	oo<<"ylabel('MFLOPS');"<<endl;
	oo<<"legend('Vector', 'f77', 'BLAS', 'valarray', 'C array', 'vec class');"<<endl;




}

