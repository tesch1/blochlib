
#include "blochlib.h"
#include <valarray>

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

#define F77LOOPER(len,  doto, NOPs, mess, VECT) \
	stopwatch.reset();	\
	doto \
	calctime=stopwatch(); \
 	stopwatch.reset(); 	\
  fidaxpyo(con,x,y,len,loops); 	\
  ltime=stopwatch();	\
	tottime=calctime-ltime;	\
	cout<<"** "<<mess<<" **"<<std::endl;	\
	cout<<"time taken: "<<tottime*1.e6/double(loops)<<" microseconds"<<std::endl;	\
	cout<<"MFLOPS: "<<NOPs/(tottime/double(loops))/1.e6<<std::endl;	\
	VECT.push_back(NOPs/(tottime/double(loops))/1.e6); \
	
#define F77BLASLOOPER(len,  doto, NOPs, mess, VECT) \
	stopwatch.reset();	\
	for(int i=0;i<loops;++i){	\
		doto \
	}	\
	calctime=stopwatch(); \
 	stopwatch.reset(); 	\
  fidaxpyo(con,x,y,len,loops); 	\
  ltime=stopwatch();	\
	tottime=calctime-ltime;	\
	cout<<"** "<<mess<<" **"<<std::endl;	\
	cout<<"time taken: "<<tottime*1.e6/double(loops)<<" microseconds"<<std::endl;	\
	cout<<"MFLOPS: "<<NOPs/(tottime/double(loops))/1.e6<<std::endl;	\
	VECT.push_back(NOPs/(tottime/double(loops))/1.e6); \

template<class T>
void dummyf(T &in){}

#define VALLOOPER(len,  doto, NOPs, mess, VECT) \
	stopwatch.reset();	\
	for(int i=0;i<loops;++i){	\
		doto \
	}	\
	calctime=stopwatch(); \
 	stopwatch.reset(); 	\
 for(int i=0;i<loops;++i){	\
		dummyf(x); \
	}	\
	ltime=stopwatch();	\
	tottime=calctime-ltime;	\
	cout<<"** "<<mess<<" **"<<std::endl;	\
	cout<<"time taken: "<<tottime*1.e6/double(loops)<<" microseconds"<<std::endl;	\
	cout<<"MFLOPS: "<<NOPs/(tottime/double(loops))/1.e6<<std::endl;	\
	VECT.push_back(NOPs/(tottime/double(loops))/1.e6); \
	

	
int main(int argc, char **argv){
	
	int loops=0;
	int q=1;
	int Nmax;
	int whx=0;
	query_parameter(argc, argv, q++, "Enter number of loops for each test: ", loops);
	int oldloops=loops;
	query_parameter(argc, argv, q++, "Enter max exponent for vector length (10^N): ", Nmax);
	query_parameter(argc, argv, q++, "F77 Blas[0] or F77 'optim'[1] or valarray[2] or all[3]: ", whx);
	Vector<double>  f77daxpyF, valdaxpyF, blasdaxpyF;
		
	double calctime, ltime, tottime;
	Vector<int> logiter(Nmax*3,0);
	int ct=1;
	for(int i=0;i<logiter.size();i+=3){
		 logiter[i]=ct;
		 logiter[i+1]=ct*2;
		 logiter[i+2]=ct*5;
		 ct*=10;
	}

if(whx==3 || whx==1){
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);	
		std::string presstr= " precision: double";	
		double* x = new double[N];
  	double* y = new double[N];
		double con=Rand();	
		for(int i=0;i<N;++i){	x[i]=Rand();y[i]=Rand();	}	
		double pars=pow(10.0, (N+1)/4.0);
	 loops=oldloops;
		F77LOOPER(N, 	fidaxpy(con, x, y, N, loops);, 2*N*4,  "a=a+const*b..."+lenstr+presstr, f77daxpyF)	
		delete [] x;
		delete [] y;
	}
}else if(whx==0 || whx==3){
	loops=oldloops;
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);	
		std::string presstr= " precision: double";	
		double* x = new double[N];
  	double* y = new double[N];
		double con=Rand();	
		for(int i=0;i<N;++i){	x[i]=Rand();y[i]=Rand();	}	
	  int st=1;
		double pars=pow(10.0, (N+1)/4.0);
		loops=oldloops;
		F77BLASLOOPER(N, 	fdaxpy(N,con, x,st,y, st);, 4*N*2,  "a=a+const*b..."+lenstr+presstr, blasdaxpyF)	
		delete [] x;
		delete [] y;
	}
}else if(whx==2 || whx==3){
	loops=oldloops;
	for(int i=0;i<logiter.size();++i){
		int N=logiter[i];
		std::string lenstr="length: "+itost(N);	
		std::string presstr= " precision: double";	
		valarray<double> x(N), y(N);
		double con=Rand();	
		for(int i=0;i<N;++i){	x[i]=Rand();y[i]=Rand();	}	
	  int st=1;
		double pars=pow(10.0, (N+1)/4.0);
		loops=oldloops;
		VALLOOPER(N, 	y+=con*x; , 4*N*2,  "a+=const*b..."+lenstr+presstr, valdaxpyF)	
	}
}


	ofstream oo("fortranspeed");
	oo<<"R=["<<logiter<<"]"<<endl;
	oo<<"f77daxpyF=["<<f77daxpyF<<"]"<<endl;
	oo<<"blasdaxpyF=["<<blasdaxpyF<<"]"<<endl;
	oo<<"valdaxpyF=["<<valdaxpyF<<"]"<<endl;



}
