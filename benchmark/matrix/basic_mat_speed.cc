/*
  This is the matrix multiply driver

  The output format: "Size: %u\tmflop/s: %g\n"
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <float.h>
#include <time.h>
#include <math.h>

#include <sys/types.h>

#include "blochlib.h"

using namespace BlochLib;
using namespace std;


timer stopwatch;
/*
  We try to run enough iterations to get reasonable timings.  The matrices
  are multiplied at least MIN_RUNS times.  If that doesn't take MIN_CPU_SECS
  seconds, then we double the number of iterations and try again.

  You may want to modify these to speed debugging...
*/
#define MIN_RUNS     4
#define MIN_CPU_SECS 1

/*
  Note the strange sizes...  You'll see some interesting effects
  around some of the powers-of-two.
*/
#ifndef FASTTEST
/*
  Change that to zero for an abbreviated list while debugging...
*/
const int test_sizes[] = {
     2,
     5,
     20,
     24,
     31,
     32,
     36,
     48,
     64,
     73,
     96,
     97,
     127,
     128,
     129,
     163,
     191,
     192,
     229,
     255,
     256,
     257,
     319,
     320,
     321,
     417,
     479,
     480,
     511,
     512,
};

# define MAX_SIZE 512u
#else
const int test_sizes[] = {
     20,
     24,
     31,
     32,
     36,
     48,
     127,
     128,
     129,
     255,
     256,
     257,
     319,
     320,
     321,
     511,
     512,
};
#  define MAX_SIZE 512u
#endif

#define TTYPE Complex<double>

#define N_SIZES ((int) sizeof (test_sizes) / sizeof (int))


template<class Mat_T>
Mat_T mulmat(Mat_T &a, Mat_T &b)
{
 int i,j,k;
 Mat_T c(a.rows(), b.cols());
 for(i=0;i<c.rows();++i){
   for(j=0;j<c.cols();++j){
     c(i,j)=a(i,0) * b(0,j);
     for(k=1; k<c.cols();++k){
       c(i,j)+=a(i,k) * b(k,j);
     }
   }
 }
 return c;
}


template<class Mat_T>
double
time_mul1 ( Mat_T &A,  Mat_T &B, Mat_T& C)
{
    clock_t cpu_time;
    clock_t last_clock;
    double mflops, mflop_s;
    double secs = -1;

    int num_iterations = MIN_RUNS;
    int i, M=A.rows();

    while (secs < MIN_CPU_SECS) {

        cpu_time = 0;
		last_clock = clock();
        for (i = 0; i < num_iterations; ++i) {
  			C=mulmat(A,B);
        }
        cpu_time += clock() - last_clock;

	//1 matrix multiplies (M^3)...
	// 1 assign (M)
	// if TTYPE is a double the size is 2x that of a float
        mflops  = (sizeof(TTYPE)/sizeof(float) * num_iterations * (M * M * M))/ 1.0e6;
        secs    = cpu_time / ((double) CLOCKS_PER_SEC);
        mflop_s = mflops/secs;

        num_iterations *= 2;
    }
    printf ("Size: %u\ttime: %g us\n", M, (2*secs/num_iterations)*1e6);

    return mflop_s;
}




int
main (void)
{
	int sz_i;
	double mflop_s;
	_matrix<TTYPE, FullMatrix> A,B,C;
	cout<<"Performing speed tests for mat muls with an element byte size of "<<sizeof(TTYPE)<<endl;
	Random<UniformRandom<TTYPE > > mrY(TTYPE(1.00), TTYPE(2.0));
	Vector<double>  basic_mul(N_SIZES);


    cout<<"\nPerforming Multipication Tests C=A*B"<<endl;
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);
		mflop_s = time_mul1(A, B, C);
		basic_mul[sz_i]=mflop_s;

		printf ("Size: %u\tmflop/s: %g\n", M, mflop_s);
		//validate_dgemm (A, B, C);
	}

	ofstream oo("basic_matspeed.m");
	oo<<"Infinity=Inf;"<<endl;
	oo<<"R=["<<test_sizes<<" "<<"];"<<endl;

	oo<<"basic_mul=["<<basic_mul<<" "<<"];"<<endl;

	oo<<"semilogx(R, basic_mul,'r');"<<endl;
	oo<<"title('propagator(c=a*b*adjoint(a)) 700 MHz Pentium III Xeon');"<<endl;
	oo<<"xlabel('NxN');"<<endl;
	oo<<"ylabel('MFLOPS');"<<endl;
	oo<<"legend('Basic C=A*B');"<<endl;


     return 0;
}

