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

#define TTYPE Complex<float>
#define ERRR FLT_EPSILON

#define N_SIZES ((int) sizeof (test_sizes) / sizeof (int))



/*
  Dot products satisfy the following error bound:
   float(sum a_i * b_i) = sum a_i * b_i * (1 + delta_i)
  where delta_i <= n * epsilon.  In order to check your rmatrix
  multiply, we compute each element in term and make sure that
  your product is within three times the given error bound.
  We make it three times because there are three sources of
  error:

   - the roundoff error in your multiply
   - the roundoff error in our multiply
   - the roundoff error in computing the error bound

  That last source of error is not so significant, but that's a
  story for another day.
 */

template<class Mat_T>
void
validate_dgemm (Mat_T &A, Mat_T &B, Mat_T &C)
{
    int i, j, k;

    typename Mat_T::numtype dotprod = 0,errorbound = 0;
	double err;
	for (i = 0; i < A.rows(); ++i) {
        for (j = 0; j < A.rows(); ++j) {
			dotprod = 0;
			errorbound = 0;
			for (k = 0; k < A.rows(); ++k) {
				dotprod += A(i,k) * B(k,j);
				errorbound += abs(A(i,k) * B(k,j));
			}
			errorbound *= (A.rows() * ERRR);

            err = abs(C(i,j) - dotprod);
            if (err > 3.0*abs(errorbound)) {
                printf("Matrix multiply failed.\n");
                printf("C(%d,%d) should be %g, was %g\n", i, j,
                       abs(C(i,j)), abs(dotprod));
                printf("Error of %g, acceptable limit %g\n",
                       err, 3.0*abs(errorbound));
                exit(-1);
            }
        }
    }
}

template<class Mat_T>
double
time_prop ( Mat_T &A,  Mat_T &B, Mat_T& C)
{
    clock_t cpu_time;
    clock_t last_clock;
    double mflops, mflop_s;
    double secs = -1;

    int num_iterations = MIN_RUNS;
    int i, M=A.rows(), did=0;

    while (secs < MIN_CPU_SECS) {

        cpu_time = 0;
		last_clock = clock();
        for (i = 0; i < num_iterations; ++i) {
  			C=prop(A,B);//C=A*B*adjoint(A);
  			did++;
        }
        cpu_time += clock() - last_clock;

	//2 matrix multiplies (M^3)...
	// What ever our 'element' is (TTYPE) is/size(float)
	// if TTYPE is a double the size is 2x that of a float
	// one adjoint operatation (M)
	// a complex mat mul is equivilent to 4 non-complex single matmuls
	// One assignment (M)
        mflops  = (num_iterations*(2.0*sizeof(TTYPE)/sizeof(float) *(M * M * M +M+M) ) )/ 1.0e6;
        secs    = cpu_time / ((double) CLOCKS_PER_SEC);
        mflop_s = mflops/secs;

        num_iterations *= 2;
    }
    printf ("Size: %u\ttime: %g us\n", M, (2*secs/num_iterations)*1e6);

    return mflop_s;
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
  			C=A*B;
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



template<class Mat_T>
double
time_mul2 ( Mat_T &A,  Mat_T &B, Mat_T& C)
{
    clock_t cpu_time;
    clock_t last_clock, assClock;
    double mflops, mflop_s;
    double secs = -1;

    int num_iterations = MIN_RUNS;
    int i, M=A.rows();

    while (secs < MIN_CPU_SECS) {

        cpu_time = 0;
		last_clock = clock();
        for (i = 0; i < num_iterations; ++i) {
  			C*=A;
  			C=B;
        }
        cpu_time += clock() - last_clock;

	//need to subtract out reset C=B time
	    assClock= clock();
	    for (i = 0; i < num_iterations; ++i) {
  			C=B;
        }
        cpu_time-=(clock()-assClock);


	//1 matrix multiplies (M^3)...
	// a complex value another 2x for 2 floating pt numbers
	// if TTYPE is a double the size is 2x that of a float
	// a complex mat mul is equivilent to 4 non-complex single matmuls
	//  thus another factor of 2 (we got a factor of 2 above already)
        mflops  = (sizeof(TTYPE)/sizeof(float)* num_iterations *( M * M * M))/ 1.0e6;
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
	int wh=1;
	cout<<"\nWhich test, prop[0], C=A*B[1], C*=A[2],Validate C=A*B [3], ,Validate C*=A [4], Validate C=prop(A,B) [5] or all [6]: ";
	cin>>wh;
	//  matrix_init (A);

	//  matrix_init (B);
	_matrix<TTYPE, FullMatrix> A,B,C, tmp;
	cout<<"Performing speed tests for Complex<float> with a size of "<<sizeof(TTYPE)<<endl;
	Random<UniformRandom<TTYPE > > mrY(TTYPE(1.00), TTYPE(2.0));
	Vector<double>  bllib_prop(N_SIZES);
	Vector<double>  bllib_mul1(N_SIZES);
	Vector<double>  bllib_mul2(N_SIZES);

if(wh==0 || wh>5){
	cout<<"Performing Propogator Tests C=prop(A,B)"<<endl;
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);
		mflop_s = time_prop(A, B, C);
		bllib_prop[sz_i]=mflop_s;

		printf ("Size: %u\tmflops: %g\n", M, mflop_s);
		//validate_dgemm (A, B, C);
	}
}
if(wh==1 || wh>5){
    cout<<"\nPerforming Multipication Tests C=A*B"<<endl;
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);
		mflop_s = time_mul1(A, B, C);
		bllib_mul1[sz_i]=mflop_s;

		printf ("Size: %u\tmflop/s: %g\n", M, mflop_s);
		//validate_dgemm (A, B, C);
	}
}
if(wh==3 || wh>5){
    cout<<"\nPerforming Multipication VALIDATION Tests C=A*B"<<endl;
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);
		C = A*B;
		validate_dgemm(A, B, C);
		printf ("Size: %u\t C= A*B Validated\n", M);
		//validate_dgemm (A, B, C);
	}
}

if(wh==4 || wh>5){
    cout<<"\nPerforming Multipication VALIDATION Tests C*=A"<<endl;
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);
		tmp = C;
		C *= A;
		validate_dgemm(A, tmp, C);
		printf ("Size: %u\t C*=A Validated\n", M);
		//validate_dgemm (A, B, C);
	}
}
if(wh==5 || wh>5){
    cout<<"\nPerforming Multipication VALIDATION Tests C=prop(A,B)"<<endl;
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);
		tmp = A*adjoint(B);
		C = prop(B,A);
		validate_dgemm(B, tmp, C);
		printf ("Size: %u\t C=prop(A,B) Validated\n", M);
		//validate_dgemm (A, B, C);
	}
}


if(wh==2 || wh>5){
	cout<<"\nPerforming Self Multipication Tests C*=A (C=A*C)"<<endl;
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);
		mflop_s = time_mul2(A, B, C);
		bllib_mul2[sz_i]=mflop_s;

		printf ("Size: %u\tmflop/s: %g\n", M, mflop_s);
		//validate_dgemm (A, B, C);
	}
}
	ofstream oo("bl_matspeed.m");
	oo<<"Infinity=Inf;"<<endl;
	oo<<"R=["<<test_sizes<<" "<<"];"<<endl;

	oo<<"bllib_prop=["<<bllib_prop<<" "<<"];"<<endl;
	oo<<"bllib_mul1=["<<bllib_mul1<<" "<<"];"<<endl;
	oo<<"bllib_mul2=["<<bllib_mul2<<" "<<"];"<<endl;

	oo<<"semilogx(R, bllib_prop,'r');"<<endl;
	oo<<"title('propagator(c=a*b*adjoint(a)) 700 MHz Pentium III Xeon');"<<endl;
	oo<<"xlabel('log(N)');"<<endl;
	oo<<"ylabel('MFLOPS');"<<endl;
	oo<<"legend('BlochLib 1.0--ATLAS-prop','BlochLib 1.0--ATLAS-C=A*B','BlochLib 1.0--ATLAS-C*=A');"<<endl;


     return 0;
}

