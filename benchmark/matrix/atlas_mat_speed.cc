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
#define MULTYPE DGEMM

#define N_SIZES ((int) sizeof (test_sizes) / sizeof (int))

TTYPE A[MAX_SIZE * MAX_SIZE];
TTYPE B[MAX_SIZE * MAX_SIZE];
TTYPE C[MAX_SIZE * MAX_SIZE];

void
matrix_init (TTYPE *A)
{
     int i;

     for (i = 0; i < MAX_SIZE*MAX_SIZE; ++i) {
          Re(A[i],drand48 ());
          Im(A[i],drand48 ());
     }
}

void
matrix_clear (TTYPE *C)
{
     memset (C, 0, MAX_SIZE * MAX_SIZE * sizeof (TTYPE));
}

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
 template<class TT>
void validate_dgemm ( int M,TT *A,  TT *B, TT *C)
{
    int i, j, k;

    matrix_clear (C);
 	TT one=1.0, zro=0.0;
  	MULTYPE("N","N", &M,&M,&M,&one,A,&M,B,&M,&zro, C, &M);

    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j) {

            TT dotprod = 0;
            TT errorbound = 0;
            double err;

            for (k = 0; k < M; ++k) {
                TT prod = A[k*M + i] * B[j*M + k];
                dotprod += prod;
                errorbound += abs(prod);
            }
            errorbound *= (M * DBL_EPSILON);

            err = abs(C[j*M + i] - dotprod);
            if (err > 3*abs(errorbound)) {
                printf("Matrix multiply failed.\n");
                printf("C(%d,%d) should be (%g, %g), was (%g, %g)\n", i, j,
                       Re(C[j*M + i]), Im(C[j*M + i]), Re(dotprod),Im(dotprod));
                printf("Error of %g, acceptable limit %g\n",
                       err, 3*abs(errorbound));
                exit(-1);
            }
        }
    }
}


template<class TT>
double time_prop ( int M, TT *A,  TT *B, TT *C)
{

    clock_t cpu_time;
    clock_t last_clock;
    double mflops, mflop_s;
    double secs = -1;

    int num_iterations = MIN_RUNS;
    int i;
	TT one=1.0, zro=0.0;
    while (secs < MIN_CPU_SECS) {

        cpu_time = 0;

        matrix_clear (C);
        last_clock = clock();
        for (i = 0; i < num_iterations; ++i) {
   			TTYPE *TMMMAT; TMMMAT=new TTYPE[M * M];
	    	MULTYPE("N","C", &M,&M,&M,&one,B,&M,A,&M,&zro, TMMMAT, &M);
   			MULTYPE("N","N", &M,&M,&M,&one,A,&M,TMMMAT,&M,&zro, C, &M);
   			delete [] TMMMAT;
        }
        cpu_time += clock() - last_clock;

       	//2-->doubles
       	//2-->2 matrix muls
       	//4-->one complex mat mul is 4 double mat muls
       	// MxMxM--> normal N^3 behavior
       	// M the adjoint operation
       	mflops  =( sizeof(TTYPE)/sizeof(float) *2.0* num_iterations * (M * M * M +M))/ 1.0e6;
        secs    = cpu_time / ((double) CLOCKS_PER_SEC);
        mflop_s = mflops/secs;

        num_iterations *= 2;
    }
    printf ("Size: %u\ttime: %g us\n", M, (2*secs/num_iterations)*1e6);

    return mflop_s;
}

template<class TT>
double time_mul1 ( int M, TT *A,  TT *B, TT *C)
{
    clock_t cpu_time;
    clock_t last_clock;
    double mflops, mflop_s;
    double secs = -1;

    int num_iterations = MIN_RUNS;
    int i;
	TT one=1.0, zro=0.0;
    while (secs < MIN_CPU_SECS) {

        cpu_time = 0;

        matrix_clear (C);
        last_clock = clock();
        for (i = 0; i < num_iterations; ++i) {
   			MULTYPE("N","N", &M,&M,&M,&one,A,&M,B,&M,&zro, C, &M);
        }
        cpu_time += clock() - last_clock;

       	//MxMxM-->1 matrix muls
       	mflops  =( sizeof(TTYPE)/sizeof(float)* num_iterations * (M * M * M ))/ 1.0e6;
        secs    = cpu_time / ((double) CLOCKS_PER_SEC);
        mflop_s = mflops/secs;

        num_iterations *= 2;
    }
    printf ("Size: %u\ttime: %g us\n", M, (secs/num_iterations)*1e6);

    return mflop_s;
}
int main (void)
{
	int sz_i;
	double mflop_s;

	matrix_init (A);
	matrix_init (B);
	Random<UniformRandom<TTYPE> > mrY(TTYPE(-1), TTYPE(1.0));
	Vector<double>  atlas_prop(N_SIZES), atlas_mul1(N_SIZES);

  	int wh=1;
	cout<<"Which test, prop[0], C=A*B[1], Validate [2] or all [3]: ";
	cin>>wh;

 if(wh==0 || wh>3){
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		mflop_s = time_prop(M,A, B, C);
		atlas_prop[sz_i]=mflop_s;
		printf ("Size: %u\tmflop/s: %g\n", M, mflop_s);
	}
}

if(wh==2 || wh>3){
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		validate_dgemm(M,A, B, C);
		printf ("Size: %u \t Validated\n", M, mflop_s);
	}
	exit(0);
}

if(wh==1 || wh>3){
    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		mflop_s = time_mul1(M,A, B, C);
		atlas_mul1[sz_i]=mflop_s;
		printf ("Size: %u\tmflop/s: %g\n", M, mflop_s);
	}
}


	ofstream oo("at_matspeed.m");
	oo<<"Infinity=Inf;"<<endl;
	oo<<"R=[";
	for(int sz_i=0;sz_i< N_SIZES; ++sz_i)	oo<<test_sizes[sz_i]<<" ";
	oo<<"];"<<endl;

	oo<<"at_prop=["<<atlas_prop<<"];"<<endl;
	oo<<"at_mul1=["<<atlas_mul1<<"];"<<endl;

	oo<<"semilogx(R, at_prop,'r',R, at_mul,'r*');"<<endl;
	oo<<"title('propagator(c=a*b*adjoint(a)) 700 MHz Pentium III Xeon');"<<endl;
	oo<<"xlabel('log(N)');"<<endl;
	oo<<"ylabel('MFLOPS');"<<endl;
	oo<<"legend('ATLAS-3.4.2--C=A.B.A''',ATLAS-3.4.2--C=A.B' );"<<endl;

	matrix_clear(A);
	matrix_clear(B);
	matrix_clear(C);
	//matrix_clear(TMMMAT);


     return 0;
}

