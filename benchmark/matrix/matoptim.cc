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


timer stopwatch;
/*
  We try to run enough iterations to get reasonable timings.  The matrices
  are multiplied at least MIN_RUNS times.  If that doesn't take MIN_CPU_SECS
  seconds, then we double the number of iterations and try again.

  You may want to modify these to speed debugging...
*/
#define MIN_RUNS     4
#define MIN_CPU_SECS 0.25

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

#define TTYPE double

#define N_SIZES ((int) sizeof (test_sizes) / sizeof (int))

TTYPE A[MAX_SIZE * MAX_SIZE];
TTYPE B[MAX_SIZE * MAX_SIZE];
TTYPE C[MAX_SIZE * MAX_SIZE];

void
matrix_init (TTYPE *A)
{
     int i;

     for (i = 0; i < MAX_SIZE*MAX_SIZE; ++i) {
          A[i] = drand48 ();
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
void
validate_dgemm ( int M,
                 TT *A,  TT *B, TT *C)
{
    int i, j, k;

    matrix_clear (C);
 	TT one=1.0, zro=0.0;
  	DGEMM("N","N", &M,&M,&M,&one,A,&M,B,&M,&zro, C, &M);

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
			errorbound *= (A.rows() * DBL_EPSILON);

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

template<class TT>
double
time_dgemm ( int M,
             TT *A,  TT *B, TT *C)
{
    /*
      clock() normally measures milliseconds.  In non-Alpha Linux,
      though, it's limited by the HZ in include/asm-i386/param.h.
      The setting is, imho, unreasonably low for faster machines,
      but...  It's a subject of flame wars.  sigh.

      Timing under Linux is not fun.
    */

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
            //_bl_square_dgemm (M, A, B, C);
	//  DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
			DGEMM("N","N", &M,&M,&M,&one,A,&M,B,&M,&zro, C, &M);
        }
        cpu_time += clock() - last_clock;

        mflops  = 2.0 * num_iterations * M * M * M / 1.0e6;
        secs    = cpu_time / ((double) CLOCKS_PER_SEC);
        mflop_s = mflops/secs;

        num_iterations *= 2;
    }
    printf ("Size: %u\ttime: %g us\n", M, (2*secs/num_iterations)*1e6);

    return mflop_s;
}

template<class Mat_T>
double
time_dgemm ( Mat_T &A,  Mat_T &B, Mat_T& C)
{
    /*
      clock() normally measures milliseconds.  In non-Alpha Linux,
      though, it's limited by the HZ in include/asm-i386/param.h.
      The setting is, imho, unreasonably low for faster machines,
      but...  It's a subject of flame wars.  sigh.

      Timing under Linux is not fun.
    */

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
  			C=A*B*adjoint(A);
  			did++;
        }
        cpu_time += clock() - last_clock;

        mflops  = (2.0*8.0 * num_iterations * M * M * M +M)/ 1.0e6;
        secs    = cpu_time / ((double) CLOCKS_PER_SEC);
        mflop_s = mflops/secs;

        num_iterations *= 2;
    }
    printf ("Size: %u\ttime: %g us\n", M, (2*secs/did)*1e6);

    return mflop_s;
}

#define VALIDATE 1

int
main (void)
{
	int sz_i;
	double mflop_s;

	//  matrix_init (A);

	//  matrix_init (B);
	matrix A,B,C;
	Random<UniformRandom<complex> > mrY(complex(0,2.0), complex(2.0,2.0));
	//Random<UniformRandom<> > mrY(0.0, 2.0);

	double calctime, ltime, tottime;
	int size=((Nmax-Mmax)/5);
	int *logiter; logiter=new int[size];
	int ct=Mmax;
	for(int i=0;i<size;i++){
		 logiter[i]=ct;
		 ct+=5;
	}

	double  *gamma_dsymmF, *c_dsymmF;
	gamma_dsymmF=new double[size];
	c_dsymmF=new double[size];
	query_parameter(argc, argv, q++, "Enter number min Rows..", Mmax);
	query_parameter(argc, argv, q++, "Enter number max Rows..", Nmax);

    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);

/*		int did=0;
		double tt1=0.0, now;
		double one=1.0, zero=0.0;
		while(tt1<1.0){
			now=stopwatch();
			C=A*B;
			tt1+=stopwatch()-now;
			++did;
		}
		cout<<" General Matrix Multiply:: C=A*B; time: "<<tt1/did *1e6<<" us"<<endl;
*/

		mflop_s = time_dgemm(A, B, C);
		printf ("Size: %u\tmflop/s: %g\n", M, mflop_s);
		//validate_dgemm (A, B, C);
	}

     return 0;
}

