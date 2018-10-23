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

#include "gamma.h"


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

#define TTYPE double

#define N_SIZES ((int) sizeof (test_sizes) / sizeof (int))


void matrix_init (matrix A, int on)
{
     int i,j;
     for (i = 0; i < on; ++i) {
     	for (j = 0; j < on; ++j) {
          Re(A(i,j),drand48 ());
          Im(A(i,j),drand48 ());
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


template<class Mat_T>
double
time_dgemm ( Mat_T &A,  Mat_T &B, Mat_T& C)
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


int
main (void)
{
	int sz_i;
	double mflop_s;

	double  *gamma_prop; gamma_prop=new double[N_SIZES];

    for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {
		int M = test_sizes[sz_i];
		matrix A(M,M),B(M,M),C(M,M);
		matrix_init(A,M);
		matrix_init(B,M);
		matrix_init(C,M);
		mflop_s = time_dgemm(A, B, C);
		gamma_prop[sz_i]=mflop_s;

		printf ("Size: %u\tmflop/s: %g\n", M, mflop_s);
		//validate_dgemm (A, B, C);
	}

	ofstream oo("gam_matspeed.m");
	oo<<"Infinity=Inf;"<<endl;
	oo<<"R=[";
	for(int sz_i=0;sz_i< N_SIZES; ++sz_i)	oo<<test_sizes[sz_i]<<" ";
	oo<<"];"<<endl;

	oo<<"gamma_prop=[";
	for(int k=0;k<N_SIZES;++k)	oo<<gamma_prop[k]<<" ";
	oo<<"];"<<endl;

	oo<<"semilogx(R, gamma_prop,'r');"<<endl;
	oo<<"title('propagator(c=a*b*adjoint(a)) 700 MHz Pentium III Xeon');"<<endl;
	oo<<"xlabel('log(N)');"<<endl;
	oo<<"ylabel('MFLOPS');"<<endl;
	oo<<"legend('Gamma 4.0.5--matrix');"<<endl;
	delete gamma_prop;

     return 0;
}

