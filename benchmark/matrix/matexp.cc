/*
  This is the matrix exponent driver

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

#define N_SIZES ((int) sizeof (test_sizes) / sizeof (int))

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
  			//C=A*B;
  			C=Mexp(A,1.0);
  			did++;
        }
        cpu_time += clock() - last_clock;

        mflops  = 8.0 * num_iterations * M * M * M / 1.0e6;
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
     hmatrix A,B,C;
   Random<UniformRandom<complex> > mrY(complex(0,2.0), complex(2.0,2.0));
	//Random<UniformRandom<> > mrY(0.0, 2.0);

     for (sz_i = 0; sz_i < N_SIZES; ++sz_i) {

           int M = test_sizes[sz_i];
		A.resize(M,M); A.apply(mrY);
		B.resize(M,M);	B.apply(mrY);
		C.resize(M,M);	C.apply(mrY);

		int did=0;
		double tt1=0.0, now;
		double one=1.0, zero=0.0;
		while(tt1<1.0){
			now=stopwatch();
			C=Mexp(A, 1.0);
			tt1+=stopwatch()-now;
			++did;
		}
		cout<<" General Matrix Multiply:: C=A*B; time: "<<tt1/did *1e6<<" us"<<endl;


         //validate_dgemm (M, A, B, C);
         // mflop_s = time_dgemm(M, A, B, C);
		 mflop_s = time_dgemm(A, B, C);
          printf ("Size: %u\tmflop/s: %g\n", M, mflop_s);
         //validate_dgemm (A, B, C);
     }

     return 0;
}

