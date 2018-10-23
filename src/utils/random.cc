/* random.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-8-01
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
 	random.cc-->a collection of random number generators and classes
 */


 #ifndef _bloch_random_cc_
 #define _bloch_random_cc_


#include "utils/random.h"
#include "utils/constants.h"
#  include <time.h>
#include<math.h>

BEGIN_BL_NAMESPACE


//A small constant class to generate
// random numbers
//This generator came from
/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* Any feedback is very welcome. For any question, comments,       */
/* see http://www.math.keio.ac.jp/matumoto/emt.html or email       */
/* matumoto@math.keio.ac.jp */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */



//make the 'static values' linkable...
int BL_Random::mti=BL_Random::N+1;
unsigned long BL_Random::mt[N];
time_t BL_Random::rand_timer;

//seed it initially with a timer value
BL_Random::BL_Random()
{
//the 'plus' one here is because the abve generator cannot handle
// a '0' starting value
	mti=N+1;
	seed((unsigned long)(time(&rand_timer))+1);
}


void BL_Random::seed(unsigned long seedr)
{
	/* setting initial seeds to mt[N] using         */
	/* the generator Line 25 of Table 1 in          */
	/* [KNUTH 1981, The Art of Computer Programming */
	/*    Vol. 2 (2nd Ed.), pp102]                  */
	mt[0]= seedr & 0xffffffff;
	for (mti=1; mti<N; mti++)
		mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;

    /* Using seeds that are derived from the above
       generator to seed sgenrand leads to shifted
       series. Furthermore, the least-significant bit
       is always odd, or always even.
       The following code should prevent that. The
       random generator used is more or less arbitrary, but
       it has a reasonably long period (1825731182) and
       should generate well-mixed bit-streams. */
    unsigned long s = 373737;
    for (mti=1; mti<N; mti++)
    {
		mt[mti] ^= s;
		s = s * 5531 + 81547;
		s ^= (s >> 9) ^ (s << 19);
    }
}


//integrer random numbers
unsigned long BL_Random::randomI()
{
	unsigned long y;
	static unsigned long mag01[2]={0x0, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N) { /* generate N words at one time */
		int kk;

		if (mti == N+1)   /* if seed() has not been called, */
			seed((unsigned long)(time(&rand_timer))+1);

		for (kk=0;kk<N-M;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for (;kk<N-1;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

		mti = 0;
	}

	y = mt[mti++];
	y ^= y >> 11;
	y ^= y << 7 & TEMPERING_MASK_B;
	y ^= y << 15 & TEMPERING_MASK_C;
	y ^= y >> 18;
	return y;
}

//double random generator
double BL_Random::randomD()
{
	return ( (double)randomI() / (unsigned long)0xffffffff ); /* reals */
}

double Rand()
{
	return BL_Random::randomD();
}


//Gaussian Random number.. with unit variance and 0 mean
double GaussianRand()
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	if  (iset == 0) {
		do {
			v1=2.0*Rand()-1.0;
			v2=2.0*Rand()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0);
		fac=pow(-2.0* log(rsq)/rsq, 0.5);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

double GammaLog(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

//Poisson Distributed Random Number...distributed AROUND the mean 'xm'
double PoissonRand(double xm)
{
	static double sq,alxm,g,oldm=(-1.0);
	double em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= Rand();
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=pow(2.0*xm, 0.5);
			alxm=log(xm);
			g=xm*alxm-GammaLog(xm+1.0);
		}
		do {
			do {
				y=tan(PI*Rand());
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-GammaLog(em+1.0)-g);
		} while (Rand() > t);
	}
	return em;
}


END_BL_NAMESPACE



#endif


