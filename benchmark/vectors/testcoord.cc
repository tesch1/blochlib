
/* testcoord.h ********/


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
	testcoord.h--> performs speed tests for 'coord<>' class
	
	to compile you must set "-ftemplate-depth-NN " where 'NN' must be larger then the 
	Coord length...'N'
	
*/


#include "blochlib.h"


timer stopwatch;
void printTime(int nrounds, std::string mess=""){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds)*1.e6 << " microseconds for "<<mess<<"\n";
}

void Info(std::string mess)
{
	cout<<mess;
	cout.flush();
}

template<class T>
void dummyf(T &oo){ 	}

#define LOOPER(len,  doto, NOPs, mess, VECT) \
	stopwatch.reset();	\
	for(int i=0;i<loops;++i){	\
		doto \
	}	\
	calctime=stopwatch(); \
 	stopwatch.reset(); 	\
   for(int i=0;i<loops;++i){ dummyf(a ## len); }	\
  ltime=stopwatch();	\
	tottime=calctime-ltime;	\
	cout<<"** "<<mess<<" **"<<std::endl;	\
	cout<<"time taken: "<<tottime*1.e6/double(loops)<<" microseconds"<<std::endl;	\
	cout<<"MFLOPS: "<<NOPs*double(len)/(tottime/double(loops))/1.e6<<std::endl;	\
	VECT.push_back(NOPs*double(len)/(tottime/double(loops))/1.e6); \
	

#define NITERATE(N) \
	std::string lenstr ## N="length: "+itost(N);	\
	std::string presstr ## N= " precision: double";	\
	coord<double, N> a ## N, b ## N, c ## N;	\
	for(int i=0;i<N;++i){	a ## N(i)=Rand();b ## N(i)=Rand();c ## N(i)=Rand();	}	\
	LOOPER(N, a ## N =prod(b ## N);,1,  "a=prod(b)..."+lenstr ## N+presstr ## N, prodF)	\
	LOOPER(N, a ## N =sum(b ## N);,1, "a=sum(b)..."+lenstr ## N+presstr ## N, sumF)	\
	LOOPER(N, a ## N =b ## N * c ## N;,1, "a=b*c..."+lenstr ## N+presstr ## N, mulF)	\
	LOOPER(N, a ## N =b ## N + c ## N;,1, "a=b+c..."+lenstr ## N+presstr ## N, addF)	\
	LOOPER(N, a ## N +=b ## N;,1, "a+=b..."+lenstr ## N+presstr ## N, saddF)	\
	LOOPER(N, a ## N *= c ## N + b ## N;,2, "a*=(c+b)..."+lenstr ## N+presstr ## N, muladdF)	\
	LOOPER(N, a ## N =dot(c ## N,b ## N);,2, "a=dot(c,b)..."+lenstr ## N+presstr ## N, dotF)	\

		
#define DAXPY(N) \
	std::string lenstr ## N="length: "+itost(N);	\
	std::string presstr ## N= " precision: double";	\
	coord<double, N> a ## N, b ## N, c ## N;	\
	double con ## N=Rand();	\
	for(int i=0;i<N;++i){	a ## N(i)=Rand();b ## N(i)=Rand();c ## N(i)=Rand();	}	\
	LOOPER(N, a ## N =a ## N + con ## N * b ## N ; , 8,  "a=a+const*b..."+lenstr ## N+presstr ## N, daxpyF)	\
	LOOPER(N, a ## N =a ## N + con ## N * b ## N ; , 8,  "a+=const*b..."+lenstr ## N+presstr ## N, daxpy2F)	\
			

int main(int argc, char **argv){
	
	int loops=0;
	int q=1;
	query_parameter(argc, argv, q++, "Enter number of loops for each test: ", loops);
	Vector<double> prodF, sumF, mulF, addF, saddF, muladdF, dotF, daxpyF, daxpy2F;
	
	double calctime, ltime, tottime;
//	NITERATE(1)
//	NITERATE(2)
//	NITERATE(5)
//	NITERATE(10)
//	NITERATE(20)
//	NITERATE(50)
//	NITERATE(100)
//	NITERATE(200)

	DAXPY(1)
	DAXPY(2)
	DAXPY(3)
	DAXPY(4)
	DAXPY(5)
	DAXPY(6)
	DAXPY(7)
	DAXPY(8)
	DAXPY(9)
	DAXPY(10)
	DAXPY(20)
	DAXPY(50)
	DAXPY(100)				

	ofstream oo("coordspeed");
	oo<<"R=[1 2 3 4 5 6 7 8 9 10 20 50 100]"<<endl;
	oo<<"prodF=["<<prodF<<"]"<<endl;
	oo<<"sumF=["<<sumF<<"]"<<endl;
	oo<<"mulF=["<<mulF<<"]"<<endl;
	oo<<"addF=["<<addF<<"]"<<endl;
	oo<<"saddF=["<<saddF<<"]"<<endl;
	oo<<"muladdF=["<<muladdF<<"]"<<endl;
	oo<<"dotF=["<<dotF<<"]"<<endl;
	oo<<"daxpyF=["<<daxpyF<<"]"<<endl;
	oo<<"daxpy2F=["<<daxpy2F<<"]"<<endl;




}
