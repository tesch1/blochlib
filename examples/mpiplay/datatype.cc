

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/**
This file is to test the workings of the MPI functions
with the various data types

**/

#define TEST_MPI_B(T) \
	T rr ## T;	\
	if(MPIworld.master()){	\


#define TEST_MPI_E(T) \
		for(int i=1;i<MPIworld.size();++i) {	\
			std::cout<<"Sending "<<#T<<" To proc #"<<i<<std::endl<<" [ "<<rr ## T<<" ]"<<std::endl;	\
			MPIworld.put(rr ## T, i);	\
		}	\
	}else{	\
		std::cout<<"Getting "<<#T<<" On proc #"<<MPIworld.rank();	\
		MPIworld.get(rr ## T, 0);	\
		std::cout<<std::endl<<" [ "<<rr ## T<<" ]"<<std::endl;	\
	}	\
	MPIworld.advanceTag(); \


#define TEST_MPI_B2(T, T1) \
	T ## < ## T1 ## > rr ## T ## T1;	\
	if(MPIworld.master()){	\


#define TEST_MPI_E2(T, T1) \
		for(int i=1;i<MPIworld.size();++i) {	\
			std::cout<<"Sending "<<#T<<"<"<<#T1<<"> To proc #"<<i<<std::endl<<" [ "<<rr ## T ## T1<<" ]"<<std::endl;	\
			MPIworld.put(rr ## T ## T1, i);	\
		}	\
	}else{	\
		std::cout<<"Getting "<<#T<<"<"<<#T1<<"> On proc #"<<MPIworld.rank();	\
		MPIworld.get(rr ## T ## T1, 0);	\
		std::cout<<std::endl<<" [ "<<rr ## T ## T1<<" ]"<<std::endl;	\
	}	\
	MPIworld.advanceTag();	\



#define TEST_MPI_B3(T, T1, T2) \
	T ## < ## T1 ## , ## T2 ##> rr ## T ## T1 ## T2;	\
	if(MPIworld.master()){	\


#define TEST_MPI_E3(T, T1, T2) \
		for(int i=1;i<MPIworld.size();++i) {	\
			std::cout<<"Sending "<<#T<<"<"<<#T1<<","<<#T2<<"> To proc #"<<i<<std::endl<<" [ "<<rr ## T ## T1 ## T2<<" ]"<<std::endl;	\
			MPIworld.put(rr ## T ## T1 ## T2, i);	\
		}	\
	}else{	\
		std::cout<<"Getting "<<#T<<"<"<<#T1<<","<<#T2<<"> On proc #"<<MPIworld.rank();	\
		MPIworld.get(rr ## T ## T1, 0);	\
		std::cout<<std::endl<<" [ "<<rr ## T ## T1 ## T2<<" ]"<<std::endl;	\
	}	\
	MPIworld.advanceTag();	\


int main(int argc,char* argv[])
{
	if(argc==1){
		std::cout<<std::endl<<" This program is made to test the MPI capabilites"<<std::endl
		<<" Of this machine...so perhaps you should run this program using "<<std::endl
		<<" 'mpirun -np 2 "<<argv[0]<<"'"<<std::endl<<std::endl;
	}

#ifdef HAVE_MPI
//	MPI_Init(&argc, &argv);
	MPIworld.start(argc, argv);
	std::cout<<MPIworld.name()<<"::"<<MPIworld.rank()<<"/"<<MPIworld.size()<<std::endl;


//int test
	TEST_MPI_B(int)
		rrint=MPIworld.rank();
	TEST_MPI_E(int)

//float test
	TEST_MPI_B(float)
		rrfloat=MPIworld.rank()+0.6;
	TEST_MPI_E(float)

//double test
	TEST_MPI_B(double)
		rrdouble=MPIworld.rank()+0.8;
	TEST_MPI_E(double)

//complex test
	TEST_MPI_B(complex)
		rrcomplex=complex(MPIworld.rank(),2);
	TEST_MPI_E(complex)

//coord<int,3> test
	TEST_MPI_B2(coord,int)
		rrcoordint(0, 1,2);
	TEST_MPI_E2(coord,int)

//coord<double,3> test
	TEST_MPI_B2(coord,double)
		rrcoorddouble(0, 5,10);
	TEST_MPI_E2(coord,double)

//coord<complex,3> test
	TEST_MPI_B2(coord,complex)
		rrcoordcomplex=coord<>(0, 5,10)+complex(3,2);
	TEST_MPI_E2(coord,complex)

//Vector<int> test
	TEST_MPI_B2(Vector,int)
		rrVectorint.resize(5);
		for(int i=0;i<5;++i) rrVectorint[i]=i;
	TEST_MPI_E2(Vector,int)

//Vector<double> test
	TEST_MPI_B2(Vector,double)
		rrVectordouble.resize(5);
		for(int i=0;i<5;++i) rrVectordouble[i]=double(i)+0.5;
	TEST_MPI_E2(Vector,double)

//Vector<complex> test
	TEST_MPI_B2(Vector,complex)
		rrVectorcomplex.resize(5);
		for(int i=0;i<5;++i) rrVectorcomplex[i]=complex(i,i+1);
	TEST_MPI_E2(Vector,complex)

//Vector<coord<> > test
	Vector<coord<> > CCvec(5);
	if(MPIworld.master()){
		for(int i=0;i<5;++i) CCvec[i]=coord<>(i,i+1,i+8);
		for(int i=1;i<MPIworld.size();++i) {
			std::cout<<"Sending Vector<coord<> > To proc #"<<i<<std::endl<<" [ "<<CCvec<<"]"<<std::endl;
			MPIworld.put(CCvec, i);
		}
	}else{
		std::cout<<"Getting Vector<coord<> > On proc #"<<MPIworld.rank();
		MPIworld.get(CCvec, 0);
		std::cout<<std::endl<<" [ "<<CCvec<<" ]"<<std::endl;
	}
	MPIworld.advanceTag();

//rmatrix test
	TEST_MPI_B(rmatrix)
		rrrmatrix.resize(5,4);
		for(int i=0;i<5;++i)
			for(int j=0;j<4;++j)
				rrrmatrix(i,j)=i+j;
	TEST_MPI_E(rmatrix)

//matrix test
	TEST_MPI_B(matrix)
		rrmatrix.resize(5,2);
		for(int i=0;i<5;++i)
			for(int j=0;j<2;++j)
				rrmatrix(i,j)=complex(i+j, i);
	TEST_MPI_E(matrix)

//hmatrix test
	TEST_MPI_B(hmatrix)
		rrhmatrix.resize(4,4);
		for(int i=0;i<4;++i)
			for(int j=0;j<4;++j)
				rrhmatrix(i,j)=complex(i+j, j);
	TEST_MPI_E(hmatrix)

//smatrix test
	TEST_MPI_B(smatrix)
		rrsmatrix.resize(4,4);
		for(int i=0;i<4;++i)
			for(int j=0;j<4;++j)
				rrsmatrix(i,j)=i*j;
	TEST_MPI_E(smatrix)

//dmatrix test
	TEST_MPI_B(dmatrix)
		rrdmatrix.resize(4,4);
		for(int i=0;i<4;++i)
			rrdmatrix(i,i)=i*5;
	TEST_MPI_E(dmatrix)

//imatrix test
	TEST_MPI_B(imatrix)
		rrimatrix.resize(4,4);
	TEST_MPI_E(imatrix)

	MPIworld.end();

#endif


}
