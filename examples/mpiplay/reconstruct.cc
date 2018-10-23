

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/**
This file is to test the workings of the MPI functions
with the various data types

**/



int main(int argc,char* argv[])
{

	MPIstart(argc, argv);
	std::cout<<MPIworld.name()<<"::"
			 <<MPIworld.rank()<<"/"<<MPIworld.size()<<std::endl;

//Vector<coord<> > test
	Vector<double > CCvec(10);
	for(int i=0;i<10;++i) CCvec[i]=i;
	int b=0, e=10, div=1;
	Range toadd=MPIworld.splitLoop(b,e,div);
	cout<<"Original rank: "<<MPIworld.rank()<<" [ "<<CCvec<<"]"<<endl;
	cout<<" Split Parameters: Begin: "<<b<<" End: "<<e<<" data size: "<<CCvec.size()<<endl;
	CCvec(toadd)+=(MPIworld.rank()+1)*100;
	cout<<"After Addtion: rank: "<<MPIworld.rank()<<" [ "<<CCvec<<"]"<<endl;
	MPIworld.reconstruct(CCvec, b,e);
	cout<<"After Reconstruct: rank: "<<MPIworld.rank()<<" [ "<<CCvec<<"]"<<endl;

	MPIworld.end();
}
