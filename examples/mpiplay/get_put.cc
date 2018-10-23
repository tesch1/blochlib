
#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

//a simple example fo get, put and scatter using the 'MPIworld' controller

int main(int argc,char* argv[])
{
//start up MPI subsystem
	MPIworld.start(argc, argv);

//print out some info....
	std::cout<<MPIworld.name()
	         <<" rank: "<<MPIworld.rank()
	         <<" num procs: "<<MPIworld.size()<<std::endl;

//create some vector to move around
	Vector<double> tmD(Spread<double>(1.0,4.0,1));
	Range All(Range::Start, Range::End);
	if(MPIworld.master()) tmD(All)=tmD(All)*tmD(All);
	std::cout<<"rank: "<<MPIworld.rank()<<" [ "<<tmD<<"] "<<std::endl;

	if(MPIworld.master()){
//put the manipulated vector to the other procs
		for(int i=1;i<MPIworld.size();++i) MPIworld.put(tmD, i);
	}else{
//get the vector from proc 0
		MPIworld.get(tmD,0);
	}
	std::cout<<"After put/get--rank: "<<MPIworld.rank()<<" [ "<<tmD<<"]"<<std::endl;

	// Or use the 'MPIscatter' to distribute the new value
	if(MPIworld.master()) tmD*=9;
	MPIworld.scatter(tmD);
	std::cout<<"After scatter--rank: "<<MPIworld.rank()<<" [ "<<tmD<<"]"<<std::endl;



	//end this very boring mpi session
	MPIworld.end();
	return 0;
}
