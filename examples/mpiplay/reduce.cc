
#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

//An Examples on how to use 'reduce' for a simple sum

int main(int argv, char **argc){
	MPIworld.start(argv, argc); //start up MPI subsystem

	//print out some info....
	std::cout<<MPIworld.name()
	       <<" rank: "<<MPIworld.rank()
	       <<" num procs: "<<MPIworld.size()<<std::endl;

	Vector<double> vect(10,3.0);
	int begin=0, end=vect.size(), div=0;
	//for MPIrank==0, begin=0, end=5, div=5, splitR=Range(0,5)
	//for MPIrank==1, begin=5, end=10, div=5, splitR=Range(5,10)
	Range splitR=MPIworld.splitLoop(begin, end, div);
	cout<<"rank--"<<MPIworld.rank()<<" begin: "<<begin<<" end: "<<end<<" range: "<<splitR<<endl;
	//perform the sum
	double summ=sum(vect(splitR));
	std::cout<<"partial sum- "<<summ<<" -rank: "<<MPIworld.rank()<<" "<<std::endl;
	//reconstruct the master sum on the MPIrank=0
	MPIworld.reduce(summ, Reduce::Add);
	std::cout<<"partial sum- "<<summ<<" -rank: "<<MPIworld.rank()<<" "<<std::endl;

	//perform the sum
	double prodd=prod(vect(splitR));
	std::cout<<"partial Product- "<<prodd<<" -rank: "<<MPIworld.rank()<<" "<<std::endl;
	//reconstruct the master sum on the MPIrank=0
	MPIworld.reduce(prodd, Reduce::Multiply);
	std::cout<<"Total Product- "<<prodd<<" -rank: "<<MPIworld.rank()<<" "<<std::endl;

	if(prod(vect)!=prodd){
		if(MPIworld.master()) std::cout<<"MPIreduce, Multiply, FAILED!!!"<<std::endl;
	}else{
		if(MPIworld.master()) std::cout<<"MPIreduce, Multiply, PASSED!!!"<<std::endl;
	}

	if(sum(vect)!=summ){
		if(MPIworld.master()) std::cout<<"MPIreduce, Add, FAILED!!!"<<std::endl;
	}else{
		if(MPIworld.master()) std::cout<<"MPIreduce, Add, PASSED!!!"<<std::endl;
	}

	MPIworld.end();
	return 0;
}
