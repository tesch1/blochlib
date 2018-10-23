

#include <string>
#include "blochlib.h"
#include "spindex.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/*
Generates all the permutation sequences for Spin matrices (either Spherical order Cartisian)
for multiple spins and saves them to a binary file so other programs can read them in

this can the permutation generate ~4 million spins if both T1s ans T2s permutations
are desired...and will take a LONG time to generate them with no duplicates

be very careful

*/


int main(int argc,char* argv[])
{

	int q=1, intspin=0;
	std::string fout, foutlog, params, namef;
	query_parameter(argc, argv, q++,"File to parse... ", params);
	Parameters pset(params);
	pset.addSection("tensors");
	int order=pset.getParamI("order", "tensors");
	std::string corss=pset.getParamS("axis", "tensors");
	int cors=0;

//the Multipication type
	std::string multps=pset.getParamS("operator", "tensors");
	int multp;
	if(multps=="commutate"){ multp=TensorGenOps::Commutators;	}
	else if(multps=="add"){	multp=TensorGenOps::Add;	}
	else{	multp=TensorGenOps::Multiply;	}

	std::string sssp=pset.getParamS("spinops", "tensors");
	int ssp;

	Vector<int>  spins, spins2;
	if(sssp=="single"){
		spins.resize(1);
		spins[0]=pset.getParamI("singlespin", "tensors");
	}else if(sssp=="choose"){
		spins=pset.getParamVectorI("choosespins1", "tensors");
		if(corss=="T2" || corss=="T1T2"){
			spins2=pset.getParamVectorI("choosespins2", "tensors");
			if(spins.size() != spins2.size()){
				std::cerr<<std::endl<<"***Error: Both sets of spins MUST be the same length" <<std::endl;
				exit(-1);
			}
		}
	}else{
		spins.resize(1,-1);
	}

	fout=pset.getParamS("binout", "tensors");
	foutlog=pset.getParamS("logout", "tensors");
	namef=pset.getParamS("nameout", "tensors");

	ofstream logf(foutlog.c_str());

	logf<<" Tensor Generation log file..."<<std::endl;
	std::string h= " Created on: ";
	time_t now;
	time(&now);
	h+=ctime(&now);
	logf<<h<<std::endl;


	logf<<"Options:: "<<endl;
	if(corss=="xyz"){
		cors=TensorGenOps::I;
		logf<<"--Using Cartesian Tensors (Ix, Iy, Iz)"<<std::endl;
		logf<<std::endl<<"*There will be more then "<<sum_factorial(order)<<" unique permutations"<<std::endl;
	}else if(corss=="T1"){
		cors=TensorGenOps::T1;
		logf<<"--Using Spherical T1 Tensors (T1m1, T10, T11)"<<std::endl;
		logf<<std::endl<<"*There will be more then "<<sum_factorial(order)<<" unique permutations"<<std::endl;
	}else if(corss=="T2"){
		cors=TensorGenOps::T2;
		logf<<"--Using Spherical T2 Tensors (T2m2,T2m1, T20, T21, T22)"<<std::endl;
		logf<<std::endl<<"*There will be more then "<<sum_factorial(order)<<" unique permutations"<<std::endl;
	}else if(corss=="T1T2"){
		cors=TensorGenOps::T;
		logf<<"--Using Spherical T2 & T1 Tensors (T1m1, T10, T11, T2m2,T2m1, T20, T21, T22)"<<std::endl;
		logf<<std::endl<<"*There will be more then "<<sum_factorial(order)<<" unique permutations"<<std::endl;
	}else if(corss=="hamiltonians"){
		cors=TensorGenOps::Hamiltonians;
		logf<<"--Using Solid Systems Spin interations"<<std::endl;
		logf<<std::endl<<"*There will be more then "<<sum_factorial(order)<<" unique permutations"<<std::endl;
	}

	TensorGen tens;
	tens.setOrder(order);
	if(sssp=="single"){
		ssp=TensorGenOps::SingleSpin;
		logf<<"--Only Single Spin permutations"<<std::endl;
		logf<<"----Using spin: "<<intspin<<endl;
		tens.setSpin(intspin);
	}else if(sssp=="choose"){
		ssp=TensorGenOps::ChooseSpins;
		logf<<"--Choosing Spins..."<<std::endl;
		if(cors!=TensorGenOps::T2 && cors!=TensorGenOps::T ){
			logf<<"----Using spins: "<<spins<<endl;
			tens.setSpins(spins);
		}else{
			logf<<"----Using spins set 1: "<<spins<<endl;
			logf<<"----Using spins set 2: "<<spins2<<endl;
			tens.setSpins(spins, spins2);
		}
	}else{
		ssp=TensorGenOps::TotalSpace;
		logf<<"--Using Total Spin Space..."<<std::endl;
		tens.setSpins(-1);
	}


	tens.setOptions(cors | ssp |multp);
	logf<<std::endl<<"*Calculating all the permutations...removeing all duplicates"<<std::endl;
	SolidSys sys;
	if(corss=="hamiltonians"){
		pset.addSection("spins");
		sys=SolidSys(pset.section("spins"));
	}
	tens.generate(sys);

	fstream oo(fout.c_str(), ios::out | ios::binary);
	logf<<std::endl<<"*writing to binary file...\""<<fout<<"\""<<std::endl;
	oo<<tens;
	oo.close();
	fstream ii(fout.c_str(), ios::in | ios::binary);
	TensorGen inner;
	logf<<std::endl<<"*varifiying binary file by reading then dumping name to the name file.."<<std::endl;
	ii>>inner;

	fstream nameout(namef.c_str(),ios::out);

	for(int i=0;i<inner.size();++i){
		logf<<inner.name(sys,i)<<std::endl;
	}
}

