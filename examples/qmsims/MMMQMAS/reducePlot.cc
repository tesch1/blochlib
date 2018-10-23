
#include "blochlib.h"
#include "propreduce.h"
using namespace BlochLib;
using namespace std;

//dumps out the presentage of matrix mul reductions
// using the rational reduction method give the obeervation factor
// up to a periodicity factor of 500

int main(int argc,char* argv[]){

	int factor;
	query_parameter(argc, argv, 1, "Enter factor: ", factor);
	//std::string fname;
	//query_parameter(argc, argv, 2, "Enter log file name: ", fname);
	//ofstream oo(fname.c_str());

	int len=120;
	int maxBase=factor+1+len;
	int stp=int(double(maxBase)/double(len));
	Vector<int> bases(len); bases[0]=factor+1;
	Vector<int> origMult(len);origMult[0]=factor*(factor+1);
	Vector<int> bestMuls(len);

	for(int i=1;i<len;++i){
		bases[i]=bases[i-1]+stp;
		origMult[i]=factor*bases[i];
	}

	for(int i=0;i<len;++i){
		cout<<"Factor: "<<factor<<" On Base: "<<bases[i]; cout.flush();
		if(bases[i]%factor!=0){
			PropReduce myReduce(bases[i], factor);
			myReduce.reduce();
			bestMuls[i]=myReduce.bestMultiplications();
		}else{
			bestMuls[i]=bases[i];
		}
		cout<<" orig MMs: "<<origMult[i]<<" best MM: "<<bestMuls[i]<<endl;


	}

	ofstream mmout("plotBest.m");
	mmout<<"n="<<factor<<";"<<endl;
	mmout<<"m=["<<bases<<"];"<<endl;
	mmout<<"origM=["<<origMult<<"];"<<endl;
	mmout<<"redM=["<<bestMuls<<"];"<<endl;
	mmout<<"pres=-100*((redM./origM)-1);"<<endl;
	mmout<<"plot(m, pres);"<<endl;
	mmout<<"xlabel('perodicity factor, m'); ylabel('reduction (%)')"<<endl;
	mmout<<"title(['percent reduction for observation factor ' num2str(n)]);"<<endl;



}

