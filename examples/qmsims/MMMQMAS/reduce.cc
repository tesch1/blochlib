
#include "blochlib.h"
#include "propreduce.h"
using namespace BlochLib;
using namespace std;

int main(int argc,char* argv[]){

	int base, factor;
	query_parameter(argc, argv, 1, "Enter Base: ", base);
	query_parameter(argc, argv, 2, "Enter factor: ", factor);
	std::string fname;
	query_parameter(argc, argv, 3, "Enter log file name: ", fname);
	ofstream oo(fname.c_str());
	int method=0;
	query_parameter(argc, argv, 4, "method Squence Reduce[0], Time Shift Reduce[1]: ", method);
	
	PropReduce myReduce(base, factor, &oo);
	if(method == 0)	myReduce.sequenceReduce();
	else myReduce.timeShiftReduce();
}

