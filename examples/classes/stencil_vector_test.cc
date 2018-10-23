

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

timer stopwatch;
void printTime(int nrounds=1){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}


int main()
{

	ofstream oo("moo"), oo2("loo");

	double div=0.1;
	Vector<double> kk(Spread<double>(-10,10,div)); kk=cos(kk);
	Vector<double> ans(Spread<double>(-10,10,div));
	ans=Derivative_2_4n( kk)/(div*div)-Derivative_3_4n( kk)/(div*div*div);
	//ans/=0.01;
	Vector<double>::iterator myit(kk);
	myit.begin(5);
	myit.end(kk.size()-5);
	while(myit)
	{
		oo<<myit()<<endl;
		oo2<<ans(myit.curpos())<<endl;
		++myit;
	}
	printTime();
}
