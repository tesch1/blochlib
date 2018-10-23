/*****************************************************************************

 */

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

//a function class to integrate (here it is sin(x))
//NOTE:: you need to define the "double operator(double)" function
class F
{
	public:
		F(){}
		double operator()(double x)
		{	return sin(x);	}
};



int main()
{


	F myF; //decalare a function


//declare the Integrate object first and then reuse it as nesseasry
	Integrate<F> myint(myF);
	ofstream out("out");
	double dt=.1, a=0, b=20, integ=0;

	//here we integrate F from t=0..20 in step sizes of dt=0.1
	while(a<b)
	{

//to Integrate from t1 to t2 simply use the "operator(t1, t2)" function
//or use the function "integrate(t1, t2)"

		//integ+=myint(a, a+dt); //the integration is ADDITIVE
		//or
		integ+=myint.integrate(a, a+dt);

		//print "<time> <orig> <integ>" to file
		out<<a<<" "<<myF(a+dt/2.0)<<" "<<integ<<endl;
		a+=dt;
	}


    return 1;
}

