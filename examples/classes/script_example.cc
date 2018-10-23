

#include "blochlib.h"


// As simple on how-to extend the ScriptParse class

using namespace BlochLib;
using namespace std;

class myScriptClass:
	public ScriptParse
{

	private:

	//these functions are the actuall 'dooers' of
	//the script when the script hits these things

	//this performs a more complex print operation with
	// the input...uses the 'myPrint' keyword
		void domyPrint(std::string inSqe);

	//this dumps a matrix out to the screen
	// uses the key word 'makematrix'
		void doMakeMatrix(std::string inSqe);


	public:
	//wrap the basic constructos back to the
	//script parse class
		myScriptClass():
		   ScriptParse()
		{}
		myScriptClass(Parameters &pset):
		   ScriptParse(pset)
		{}
		myScriptClass(const Vector<std::string> &in):
		   ScriptParse()
		{}

	//this is the new decide function
		bool decide(std::string inSqe);
};

//create the functions

//first create our new decide function
bool myScriptClass::decide(std::string inSqe)
{
	if(inSqe.find("makematrix(")<inSqe.size()){
		doMakeMatrix(inSqe);
	}else if(inSqe.find("myPrint(")<inSqe.size()){
		domyPrint(inSqe);
	}else{
	//wrap any other input back to the
	// base class so that loop and assignments
	//and if's work as well
		return ScriptParse::decide(inSqe);
	}
	return true;
}


//this parse the 'myPrint' statment
// simply prints out what is given to it
//(DOES NOT evalute the expression)
//the syntax:: myprint(monkey)
// will print 'monkey' to cout
void myScriptClass::domyPrint(std::string inSqe)
{
	if(inSqe.find("myPrint(")<inSqe.size())
	{
		std::string tmS=getInside(inSqe); //get the stuff inside the '(' ')'
		if(tmS.size()==0 )
		{
			BLEXCEPTION(std::string("Bad 'on' usage for \"")+inSqe+"\" should be \"myprint(thing)\" ")
		}
	//just dump the input to the screen
		std::cout<<tmS<<endl;
	}
}


//this parse the 'makematrix' statment
// creates a real matrix, then prints it to the screen
//the syntax:: makematrix(3,3)
// creates a 3x3 real matrix filled with 0's
//the syntax:: makematrix(3,3,5)
// creates a 3x3 real matrix filled with 5's
void myScriptClass::doMakeMatrix(std::string inSqe)
{
	Vector<std::string> ps;
//remove any white spaces in the input
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("makematrix(")<tmS.size())
	{
		tmS=getInside(tmS); //get the stuff inside the '(' ')'
		ps=parseComma(tmS); //split the commas in the function
		if(ps.size()<=1 )
		{
			BLEXCEPTION(std::string("Bad on usage for \"")+inSqe+"\" should be \"makematrix(rows, cols)\" or \"makematrix(rows, cols, filler)\"")
		}
	//our matrix object
		rmatrix mymat;
	//the input parameters could be expression
	//so we need to parse them up using the 'myParse' object
		myParse.parse(ps[0]); //the rows
		int rows=int(myParse());
		myParse.parse(ps[1]); //the cols
		int cols=int(myParse());

		double filler=0;
	//get the filler if any
		if(ps.size()>2){
			myParse.parse(ps[2]);
			filler=myParse();
		}

	//resize the matrix
		mymat.resize(rows, cols, filler);
	//just dump the matrix to the screen
		std::cout<<mymat;
	}
}

/*
#---- an example input file ----

A=2;
loop(i=1:3)
   myPrint(on the matrix tonight )
   print(A)
   if(i==2)
      makematrix(i,A+3, 5)
   else
      makematrix(i,A)
   end
end


#------end example input file-----
*/

int main(int argc, char **argv)
{

	std::string fname;
	query_parameter(argc,argv,1, "Input File Name: ", fname);

	Parameters pset(fname);

	myScriptClass myScr;
	myScr.parse(pset);
}




