
#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/* A Test of the matlab io bits

requires a binary matlab file to work..

*/

int main(int argc,char* argv[]){

	int q=1;
	std::string fname;
	query_parameter(argc,argv,q++, "Enter File Name to read: ", fname);

	Matlab5 mat(fname, ios::in | ios::binary);
	mat.whos(cout);
}
