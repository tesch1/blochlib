
#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/* A Test of the matlab io bits

requires a binary matlab file to work..

*/

int main(int argc,char* argv[]){ //get the file name of the matlab file
	int q=1;
	std::string fname;
	query_parameter(argc,argv,q++, "Enter File Name to read: ", fname);

	//open a matlab file for reading
	matstream mat(fname, ios::in | ios::binary);
	//print to the console what is in the file...
	//(this prints things exactly like 'whos' in matlab itself)
	mat.whos(cout);

	//get a complex vector from the inputfile named 'vdat'
	Vector<complex> loo;
	mat.get("vdat", loo);

	//get the simple variable 'zfill' from the file
	double zf;
	mat.get("zfill", zf);

	//close the input file...
	mat.close();

	//open a file for output...
	matstream outt("out.mat", ios::out | ios::binary);

	//define several vars to write...
	complex poo(7,8), hhh(3,4);
	Vector<complex> outmat(loo.size(), complex(6,6));

	//put these containers giveing them matlab names...
	outt.put("vdat", outmat);
	outt.put("poo",poo);
	outt.put("hhh",hhh);

	Vector<complex> hj(7, complex(3,2));
	Vector<int> kj(5,4);

	//puts matrices, vectors, numbers, etc..
	outt.put("loo", loo);
	outt.put("hj",hj);
	outt.put("kj",kj);
	outt.put("jjj",kj);
	matrix kji(5,5,3);
	outt.put("kji", kji);

	//put a more complex data strcuture into the file...
	//the output in the matlab file will be 'VecC0', 'VecC1', VecC2'
	//corresponding to each direction in the coord (x,y,z)
	Vector<coord<> > jk(35, 21);
	outt.put("VecC", jk);

	//close the output file
	outt.close();

	//open the ouput file and check output its contents to the console
	matstream mat2("out.mat", ios::in | ios::binary);
	mat2.whos(cout);
}
