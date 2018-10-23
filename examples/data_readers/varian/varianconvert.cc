
#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

// Converts binary Varian INOVA Binary data
//(COLLECTED ON A SUN ie the file has BigEndian ordering)
// to a matlab file for a 1D spectrum

int main(int argc, char **argv){
   	int q=1;
	std::string fname="";
	std::string matname;
	std::cout<<std::endl<<"**BlochLib VNMR (Varian) 'directory' --> Matlab converter"<<std::endl;
	query_parameter(argc, argv, q++, "  --Enter Inova Directory: ", fname);
	query_parameter(argc, argv, q++, "  --Enter Output file name: ", matname);
	int ch=0;
	query_parameter(argc, argv, q++, "  --Output Text[0] or Matlab[1]?: ", ch);


	matrix fids;

	VNMRdata vfile;
	vfile.open(fname,ios::binary | ios::in);
	vfile.read(fids);


	if(ch==1){
		std::string vname=matname;
		if(matname.find(".mat")>matname.size()){
			matname+=".mat";
		}else{
			vname=matname.substr(0, matname.find(".mat"));
		}
		matstream matout(matname.c_str(), ios::binary | ios::out);
		matout.put(vname, fids);
		matout.close();
		std::cout<<"  *Data saved in Matlab file "<<matname
		         <<" under variable '"<<vname<<"'"<<std::endl;
	}else{
		ofstream oo(matname.c_str());
		if(!vfile.is2D()){
			for(int i=0;i<fids.rows();++i){
				oo<<fids(i,0).Re()<<" "<<fids(i,0).Im()<<std::endl;
			}
		}else{
			for(int i=0;i<fids.rows();++i){
				for(int j=0;j<fids.cols();++j){
					oo<<i<<" "<<j<<" "<<fids(i,j).Re()<<" "<<fids(i,j).Im()<<std::endl;
				}
			}
		}
		std::cout<<" *Data saved in ASCII file "<<matname
		         <<"  (cols as <row> <col> <real> <imag>)"<<std::endl;
	}
}
