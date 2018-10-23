
#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

// Converts binary SpinSight NMR data (BIG ENDIAN)
// to a matlab file for a 1D spectrum


int main(int argc, char **argv)
{
    int q=1;
	std::string fname="", matname;
	std::cout<<std::endl<<"**BlochLib SpinSight (Chemmagnetics) 'directory' --> Matlab converter"<<std::endl;
	query_parameter(argc, argv, q++, "  --Enter SpinSight Data Directory: ", fname);
	query_parameter(argc, argv, q++, "  --Enter Output File Name: ", matname);
	int txt=0;
	query_parameter(argc, argv, q++, "  --Text[0] or Matlab[1] ", txt);


	SpinSightStream loadF(fname);
	double sw,dw;
	int TD, TD2;
	std::string nnuu;


	loadF.acq1D.get("al", TD);
	cout<<"  *points in 1D: "<<TD<<endl;

	loadF.acq1D.get("dw", dw);
	sw=1.0/dw;
	cout<<"  *Sweep width in 1D: "<<sw<<"Hz"<<endl;

	if(loadF.is2D()){
		loadF.acq2D.get("al2", TD2);
		cout<<"  *points in 2D: "<<TD2<<endl;
	}

	std::string vname=matname;
	matstream mm;

	if(txt==1){
		if(matname.find(".mat")>matname.size()){
			matname+=".mat";
		}else{
			vname=matname.substr(0, matname.find(".mat"));
		}
		mm.open(matname, ios::out | ios::binary);
	}

	matrix fid2d;
	Vector<complex> fid;
	if(!loadF.is2D()){
		loadF.read(fid);
		if(txt==1){
			mm.put(vname, fid);
			mm.put("sw1", sw);
			std::cout<<"  *Data Saved in '"<<matname
	                 <<"' under the variable '"<<vname<<"'"<<std::endl;
			std::cout<<"  *The Sweep Width information is also stored as variables"
	                 <<"\n   'sw1' for 1D info"<<std::endl;
		}else{
			std::ofstream of(matname.c_str());
			for(int i=0;i<fid.size();++i){
				of<<fid[i].Re()<<" "<<fid[i].Im()<<std::endl;
			}
			std::cout<<"  *Data Saved in '"<<matname
	                 <<"' (in cols <real> <imag>) "<<std::endl;
		}

	}else{
		loadF.read(fid2d);

		if(txt==1){
			mm.put(vname, fid2d);
			mm.put("sw1", sw);
			std::cout<<"  *Data Saved in '"<<matname
					 <<"' under the variable '"<<vname<<"'"<<std::endl;
			std::cout<<"  *The Sweep Width information is also stored as variables"
					 <<"\n   'sw1' for 1D info"
					 <<"\n   'sw2' is NOT known"<<std::endl;
		}else{
			std::ofstream of(matname.c_str());
			for(int i=0;i<fid.size();++i){
				for(int j=0;j<fid.size();++j){
					of<<i<<" "<<j<<" "<<fid2d(i,j).Re()<<" "<<fid2d(i,j).Im()<<std::endl;
				}
			}
			std::cout<<"  *Data Saved in '"<<matname
	                 <<"' (in cols <row> <col> <real> <imag>) "<<std::endl;
		}
	}
	return(0);
}
