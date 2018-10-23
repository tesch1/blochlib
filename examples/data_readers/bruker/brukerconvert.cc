
#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace BL_endians; //for the endian conversions
using namespace std;

// Converts binary Bruker Xwin NMR data (BIG ENDIAN)
// to a matlab file for a 1D spectrum


int main(int argc, char **argv)
{
    int q=1;
	std::string fname="", matname;
	int choppts;
	std::cout<<std::endl<<"**BlochLib XWINNMR (Bruker) 'directory' --> Matlab converter"<<std::endl;
	query_parameter(argc, argv, q++, "  --Enter Burker Data Directory: ", fname);
	query_parameter(argc, argv, q++, "  --Enter Output File Name: ", matname);
	int txt=0;
	query_parameter(argc, argv, q++, "  --Text[0] or Matlab[1] ", txt);

	std::string Brmessage=" --Bruker tends to add ~70 pts of bad \n";
	Brmessage+="   'Digitizer garbage' at the begining of \n";
	Brmessage+="   each FID, should i remove these points? \n";
	Brmessage+="  -if so how Many pts to remove? (enter '0' for none):";

	query_parameter(argc, argv, q++,Brmessage, choppts);

	XWINNMRstream loadF(fname);
	double sw, sw2;
	int TD, TD2;
	std::string nnuu;

	loadF.acq1D.get("SW_h", sw);
	cout<<"  *Sweep width in 1D: "<<sw<<"Hz"<<endl;

	loadF.acq1D.get("TD", TD);
	cout<<"  *points in 1D: "<<TD/2<<endl;

	if(loadF.is2D()){
		loadF.acq2D.get("SW_h", sw2);
		cout<<"  *Sweep width in 2D: "<<sw2<<"Hz"<<endl;

		loadF.acq2D.get("TD", TD2);
		cout<<"  *points in 2D: "<<TD2<<endl;
	}

	matstream mm;
	std::string vname=matname;

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
		if(choppts>0)
		{	fid=fid(Range(choppts, fid.size()-1));	}

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
		matrix tmFid=fid2d;
		if(choppts>0)
		{	tmFid=fid2d(Range(0, fid2d.rows()-1),Range(choppts, fid2d.cols()-1)) ;	}

		if(txt==1){
			mm.put(vname, tmFid);
			mm.put("sw1", sw);
			mm.put("sw2", sw2);
			std::cout<<"  *Data Saved in '"<<matname
					 <<"' under the variable '"<<vname<<"'"<<std::endl;
			std::cout<<"  *The Sweep Width information is also stored as variables"
					 <<"\n   'sw1' for 1D info"
					 <<"\n   'sw2' for 2D info"<<std::endl;
		}else{
			std::ofstream of(matname.c_str());
			for(int i=0;i<fid2d.rows();++i){
				for(int j=0;j<fid2d.cols();++j){
					of<<i+1<<" "<<j+1<<" "<<fid2d(i,j).Re()<<" "<<fid2d(i,j).Im()<<std::endl;
				}
			}
			std::cout<<"  *Data Saved in '"<<matname
	                 <<"' (in cols <row> <col> <real> <imag>) "<<std::endl;
		}
	}
	std::cout<<std::endl;

	return(0);
}
