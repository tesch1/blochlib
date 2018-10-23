


#include "blochlib.h"

using namespace BlochLib;
using namespace std;

template<class NMRReader>
bool ReadData(NMRReader &reader, std::string fname, matrix &fids, Vector<complex> &fid, bool &is2D){
	if(!reader.open(fname, false)){
		return false;
	}
	if(reader.is2D()){
		is2D=true;
		return reader.read(fids);
	}
	reader.read(fid);
	is2D=false;
	return true;
}


int main(int argc, char* argv[])
{
	if (argc <3){
		cout << "reads and writes a .WAVE (audio) file from NMR data" << endl;
		cout << " can read 'SpinSight', 'XWINNMR', or 'VNMR' data" << endl;
		cout << " if the data is 2D, it will create multiple files" << endl;
		cout << " named {out}i.wav, where {out} is the base, and 'i'" << endl;
		cout << " is the fid index." << endl;
		cout << " ----------------------------" << endl;
		cout << " INPUT: {nmrdir} {outname} {rate=11025} {repeat=1} {concat=0}" << endl;
		cout << "  {nmrdir}...nmr Data Directory for your data" << endl;
		cout << "  {Rate}.....is default of 11025 Hz" << endl;
		cout << "  {repeat}...the FID 'repeat' times in the out file" << endl;
		cout << "  {concat}...if the file is a 2D, the a non 0 here " << endl;
		cout << "             will attach all the fids in one file " << endl;
	}else {

		std::string inF=argv[1];
		matrix FIDs;
		Vector<complex> FID;

		bool is2D;
		XWINNMRstream xwin;
		VNMRdata vdat;
		SpinSightStream sdat;


		if(!ReadData(xwin, inF, FIDs, FID, is2D)){
			if(!ReadData(vdat, inF, FIDs, FID, is2D)){
				if(!ReadData(sdat, inF, FIDs, FID, is2D)){
					std::cerr<<"Error: Cannot open any FIDs..."<<std::endl;
					return 0;
				}
			}
		}

		std::string waveout=argv[2];
		if(waveout.find(".wav") < waveout.size()){
			waveout=waveout.substr(0,waveout.size()-4);
		}
		int samples=11025;
		if(argc>=4){	samples=atoi(argv[3]);	}
		int repeat=1;
		if(argc>=5){	repeat=atoi(argv[4]);	}
		int concat=0;
		if(argc>=6){	concat=atoi(argv[5]);	}


		matstream matOut("matv.mat", ios::out);
		if(!is2D){
			wavestream myW;
			waveout+=".wav";
			myW.open(waveout.c_str(), ios::out);
			myW.setupFormat(samples, 16, 2); //samples, bits, channels
			FID*=1.0/max(abs(FID));
			int i=0;
			while(i++<repeat)	myW<<FID;
			myW.close();
			matOut.put("vdat", FID);
		}else{
			if(concat!=0){
				wavestream myW;
				std::string tmF=waveout+".wav";
				myW.open(tmF.c_str(), ios::out);
				myW.setupFormat(samples, 16, 2); //samples, bits, channels
				if(FIDs.cols()>FIDs.rows()){
					for(int i=0;i<FIDs.rows();++i)
					{
						FID=FIDs.row(i);
						FID*=1.0/max(abs(FID));
						int i=0;
						while(i++<repeat)	myW<<FID;
					}
				}else{
					for(int i=0;i<FIDs.cols();++i)
					{
						FID=FIDs.col(i);
						FID*=1.0/max(abs(FID));
						int i=0;
						while(i++<repeat)	myW<<FID;
					}
				}
				myW.close();
			}else{
				if(FIDs.cols()>FIDs.rows()){
					for(int i=0;i<FIDs.rows();++i)
					{
						FID=FIDs.row(i);
						FID*=1.0/max(abs(FID));
						std::string tmF=waveout+itost(i)+".wav";
						wavestream myW;
						myW.open(tmF.c_str(), ios::out);
						myW.setupFormat(samples, 16, 2); //samples, bits, channels
						int i=0;
						while(i++<repeat)	myW<<FID;
						myW.close();
					}
				}else{
					for(int i=0;i<FIDs.cols();++i)
					{
						FID=FIDs.col(i);
						FID*=1.0/max(abs(FID));
						std::string tmF=waveout+itost(i)+".wav";
						wavestream myW;
						myW.open(tmF.c_str(), ios::out);
						myW.setupFormat(samples, 16, 2); //samples, bits, channels
						int i=0;
						while(i++<repeat)	myW<<FID;
						myW.close();
					}

				}
			}
			matOut.put("vdat", FIDs);
		}
		matOut.close();
	}
	return 0;
}
