
#include "RecoupleContent.h"
#include "spindex.h"
#include "blochlib.h"


//the required 2 namespaces
using namespace BlochLib;
using namespace std;

int main(int argc,char* argv[])
{
	int q=1, len;
	std::string fname;

	MPIworld.start(argc, argv);
	query_parameter(argc, argv, q++,"Enter File Name to parse ", fname);
	if(MPIworld.master()) std::cout<<std::endl<<"Parameter File "<<fname<<std::endl;
	MPIworld.scatter(fname);

	Parameters pset(fname);
	pset.addSection("spins");
	pset.addSection("params");
	pset.addSection("pulses");
	pset.addSection("permutate");

//our main container for the trains and associated data
	RecoupleContent chunk;

	len=pset.getParamI("permutations", "permutate");
	if(len%2 !=0 && len!=1){ cerr<<" the permutations must be an even number..."<<endl; exit(0);	}

	chunk.maxtstep=pset.getParamD("maxtimestep", "params");

//load in the pregenerated tensor file
	std::string tensorname=pset.getParamS("tensorfile", "params", false);
	Vector<std::string> fcut=parse_param(tensorname, ',');
	if(fcut.size()<=0){	cerr<<" Need one or more tensor files ..."<<endl; exit(0);	}
	chunk.tensors.resize(fcut.size());
	for(int i=0;i<fcut.size();++i){
		fstream intf(fcut[i].c_str(), ios::binary | ios::in);
		if(intf.fail()){  cerr<<" the tensor file cannot be opened..."<<endl; exit(0);	}
		intf>>chunk.tensors[i];
		intf.close();
	}

//powder averageing container
	std::string powname=pset.getParamS("powder", "params");
	int thetast=1, phist=1;
	if(powname!="liquid"){
		thetast=pset.getParamI("thetasteps", "params");
		phist=pset.getParamI("phisteps", "params");
		chunk.pows=powder(powname, thetast, phist);

	}else{
		double oneth=pset.getParamD("singletheta", "params", false, 0.0)*PI/180.0;
		double oneph=pset.getParamD("singlephi", "params", false, 0.0)*PI/180.0;
		double onega=pset.getParamD("singlegamma", "params", false, 0.0)*PI/180.0;
		chunk.pows=powder(oneth, oneph, onega);
	}

//the spin system
	chunk.sys=SolidSys(pset.section("spins"));

	if(MPIworld.master()){ //print once if use MPI
		std::cout<<chunk.sys<<std::endl;
		std::cout<<"Running on : "<<MPIworld.size()<<" Processrs..."<<std::endl;
	}

//get the 'sub unit' info and parse it up
	chunk.SubUnits.setParams(pset.getParamS("subUnitIds", "permutate"),
				  pset.getParamS("typeIds", "permutate"),
				  pset.getParamS("units", "permutate"));

	chunk.SubUnits.debugflag=pset.getParamI("debuglevel", "params", false, 0);
	chunk.progress=(pset.getParamI("progress", "params", false, 0));

//see if we are caculating anything, or just a generate
// ans save operations
	std::string genOps=pset.getParamS("generationOps", "params", false,"GenerateSaveAndRun");
//file to dump out (or Read) the data
	std::string fsave=pset.getParamS("trainFile", "params", false, "rctrain.bin");

//this is the pulses sub parameters section
	Parameters myPulses(pset.section("pulses"));

//read the items in if that is what we do
//generate all the pulse trains in both the
// Sub Units AND the master (note:: these are simply
// integers and not large data strucutres so memory
// should not be a problem)
	if(genOps=="ReadAndRun"){
		chunk.wr=chunk.SubUnits.generateTrains(myPulses);
		chunk.read(fsave);
	}else if(genOps=="GenerateAndRun"){
		chunk.generateTrains(len,myPulses);
	}else{
		chunk.generateTrains(len,myPulses, fsave);
	}
//calculating the traces options

	int napps=pset.getParamI("napps", "permutate", false, 1);
	int tracidx=pset.getParamI("traceForIndex", "params", false, -1);
	std::string datout=pset.getParamS("traceout", "params", false, "");

//serveral aux file output names
	std::string seqout=pset.getParamS("seqout", "params", false, "seqnames.m");
	std::string namef=pset.getParamS("nameout", "params", false, "names.m");


//all the trace desired...
	if(tracidx==-1 && genOps!="GenerateAndSave"){
		chunk.calcTrace(napps);

	//dump our a 'logf' of what was done (NOTE:: this file can get huge)
		if(MPIworld.master()){ //print out the data for ONE proc only..

		//dump out the traces data
			int ct=0;
			std::string matfname=pset.getParamS("matout", "params", false, "matout.mat");
			matstream matout(matfname, ios::out | ios::binary);
			matout.put("data", chunk.traces);
			Vector<int> tsize(chunk.tensorSize(),0);
			Vector<int> subdiv(1,0);

			int hold=chunk.tensors[0].size(0);
			ct=0;
			for(int kk=0;kk<chunk.tensors.size();++kk){
				for(int i=0;i<chunk.tensors[kk].size();++i){
					tsize[ct]=chunk.tensors[kk].size(i);	//fill a vector where each element corresponds to an order
					if(hold!=tsize[ct]){
						subdiv.push_back(ct); //get the index where the order changes
						hold=tsize[ct];
						++ct;
					}
				}
			}
			subdiv.push_back(tsize.size()-1);
			matout.put("tsize", tsize);

		//add all the same orders together and place them inside a
		//new matrix such that we can display them
			matrix sumorders(subdiv.size()-1, chunk.trains.size());
			for(int j=0;j<chunk.trains.size();++j){
				Vector<complex> tmpc=chunk.traces.col(j);
				for(int i=0;i<subdiv.size()-1;++i){
					sumorders(i,j)=sum(tmpc(Range(subdiv(i), subdiv(i+1))))/(subdiv(i+1)-subdiv(i));
				}
			}
			matout.put("orders", sumorders);
			matout.put("subdiv", subdiv);
			matout.close();

			ofstream tout;
			if(datout!=""){
				tout.open(datout.c_str());

				tout<<std::string(50, ' ');
				//the header
				std::string tmp;
				for(int i=0;i<chunk.trains.size();++i)
				{
					tmp="Train "+itost(i);
					tout<<tmp;
					if(tmp.size()<20){	tout<<string(20-tmp.size(), ' ');	}
				}
				tout<<endl;

				for(int kk=0;kk<chunk.tensors.size();++kk){
					for(int i=0;i<chunk.tensors[kk].size();++i)
					{
						tmp=chunk.tensors[kk].name(chunk.sys,i);
						tout<<tmp;
						//tout<<endl<<chunk.tensors[kk].tensor(chunk.sys, i)<<endl;

						if(tmp.size()<50){	tout<<std::string(50-tmp.size(), ' ');	}
						for(int j=0;j<chunk.traces.cols();++j)
						{
							tmp=dbtost_form("%5.4f", Re(chunk.traces(ct,j)));
							tmp+=" ";
							if(Im(chunk.traces(i,j))>0.0) tmp+="+";
							tmp+=dbtost_form("%5.4f", Im(chunk.traces(ct,j)));
							tmp+="*i";
							tout<<tmp;
							if(tmp.size()<20){ tout<<string(20-tmp.size(), ' ');	}
						}
						++ct;
					}
					tout<<endl;
				}
			}

		//dump out the sequence names...

			ofstream seqfout(seqout.c_str());
			std::string nm=seqout,tmp;
			seqout=seqout.substr(seqout.rfind("/")+1, seqout.size());
			if(seqout.find(".m")<seqout.size())
			{	nm=seqout.substr(0,seqout.find(".m"));		}

			seqfout<<"function h="<<nm<<"()"<<endl;
			seqfout<<"h=[";

			for(int i=0;i<chunk.trains.size();++i)
			{
				tmp="Train "+itost(i);
				if(datout!="")
				{
					tout<<endl;
					tout<<tmp<<endl<<chunk.trains[i]<<endl;
				}

				std::string seqtmp=chunk.sequenceName(i);
				seqfout<<"\t'("<<seqtmp<<")"<<napps<<"'"<<endl;
			}

			seqfout<<"];"<<std::endl;
			seqfout<<"return;"<<std::endl;

		//dump out the tensor names strucutre
			fstream nameout(namef.c_str(),ios::out);
			namef=namef.substr(namef.rfind("/")+1, namef.size());

			Vector<std::string> namelist(chunk.tensorSize(),"");
			int maxlen=0;
			ct=0;
			for(int i=0;i<chunk.tensors.size();++i){
				for(int kk=0;kk<chunk.tensors[i].size();++kk){
					namelist[ct]=chunk.tensors[i].name(chunk.sys,kk);
					maxlen=std::max(maxlen, int(namelist[ct].size()));
					++ct;
				}
			}

			nm=namef;
			if(namef.find(".m")<namef.size())
			{	nm=namef.substr(0,namef.find(".m"));		}

			nameout<<"function h="<<nm<<"()"<<endl;
			nameout<<"h=[";

			for(int i=0;i<namelist.size();++i){
				nameout<<"\t'"<<std::string(maxlen-namelist[i].size(), ' ');
				nameout<<namelist[i]<<"'"<<endl;

			}
			nameout<<"];"<<std::endl;;
			nameout<<"return;"<<std::endl;
		} //end MPIrank
	}else if(tracidx!=-2 && genOps!="GenerateAndSave"){
		Vector<complex> trac=chunk.calcTrace(tracidx,napps);
		//the header
		std::string tmp;
		tmp="Train "+itost(tracidx);
		int ct=0;
		if(datout!="")
		{
			std::ofstream tout(datout.c_str());
			tout<<std::string(50, ' ');
			tout<<tmp;
			if(tmp.size()<20){	tout<<std::string(20-tmp.size(), ' ');	}
			tout<<endl;


			for(int kk=0;kk<chunk.tensors.size();++kk){
				for(int i=0;i<trac.size();++i)
				{
					tmp=chunk.tensors[kk].name(chunk.sys,i);
					tout<<tmp;
					if(tmp.size()<50){	tout<<std::string(50-tmp.size(), ' ');	}
					tmp=dbtost_form("%5.4f", Re(trac(ct)));
					tmp+=" ";
					if(Im(trac(i))>0.0) tmp+="+";
					tmp+=dbtost_form("%5.4f", Im(trac(ct)));
					tmp+="*i";
					tout<<tmp;
					if(tmp.size()<20){ tout<<string(20-tmp.size(), ' ');	}
					tout<<endl;
					++ct;
				}
			}
			tout<<endl<<endl;
			tmp="Train "+itost(tracidx);
			tout<<tmp<<endl<<chunk.trains[tracidx]<<endl;
		}
		matstream matout("matout.mat", ios::out | ios::binary);
		matout.put("data", trac);
		Vector<int> tsize(chunk.tensorSize(), 0);
		for(int kk=0;kk<chunk.tensors.size();++kk){
			for(int i=0;i<chunk.tensors.size();++i)
			{	tsize[ct]=chunk.tensors[kk].size(i);	++ct; }
		}
		matout.put("tsize", tsize);
		matout.close();
	}


//do FID collection if wanted
	pset.addSection("fids");
	int fidval=pset.getParamI("fidForIndex", "fids", false, -1);
	if(fidval!=-1  && genOps!="GenerateAndSave"){
		int npts=pset.getParamI("npts", "fids");
		std::string ro=pset.getParamS("initro", "fids", false, "Iz");
		std::string det=pset.getParamS("detect", "fids", false, "Iz");
		std::string fout=pset.getParamS("fidout", "fids", false, "dat");

		HamiltonianGen hams;
		matrix roeq=hams.Hamiltonian(chunk.sys, ro);
		matrix detect=hams.Hamiltonian(chunk.sys, det);

		Vector<complex> fid=chunk.FID(fidval, roeq, detect, npts);
		if(MPIworld.master())
			plotterFID(fid, fout, chunk.SubUnits.subPropTime[0]);
	}

//do coherence transfers
	pset.addSection("transfer");
	int dotrans=pset.getParamI("doTransfer", "transfer", false, -1);
	if(dotrans>0  && genOps!="GenerateAndSave"){
		std::string ro=pset.getParamS("initro", "transfer", false, "Iz");
		std::string det=pset.getParamS("detect", "transfer", false, "Iz");
		std::string fout=pset.getParamS("transOut", "transfer", false, "trans");

		HamiltonianGen hams;
		matrix roeq=hams.Hamiltonian(chunk.sys, ro);
		matrix detect=hams.Hamiltonian(chunk.sys, det);

	//get the trains
		Vector<std::string> TransTrains=
			parse_param(pset.getParamS("transList", "transfer"), ',');
	//get the 'best' apps
		Vector<int> Apps=pset.getParamVectorI("napps", "transfer");

		if(Apps.size()!=TransTrains.size() && MPIworld.master()){
			BLEXCEPTION(" 'napps' must be the same size as 'transList'")
		}

	//get the 1-D run through
		Vector<std::string> stC0iso=parse_param(pset.getParamS("alter1D", "transfer", false), ':');
	//get the 2-D run through
		Vector<std::string> stC1iso=parse_param(pset.getParamS("alter2D", "transfer", false), ':');
		double stC0=0.0, enC0=0.0,stpC0=0.0, curC0=0.0,
		       stC1=0.0, enC1=0.0, stpC1=0.0,curC1=0.0;
		int oDsize(0), tDsize(0);
		if(stC0iso.size()==4){
			stC0=std::atof(stC0iso[1].c_str());
			curC0=stC0;
			enC0=std::atof(stC0iso[3].c_str());
			stpC0=std::atof(stC0iso[2].c_str());
			oDsize=int(std::ceil((enC0-stC0)/stpC0)+1);
			if(MPIworld.master()){
			    std::cout<<" stepping though "<<stC0iso[0]<<"-->"
			         <<stC0<<":"<<stpC0<<":"<<enC0<<std::endl;
			}
		}

		if(stC1iso.size()==4){
			stC1=std::atof(stC0iso[1].c_str());
			curC1=stC1;
			enC1=std::atof(stC0iso[3].c_str());
			stpC1=std::atof(stC0iso[2].c_str());
			tDsize=int(std::ceil((enC1-stC1)/stpC1)+1);
			if(MPIworld.master()){
				std::cout<<" stepping though "<<stC1iso[0]<<"-->"
			         <<stC1<<":"<<stpC1<<":"<<enC1<<std::endl;
		   }
		}


	//sub unit defs
		std::string subUU=pset.getParamS("subUnitIds", "transfer");
		std::string typeIDS=pset.getParamS("typeIds", "transfer");

		Vector<complex> trans(TransTrains.size(),0);

		std::ofstream fOut;
		matstream matOut;
		for(int i=0;i<TransTrains.size();++i){
			if(MPIworld.master() && i==0 &&
			  (oDsize<=0 || tDsize<=0)){
				fOut.open(fout.c_str());
			}
		//'sub unit' info and parse  up
			chunk.SubUnits.setParams(subUU, typeIDS,TransTrains[i]);
			chunk.SubUnits.generateTrains(myPulses);
			if(oDsize<=0 || tDsize<=0){
				trans[i]=chunk.transfer(roeq, detect, Apps[i]);
				if(MPIworld.master()){
					fOut<<Apps[i]<<" "<<trans[i].Re()<<" "<<trans[i].Im()<<std::endl;
				}
			}else{
				std::string matoutname=fout+itost(i);
				matrix trans2D(oDsize, tDsize, 0.0);
				_matrix<int, FullMatrix> napps2D(oDsize, tDsize, 0.0);
				int r=0, c=0;
				for(curC0=stC0; curC0<=enC0; curC0+=stpC0){
					chunk.sys.setSpinParam(stC0iso[0], curC0);
					c=0;
					for(curC1=stC1; curC1<=enC1; curC1+=stpC1){
						chunk.sys.setSpinParam(stC1iso[0], curC1);
						if(chunk.progress>1 & MPIworld.master()){ cout<<chunk.sys; }
						trans2D(r,c)=chunk.transfer(roeq, detect, Apps[i]);
						c++;
					}
					r++;
				}
				if(MPIworld.master()){
					matOut.open(matoutname.c_str());
					matOut.put("vdat",trans2D);
					//matOut.put("napps",napps2D);
					matOut.close();
				}
			}
		}
		std::cout<<std::endl;
	}


	MPIworld.end();
	return 0;
}
