//#include "fullbinarytree.h"
#include <string>
#include "blochlib.h"
//#include "RecoupleSequence.h"
#include "RecoupleContent.h"
#include "RecoupleSubUnits.h"
#include "spindex.h"


//the required 2 namespaces
using namespace BlochLib;
using namespace std;



/****************** RECOUPLE CONTENT ****************/
int RecoupleContent::tensorSize() const
{
	int tsize=0;
	for(int i=0;i<tensors.size();++i) tsize+=tensors[i].size();
	return tsize;
}

/*
//len= the legnth of the permutation train
void RecoupleContent::generateTrains(int len, std::string type, int baseSym, int factorSym, double baseamp)
{
	//for the sequnce we only permute the phase by 0 or 180..
	//also we know, by symmetry, that the seqeunces that starts at
	//a phase =0 will give the same result if we bared the sequence
	// hense we only need to travel one side of the binary tree
	FullBinaryTree<0, 180> Perms(len-1);
	FullBinaryTree<0, 180>::iterator paths(Perms);

	//need to allow a single sequence for comparison purposes
	if(len==1){
		Vector<Recouple> pc7s(1,Recouple(type, baseSym,factorSym, baseamp, 0.0, 0.0));
		trains.push_back(RecoupleTrain(pc7s));
		return;
	}

	int maxsize=200, onsize=maxsize, ct=0;
	trains.resize(maxsize);
	Vector<RecoupleTrain> holdSubs=SubUnits.generateTrains(type,baseSym,factorSym, baseamp);
	//cout<<"HOLD SUB"<<endl<<holdSubs<<endl;
	Vector<int> tmpList(maxsize);

	while(paths){
		FullBinaryTree<0,180>::nodeiterator nodes(paths);
		int totLength=0;
		int num180=0, num0=0;
		tmpList.resize(0);
		while(nodes){
			if(nodes.value()==180){
				tmpList.push_back(1);
				totLength+=holdSubs[1].size();
			}else{
				tmpList.push_back(0);
				totLength+=holdSubs[0].size();
			}
			if(nodes.value()==180) ++num180;
			else ++num0;
			++nodes;
		}
		++paths;
		if(num180==num0){
			//to save speed on the nasty push_back operation
			if(ct>=onsize){ onsize+=maxsize; trains.resizeAndPreserve(onsize); }
			int rlen=0;
			Vector<Recouple> holdR(totLength);
			double onTime=0.0;

			for(int k=0;k<tmpList.size();++k){
				holdSubs[tmpList[k]].setBeginTime(onTime);
				for(int pp=0;pp<holdSubs[tmpList[k]].size();++pp){
					holdR[rlen]=holdSubs[tmpList[k]][pp];
					++rlen;
				}
				onTime=holdSubs[tmpList[k]].endTime();
			}
			trains[ct]=RecoupleTrain(holdR);
			++ct;
		}
	}
	trains.resizeAndPreserve(ct);
}
*/





//generates the Next K Subset of a generic list of
// N numbers (0, 1,2,3...) where we only want
// a sub set of n numbers...of course there are many posibilities
//if 'more'  is false, the list initializes, if it is true it
// calcs the next subset..this function is quite generic
void RecoupleContent::nextKSubset(int n, int k, int *subset,bool &more)
{
	static int m2, m;
	if(n<0 || k<0) return;
	if( !more ){
		m2 = 0;
		m = k;
	}else{
		if( m2 < n-m )  m = 0;
		++m;
		m2 = subset[k-m];
	}
	for(int i=0;i<m;++i) subset[k+i-m] = m2 + i+1;
	more = subset[0] != n-k+1;
}

//this calulates the 'length' of the permutation
// RecoupleSubUnits generate 'n' Propogators, but we
// typically only what N combinations of them
// to include all posible permutations the Minim length
// of the vector can only be mod SubUnits.size()
int RecoupleContent::permutationLength(int order)
{
	if(order<=0){
		std::cerr<<std::endl<<"Error: RecoupleContent::permutationLength "<<std::endl;
		std::cerr<<" Order is 0 or smaller...this cannot be.."<<std::endl;
		exit(1);
	}

	int base=SubUnits.size(), ct=base;

//for ex: if we have 5 subunits, and want to use only 4
// the permutation length is the length of the subunits.
	if(base>order) return base;

	while(base<order) base+=ct;
	return base;

}

//this fills up the 'index' list out to the appropriate length
// as given by 'permutationLength(order)'  it simply
// returns a list with int Mod SubUnits.size()
void RecoupleContent::generateLists(int permLen, int *propList)
{
	if(permLen<=0){
		std::cerr<<std::endl<<"Error: RecoupleContent::generateLists "<<std::endl;
		std::cerr<<" permutation Length is 0 or smaller...this cannot be.."<<std::endl;
		exit(1);
	}

	int ct=0;
	while(ct<permLen){ propList[ct]=ct%SubUnits.size(); ++ct; }
}

//generates the 'master' trains and the SubUnit trains
void RecoupleContent::generateTrains(int order,Parameters &in)
{
	wr=SubUnits.generateTrains(in);
	generateTrains(order);
}

//generates the 'master' trains and the SubUnit trains
//will either READ the data from a file of the 'ReadAndRun' flag
// is set, OR it will set the flag to dump
// out data when generated
void RecoupleContent::generateTrains(int order,Parameters &in, std::string outF)
{
	wr=SubUnits.generateTrains(in);
	generateTrains(order, outF);
}

//generates the 'master' trains and the SubUnit trains
//will either READ the data from a file of the 'ReadAndRun' flag
// is set, OR it will set the flag to dump
// out data when generated
void RecoupleContent::generateTrains(int order,Parameters &in, std::fstream &outF)
{
	wr=SubUnits.generateTrains(in);
	generateTrains(order);
}

//will either READ the data from a file of the 'ReadAndRun' flag
// is set, OR it will set the flag to dump
// out data when generated
void RecoupleContent::generateTrains(int order, std::string outF)
{
	outFileName=outF;
	//if(outFile)outFile->is_open()) outFile->close();
	if(killFile){
		outFile->close();
		delete outFile;
	}
	killFile=true;
	outFile=new std::fstream(outF.c_str(), ios::out | ios::binary);
	generateTrains(order, *outFile);
}

//will either READ the data from a file of the 'ReadAndRun' flag
// is set, OR it will set the flag to dump
// out data when generated
void RecoupleContent::generateTrains(int order, std::fstream &outF)
{
	if(outF.fail()){
		std::cerr<<std::endl<<"Error: RecoupleContent::generateTrains "<<std::endl;
		std::cerr<<" OutPut File cannot be written to...."<<std::endl;
		exit(1);
	}
	if(killFile && &outF!=outFile){
		outFile->close();
		delete outFile;
		killFile=false;
	}else if(&outF!=outFile){
		outFile=&outF;
	}
	generateTrains(order);

}

//will either READ the data from a file of the 'ReadAndRun' flag
// is set, OR it will set the flag to dump
// out data when generated
void RecoupleContent::generateTrains(std::string outF)
{
	read(outF);
}

//will attempt to read the train data from a file
//setting the 'read and run' flag
void RecoupleContent::generateTrains(std::fstream &outF)
{
	read(outF);
}

void RecoupleContent::generateTrains(int order)
{
	int permlength=permutationLength(order);
	int *theperms=new int[permlength];
	generateLists(permlength,theperms);

	int maxlen=10000;
	trains.resize(maxlen);

	int curorder=order;
	int mastCt=0;
	Vector<int> curp(curorder,0);

//this is the k-subset loop;
// becuase some k-subsets can be dupilcated via the 'nextksubset'
// algorithm, we need to generate all of the k-subs, and
// remove all the dupes
	Vector<Vector<int> > Ksubs;
	Vector<int> curK(curorder);
	int *subset=new int[curorder];
	bool more=false, useMe=true;
	int onK=0;
	int *perms;
	do{
		nextKSubset(permlength, curorder, subset, more);
		for(int i=0;i<curorder;++i)	curK[i]=theperms[subset[i]-1];
		int j;
		std::sort(curK.data(), curK.data()+curorder);

	//check if the current one is a dupe...
		for(j=0;j<Ksubs.size();++j){
			useMe=true;
			if(curK==Ksubs[j]){ useMe=false; break;	}
		}

	//put in the list if valid
		if(useMe || onK==0){
			Ksubs.push_back(curK);
			++onK;
		}
	}while(more); //end subset loop


//this is the permutation on the k-subset loop
	onK=0;
	int Kstart=0;
	if(SubUnits.debugflag>=1){	std::cout<<std::endl;	}
	while(onK<Ksubs.size()){
		perms=Ksubs[onK].data();
		Kstart=mastCt;
		if(SubUnits.debugflag>=1){
			std::cout<<"On the Kth subset: "<<onK<<"/"
					<<Ksubs.size()<<" subset="
					<<Ksubs[onK]<<std::endl;
		}
		do{
			for(int i=0;i<curorder;++i)	curp[i]=theperms[perms[i]];

			//need to see if this entry is just a 'mirror' of
			//an already entered sequence...if so, we do not include it
			Vector<int> revC(curorder);
			for(int i=0, j=curorder-1;i<curorder;++i){
				revC[i]=curp[j]; --j;
			}

			useMe=true;
			for(int i=Kstart;i<trains.size();++i){
				if(revC==trains[i] || curp==trains[i]){ useMe=false; break;	}
			}

			if(useMe){
				if(mastCt>=maxlen){ maxlen+=10000; trains.resizeAndPreserve(maxlen);	}
				trains[mastCt]=curp;
				if(SubUnits.debugflag>=2){
					std::cout<<"RecoupleContent info--train: "<<mastCt<<" Order: "<<curorder<<" Perms: ";
					std::cout<<sequenceName(mastCt)<<std::endl;
				}

				++mastCt;
			}

		}while(next_permutation(perms, perms+curorder)); //end the permutation loop
		++onK;
		//cout<<endl;
		//the second half of the list is
	}
	//resize the list to the prop size
	trains.resizeAndPreserve(mastCt);

	//close out the outFile
	if(outFile){	write(*outFile);	}

	delete [] subset;
	delete [] theperms;
}

//this performs the evolutions to one train in a vector of trains
// it returns the final propogator
// it requiers a powder angle
matrix RecoupleContent::generateProp(int which, powder::iterator &powit)
{
	RunTimeAssert(which<trains.size());
	matrix U=sys.Fe();
	if(SubUnits.debugflag>=2){
		std::cout<<"RecoupleContent info--current train: "<<which<<" Perms: ";
		std::cout<<sequenceName(which)<<std::endl;
	}
	for(int i=0;i<trains[which].size();++i){
		//cout<<endl<<trains[which]<<" "<<i<<endl;
		U*=SubUnits.Props[trains[which][i]]*U;
	}

	return U;
}

//gets the name of the squence 'which'as generated
// from the units of 'SubUnits'
std::string RecoupleContent::sequenceName(int which)
{
	RunTimeAssert(which<trains.size());

	std::string name;
	for(int i=0;i<trains[which].size();++i){
		name+=SubUnits.units[trains[which][i]];
	}
	return name;
}

//generate an FID from the train 'which'
// the slave runner for MPI slave
Vector<complex> RecoupleContent::FID(int which, matrix &roeq, matrix &detect,int npts, int napps, int powpt)
{
	Vector<complex> fid(npts, 0.0);
	sys.setPowderAngles(pows.theta(powpt), pows.phi(powpt), pows.gamma(powpt));
	SubUnits.generateProps(
		sys,	// the Hamiltonian generator
		wr,	    //rotor speed
		maxtstep	// maximal time step allowed
	);


	matrix U=sys.Fe();
	static complex norm=trace(detect, adjoint(detect));
	fid(0)+=trace(roeq, detect)/norm*pows.weight(powpt);
	int fidct=1;
	while(fidct<npts){
		for(int i=0;i<trains[which].size();++i){
			U*=SubUnits.Props[trains[which][i]];
			fid(fidct)+=trace(prop(U, roeq), detect)/norm*pows.weight(powpt);
			++fidct;
			if(fidct>=npts) break;
		}
	}
	return fid;

}
//generates an fid from a specific Recouple Train
Vector<complex> RecoupleContent::FID(int which, matrix &roeq, matrix &detect,int npts, int napps)
{
	Vector<complex> fid(npts, 0.0);
	if(MPIworld.master()){
		std::cout<<"Number of trains: "<<trains.size()<<", Order of trains: "<<trains[which].size()<<std::endl;
	}
	if(MPIworld.master() && which>=trains.size()){
		std::cerr<<endl<<" FID calculation for the "<<which<<" pulse train"<<std::endl;
		std::cerr<<" cannot be calculated because the train does not exsist"<<std::endl;
		return(fid);
	}
	if(MPIworld.master() && progress>=1){
		std::cout<<endl<<" FID calculation for the "<<which<<" pulse train"<<std::endl;
		std::cout<<"calculating FID for train: "<<sequenceName(which)<<std::endl;
		std::cout<<"\t----powder pt-------(theta,phi)"<<std::endl;
	}

	int done=-1,pos;
	if(!MPIworld.serial()){
		if(MPIworld.master() && pows.size()>1)
		{
			int on=0;
			for(int k=1;k<MPIworld.size();++k)
			{
				MPIworld.put(on, k);
				++on;
				if(on>=pows.size()) break;
			}
			int get, dd;
			for(int k=on;k<pows.size();++k)
			{
				get=MPIworld.getAny(dd);
				if(progress>=1)
				{
					std::cout<<"\t----"<<itost_form("%4d", k)
						 <<"/"<<itost_form("%4d", pows.size())
						 <<"-------("<<dbtost_form("%1.2f", pows.theta(k))
						 <<","<<dbtost_form("%1.2f", pows.phi(k))<<")"
						 <<"\r";
						 std::cout.flush();

				}
				MPIworld.put(k, get);
			}
			for(int k=1;k<MPIworld.size();++k)
				MPIworld.put(done, k);
		}else{ //SLAVE
			while(1)
			{
				MPIworld.get(pos, 0);
				if(pos==done) break;
				fid+=FID(which, roeq, detect, npts, napps, pos);
				MPIworld.put(pos, 0);
			}
		}
	}else{ //SERIAL
		for(int i=0;i<pows.size();++i){
			if(progress>=1)
			{
				std::cout<<"\t----"<<itost_form("%4d", i)
					 <<"/"<<itost_form("%4d", pows.size())
					 <<"-------("<<dbtost_form("%1.2f", pows.theta(i))
					 <<","<<dbtost_form("%1.2f", pows.phi(i))<<")"
					 <<"\r";
					 std::cout.flush();

			}
			fid+=FID(which, roeq, detect, npts, napps, i);
		}

	}

	MPIworld.reduce(fid, Reduce::Add);

	return fid;
}

//generate a 2Q transfer profile
// the master runner for MPI master
complex RecoupleContent::transfer(matrix &roeq, matrix &detect, int &napps)
{
	complex trans(0.0,0.0);

	if(MPIworld.master()){
		std::cout<<"Number of trains: "<<trains.size()<<", Order of trains: "<<trains[0].size()<<std::endl;
	}

	if(MPIworld.master() && progress>=1){
		std::cout<<endl<<" Transfer calculation "<<std::endl;
		std::cout<<"calculating Transfer for train: "<<sequenceName(0)<<std::endl;
		std::cout<<"\t----powder pt-------(theta,phi)"<<std::endl;
	}

	int done=-1,pos;
	if(!MPIworld.serial()){
		if(MPIworld.master() && pows.size()>1)
		{
			int on=0;
			for(int k=1;k<MPIworld.size();++k)
			{
				MPIworld.put(on, k);
				++on;
				if(on>=pows.size()) break;
			}
			int get, dd;
			for(int k=on;k<pows.size();++k)
			{
				get=MPIworld.getAny(dd);
				if(progress>=1)
				{
					std::cout<<"\t----"<<itost_form("%4d", k)
						 <<"/"<<itost_form("%4d", pows.size())
						 <<"-------("<<dbtost_form("%1.2f", pows.theta(k))
						 <<","<<dbtost_form("%1.2f", pows.phi(k))<<")"
						 <<"\r";
						 std::cout.flush();

				}
				MPIworld.put(k, get);
			}
			for(int k=1;k<MPIworld.size();++k)
				MPIworld.put(done, k);
		}else{ //SLAVE
			while(1)
			{
				MPIworld.get(pos, 0);
				if(pos==done) break;
				trans+=transfer(roeq, detect, napps, pos);
				MPIworld.put(pos, 0);
			}
		}
	}else{ //SERIAL
		for(int i=0;i<pows.size();++i){
			if(progress>=1)
			{
				std::cout<<"\t----"<<itost_form("%4d", i)
					 <<"/"<<itost_form("%4d", pows.size())
					 <<"-------("<<dbtost_form("%1.2f", pows.theta(i))
					 <<","<<dbtost_form("%1.2f", pows.phi(i))<<")"
					 <<"\r";
					 std::cout.flush();

			}
			trans+=transfer(roeq, detect, napps, i);
		}

	}
	MPIworld.reduce(trans, Reduce::Add);

	return trans;
}

//generate an 2Q from the train 'which'
// the slave runner for MPI slave
// it finds THE MAX transfer
// the max transfer occurs when the first 'dip' occurs
complex RecoupleContent::transfer(matrix &roeq, matrix &detect,int &napps, int powpt)
{
	complex trans(0.0,0.0), oldtrans=trans;
	sys.setPowderAngles(pows.theta(powpt), pows.phi(powpt), pows.gamma(powpt));
	SubUnits.generateProps(
		sys,	// the Hamiltonian generator
		wr,	    //rotor speed
		maxtstep	// maximal time step allowed
	);


	matrix U=sys.Fe(), tmU=U;
	int fidct=1;
	for(int i=0;i<trains[0].size();++i){
		tmU*=SubUnits.Props[trains[0][i]];
	}

/*	if(powpt==0){
		while(1){
			oldtrans=trace(prop(U, roeq), detect)/trace(detect, adjoint(detect))*pows.weight(powpt);
			if(abs(oldtrans)>abs(trans)){
				trans=oldtrans;
				U*=tmU;
				napps++;
			}else{
				return trans;
			}
		}

	}else{
*/
//		dmatrix Ps=exp(-sys.Fz()*complexi*pi/2.0);
//		dmatrix PsB=exp(sys.Fz()*complexi*pi/2.0);

		for(int i=0;i<napps;++i){		U*=tmU;}
//		U=Ps*U*PsB*U;
		return trace(prop(U, roeq), detect)/trace(detect, adjoint(detect))*pows.weight(powpt);
	//}
}


double RecoupleContent::sequenceTime(int which){
	RunTimeAssert(which<trains.size());
	double dt=0.0;
	for(int k=0;k<trains[which].size();++k){
		dt+=SubUnits.propTime[trains[which][k]];
	}
	return dt;
}


//calculate the tensorial components for the train 'which'
//napps is the number of sucsessive applications of the train (i.e. as in a 2D experiment)
Vector<complex> RecoupleContent::calcTrace(int which, int napps)
{

	int begin=0;
	int end=pows.size();
	int divs=1;

	Range powR=MPIworld.splitLoop(begin, end, divs);

	Vector<complex> trac(tensorSize(), 0.0);
	if(tensors.size()<=0) return trac;
	powder::iterator powit(pows, powR);
	matrix H;
	if(MPIworld.master())
		std::cout<<"\t---train   ----powder pt-------(theta,phi)"<<std::endl;


	int ct=0;
	while(powit){
		sys.setPowderAngles(powit.theta(), powit.phi(), powit.gamma());
		SubUnits.generateProps(
				sys,	// the Hamiltonian generator
				wr,	    //rotor speed
				maxtstep	// maximal time step allowed
		);

		std::cout<<"\t---"<<itost_form("%3d", which)<<"/"<<itost_form("%3d", trains.size())
				 <<"----"<<itost_form("%4d", powit.curpos())
				 <<"/"<<itost_form("%4d", pows.size())
				 <<"-------("<<dbtost_form("%1.2f", powit.theta())
				 <<","<<dbtost_form("%1.2f", powit.phi())<<")"
				 <<"\r";
		matrix U=generateProp(which, powit);
		double dt=sequenceTime(which);
		matrix tempU;
		if(napps>1) tempU=U;
		for(int i=1;i<napps;++i){	U=chop(tempU*U); }
		H=Mlog(U, 1.0/(complexi*dt*PI2));
		double nn=norm(trace(H,adjoint(H)));
		ct=0;
		for(int j=0;j<tensors.size();++j){
			for(int i=0;i<tensors[j].size();++i){
				U=tensors[j].tensor(sys, i);
				if(nn>0.0)	trac(ct)+=(norm(trace(H, U))/nn)*powit.weight();
				++ct;
			}
		}
		++powit;
	}

	MPIworld.reduce(trac, Reduce::Add);

	return trac;
}

//calculate the tensorial components for all the trains
//placing each one in the matrix 'traces'
//napps is the number of sucsessive applications of the train (i.e. as in a 2D experiment)
void RecoupleContent::calcTrace(int napps)
{
//	int begin=0;
//	int end=pows.size();
//	int divs=1;

//	Range powR=MPIworld.splitLoop(begin, end, divs);
	traces.resize(tensorSize(), trains.size());
	traces.fill(0.0);

//	Vector<complex> trac(tensorSize(), 0.0);
	if(tensorSize()<=0) return;
	if(MPIworld.master() && progress>=2)
		std::cout<<"\t---train  ----powder pt-------(theta,phi)"<<std::endl;
	if(MPIworld.master() && progress==1)
		std::cout<<"\t----powder pt-------(theta,phi)"<<std::endl;


	int ct=0, done=-1, PowTag=10;
	//matrix ro=sys.Fz();
	//matrix dete=sys.Fz();
	//int npts=128;
	//Vector<complex> fid(npts, 0.0);
	double dt=0;
	powder::iterator powit(pows);
	if(!MPIworld.serial()){
		if(MPIworld.master() && pows.size()>1)
		{
			int on=0;
			for(int k=1;k<MPIworld.size();++k)
			{
				on=powit.curpos();
				MPIworld.put(on, k, PowTag);
				++powit;
				if(!powit) break;
			}
			int get, dd;
			while(powit)
			{
				get=MPIworld.getAny(dd);
				if(progress==1)
				{
					std::cout<<"\t----"<<itost_form("%4d", powit.curpos())
					 <<"/"<<itost_form("%4d", pows.size())
					 <<"-------("<<dbtost_form("%1.2f", powit.theta())
					 <<","<<dbtost_form("%1.2f", powit.phi())<<")"
					 <<"\r";
					 std::cout.flush();
				}
				on=powit.curpos();
				MPIworld.put(on, get, PowTag);
				++powit;
			}
			for(int k=1;k<MPIworld.size();++k)
				MPIworld.put(done, k, PowTag);
		}else{
			int pos=0;
			matrix U, tempU,H;
			while(1)
			{
				MPIworld.get(pos, 0, PowTag);
				if(pos==-1) break;

				sys.setPowderAngles(pows.theta(pos), pows.phi(pos), pows.gamma(pos));
				SubUnits.generateProps(
						sys,	// the Hamiltonian generator
						wr,	    //rotor speed
						maxtstep	// maximal time step allowed
				);

				for(int which=0;which<trains.size();++which)
				{
					if(progress==2)
					{
						std::cout<<"\t---"<<itost_form("%3d", which+1)<<"/"<<itost_form("%3d", trains.size())
						 <<"----"<<itost_form("%4d", powit.curpos())
						 <<"/"<<itost_form("%4d", powit.size())
						 <<"-------("<<dbtost_form("%1.2f", powit.theta())
						 <<","<<dbtost_form("%1.2f", powit.phi())<<")"
						 <<"\r";
					}

					U=generateProp(which, powit);
					dt=sequenceTime(which);
					if(napps>1) tempU=U;
					for(int i=1;i<napps;++i){	U=tempU*U; }

					H=Mlog(U, 1.0/(complexi*dt*PI2));
					double nn=norm(trace(H,adjoint(H)));
					ct=0;
					for(int j=0;j<tensors.size();++j){
						for(int i=0;i<tensors[j].size();++i){
							U=tensors[j].tensor(sys, i);
							if(nn>0.0)	traces(ct,which)+=(norm(trace(H, U))/nn)*powit.weight();
							++ct;
						}
					}
				}
			/*	matrix U=generateProp(0, powit), tempU=U;
				dt=sequenceTime(0);
				for(int i=0;i<npts;++i){
					fid[i]+=trace(prop(U,ro), dete)*powit.weight();
					U=tempU*U;
				}
			*/
			MPIworld.put(pos, 0);
			}
		}
	}else{ //SERIAL
		matrix U, tempU,H;
		while(powit)
		{
			sys.setPowderAngles(powit.theta(), powit.phi(), powit.gamma());
				SubUnits.generateProps(
					sys,	// the Hamiltonian generator
					wr,	    //rotor speed
					maxtstep	// maximal time step allowed
			);

			for(int which=0;which<trains.size();++which)
			{
				if(progress>=2)
				{
					std::cout<<"\t---"<<itost_form("%3d", which+1)<<"/"<<itost_form("%3d", trains.size())
					 <<"----"<<itost_form("%4d", powit.curpos())
					 <<"/"<<itost_form("%4d", powit.size())
					 <<"-------("<<dbtost_form("%1.2f", powit.theta())
					 <<","<<dbtost_form("%1.2f", powit.phi())<<")"
					 <<"\r";
				}

				U=generateProp(which, powit);
				dt=sequenceTime(which);
				if(napps>1) tempU=U;
				for(int i=1;i<napps;++i){	U=tempU*U; }

				H=Mlog(U, 1.0/(complexi*dt*PI2));
				double nn=norm(trace(H,adjoint(H)));
				ct=0;
				for(int j=0;j<tensors.size();++j){
					for(int i=0;i<tensors[j].size();++i){
						U=tensors[j].tensor(sys, i);
						if(nn>0.0)	traces(ct,which)+=(norm(trace(H, U))/nn)*powit.weight();
						++ct;
					}
				}
			}
			++powit;
		}
	}
//	plotterFID(fid, "fid", dt);
//matrix lll=1e12*traces;
	MPIworld.reduce(traces, Reduce::Add);
}


//Reads in a list of 'trains' from a file
void RecoupleContent::read(std::string fname)
{
	outFileName=fname;
	if(outFile && killFile){
		outFile->close();
		delete outFile;
	}
	killFile=true;
	outFile =new std::fstream(fname.c_str(), ios::binary | ios::in);
	read(*outFile);
}

//Reads in a list of 'trains' from a file
void RecoupleContent::read(std::fstream &oo)
{

	char liner[100];
	if(oo.fail())
	{
		std::cerr<<std::endl<<"Error: RecoupleContent::read "<<std::endl;
		std::cerr<<" Cannot read from the file........"<<std::endl;
		exit(1);
	}
	while(!oo.eof()){
		oo.getline(liner, 100, '\n');
		if(std::string(liner)=="START RecoupleContent")	break;
	}

	if(oo.eof()){
		std::cerr<<std::endl<<"Error: RecoupleContent::read "<<std::endl;
		std::cerr<<" file ended before anything could be read...."<<std::endl;
		return;
	}

	int len;
	oo.read((char *)&len, sizeof(int));

	if(len==0){
		std::cerr<<std::endl<<"Error: RecoupleContent::read "<<std::endl;
		std::cerr<<" Nothing in file to read..."<<std::endl;
		return;
	}

	trains.resize(len);
	int sublen;
	Vector<int> curp;
	for(int i=0;i<len;++i){
		//read the length of the vector
		oo.read((char *)&sublen, sizeof(int));
		curp.resize(sublen);
		int j=0;
		while(j<sublen && !oo.eof()){	oo.read((char *)&curp[j], sizeof(int)); ++j;	}
		//cout<<i<<" "<<curp<<endl;
		if(oo.eof()){
			std::cerr<<std::endl<<"Error: RecoupleContent::read "<<std::endl;
			std::cerr<<" file ended before it should have..."<<std::endl;
			std::cerr<<" could not read all data..."<<std::endl;
			return;
		}
		trains[i]=curp;
	}
	while(!oo.eof()){
		oo.getline(liner, 100, '\n');
		if(std::string(liner)=="END RecoupleContent")	break;
	}
}

//write the trains into a file
void RecoupleContent::write(std::string outF)
{
	outFileName=outF;
	if(outFile && killFile){
		outFile->close();
		delete outFile;
	}
	killFile=true;
	outFile =new std::fstream(outF.c_str(), ios::binary | ios::out);
	write(*outFile);
}

void RecoupleContent::write(std::fstream &out)
{
	if(out.fail())
	{
		std::cerr<<std::endl<<"Error: RecoupleContent::read "<<std::endl;
		std::cerr<<" Cannot Write into the file........"<<std::endl;
		exit(1);
	}

	//The Header Chunk
	out<<std::endl<<"START RecoupleContent"<<std::endl;
	int csize=trains.size();

	//the number of entries
	out.write((char *)&csize, sizeof(int));

	for(int i=0;i<trains.size();++i){
		int cur=trains[i].size();
		out.write((char *)&cur, sizeof(int));
		//out<<cur<<" ";
		for(int wrct=0;wrct<trains[i].size();++wrct){
			cur=trains[i][wrct];
			out.write((char *)&cur, sizeof(int));
			//out<<cur<<" ";
		}
		//out<<endl;
	}
	out<<std::endl<<"END RecoupleContent"<<std::endl;
}

std::fstream &operator<<(std::fstream &out, RecoupleContent &oo)
{
	oo.write(out);
	return out;
}

std::fstream &operator>>(std::fstream &out, RecoupleContent &oo)
{
	oo.read(out);
	return out;
}


