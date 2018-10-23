

#include "blochlib.h"
#include "propreduce.h"
using namespace BlochLib;
using namespace std;

timer stopwatch;

class mmqmas{
	private:
		Vector<int> FromEle, ToEle;	//the start and final matrix elements
		double RealTime; 	//Realtime for the DT between FID points

/** Central Transistion FID calculators for single powder point**/
		void calcCentral(int curp);
	//calculates central Trans asuming diagonal hamiltonians
		void calcCentralDiag(int curp);
	//calculates central Trans asuming NOT diag matrix
		void calcCentralNonDiag(int curp);

/** MMQMAS Calculators for single powder point **/
		void calcMMQMAS(int curp);

	public:
		int npts; 			//number of fid points
		double maxtstep; 	//max dt for an integrations step
		double wr, rotor; 	//rotor speed and angle
		double sw; 			//sweep width for central trans calc

		SolidSys mysys; //the spin system

		matrix roeq;	//initial density matrix
		powder PowderAve;	//powder ave container
		powder::iterator powit; //the powder iterator
		PropReduce myReduce;	//the propgator reduction object

		Vector<complex> fid;	//the FID

		mmqmas():
			npts(512), maxtstep(1e-6), wr(10000), rotor(acos(1/sqrt(3.0)))
		{}

		void setElements(const Vector<int> &From,const  Vector<int> &To);

		void calcCentral();
		void calcMMQMAS();

		void writeFID(std::string fidout);
};

void mmqmas::setElements(const Vector<int> &From,const  Vector<int> &To)
{

	if(From.size() !=2 || To.size() !=2){
		std::cerr<<std::endl<<"'FromElement' and/or 'ToElement' length should be 2"<<std::endl;
		std::cerr<<" each set of 2 is a matrix element"<<std::endl;
		std::cerr<<" like 0,2..."<<std::endl;
		exit(1);
	}

	if(From.size() != To.size() ){
		std::cerr<<std::endl<<"'FromElement' should the same size 'ToElement' length should be same"<<std::endl;
		std::cerr<<" as an examples:::"<<std::endl;
		std::cerr<<" FromElement 1,2"<<std::endl;
		std::cerr<<" ToElement 0,3"<<std::endl;
		exit(1);
	}
	FromEle=From;
	ToEle=To;

//make the detect list
	if(mysys.size()>=2){
		rmatrix holdD(mysys[0].HS(), mysys[0].HS(), 0.0);
		holdD(FromEle[0], FromEle[1])=1.0;
		rmatrix holdR(mysys[0].HS(), mysys[0].HS(), 0.0);
		holdR(ToEle[0], ToEle[1])=1.0;
		rmatrix bigD1=cross(rimatrix(mysys[0].HS(),mysys[0].HS()), holdD);
		rmatrix bigD2=cross(holdD,rimatrix(mysys[1].HS(),mysys[1].HS()));
		rmatrix bigR1=cross(rimatrix(mysys[0].HS(),mysys[0].HS()), holdR);
		rmatrix bigR2=cross(holdR,rimatrix(mysys[1].HS(),mysys[1].HS()));
		Vector<int> newD;
		Vector<int> newR;
		for(int i=0;i<bigD1.rows();++i){
			for(int j=0;j<bigD1.cols();++j){
				if(bigD1(i,j)!=0.0){
					newD.push_back(i);
					newD.push_back(j);
				}
				if(bigD2(i,j)!=0.0){
					newD.push_back(i);
					newD.push_back(j);
				}
				if(bigR1(i,j)!=0.0){
					newR.push_back(i);
					newR.push_back(j);
				}
				if(bigR2(i,j)!=0.0){
					newR.push_back(i);
					newR.push_back(j);
				}
			}
		}
		roeq=bigD1+bigD2;
		FromEle=newD;
		ToEle=newR;
	}else{
		rmatrix holdD(mysys[0].HS(), mysys[0].HS(), 0.0);
		holdD(FromEle[0], FromEle[1])=1.0;
		roeq=holdD;
	}
}

/************** CENTRAL TRANSITION CALCUATION ************/
/** Master MPI function for Caclualtion of Central Trans Spectra **/
void mmqmas::calcCentral()
{
	fid.resize(npts);
	fid.fill(0.0);
	maxtstep=(1.0/wr)/std::floor((1.0/wr)/maxtstep);
	if(MPIworld.master()) std::cout<<"real Delta T="<<maxtstep<<std::endl;
	RealTime=1.0/wr;
	if(MPIworld.master()) cout<<"---(powder pt)/(total)---"<<endl;
	int pos=0, done=-1;
	powit=powder::iterator(PowderAve);
	if(MPIworld.master()){
		pos=powit.curpos();
		for(int i=1;i<MPIworld.size();++i){
			pos=powit.curpos();
			MPIworld.put(pos, i); ++powit;
			if(!powit) break;
		}
		int get, dummy;
		while(powit){
			if(MPIworld.serial()){
				calcCentral(powit.curpos());
			}else{
				get=MPIworld.getAny(dummy);
				pos=powit.curpos();
				MPIworld.put(pos, get);
			}
			cout<<"   "<<itost_form("%10d",powit.curpos())<<"\\"<<powit.size()<<"\r";
			cout.flush();
			++powit;
		}
		for(int i=1;i<MPIworld.size();++i)	MPIworld.put(done, i);

	}else if(MPIworld.parallel()){
		while(1)
		{
			MPIworld.get(pos, 0);
			if(pos==done) break;
			//powit.begin(pos);
			calcCentral(pos);
			MPIworld.put(pos, 0);
		}
	}
	MPIworld.reduce(fid, Reduce::Add);
}

void mmqmas::calcCentral(int curp)
{
	if(mysys.dip.size()==0 && mysys.jcop.size()==0) calcCentralDiag(curp);
	else calcCentralNonDiag(curp);
}

//calculates central Trans asuming diagonal hamiltonians
void mmqmas::calcCentralDiag(int curp)
{

//Diagonal Matrix propogation Q and CSA are diagonal
	static dmatrix Ut1;
	//static dmatrix Utemp;
	static matrix curro;
	Ut1=mysys.Fe();
	mysys.setPowderAngles(PowderAve.theta(curp), PowderAve.phi(curp), PowderAve.gamma(curp));
	double tstop=0.0, tbegin=0.0;

	//For CENTRAL TRans Only
	tstop=1.0/wr;
	mysys.setRotorAngles(1.0,rotor);
	for(double t=tbegin;t<tstop;t+=maxtstep){
		Ut1=Mexp(dmatrix(mysys.Hamiltonian(t, t+maxtstep, wr)), -complexi*maxtstep*PI2)*Ut1;
		//cout<<mysys.theRotations.alpha*RAD2DEG<<" "<<endl;
	}

	//Utemp=Ut1;
	curro=roeq;
	for(int k=0;k<ToEle.size();k+=2)
		fid[0]+=curro(ToEle[k], ToEle[k+1])*PowderAve.weight(curp);

	for(int j=1;j<npts;++j){
		curro.prop(Ut1);
		for(int k=0;k<ToEle.size();k+=2)
			fid[j]+=curro(ToEle[k], ToEle[k+1])*PowderAve.weight(curp);
		//Ut1=Utemp*Ut1;
	}
}

//calculates central Trans asuming NOT diag matrix
//Full Matrix propogation becuase Dip and J are off diagonal
void mmqmas::calcCentralNonDiag(int curp)
{
	static matrix Ut1;
	//static matrix Utemp;
	static matrix curro;
	Ut1=mysys.Fe();
	mysys.setPowderAngles(PowderAve.theta(curp), PowderAve.phi(curp), PowderAve.gamma(curp));
	double tstop=0.0, tbegin=0.0;

	//For CENTRAL TRans Only
	double now=stopwatch();
	tstop=1.0/wr;
	mysys.setRotorAngles(1.0,rotor);
	for(double t=tbegin;t<tstop;t+=maxtstep){
		Ut1=Mexp(mysys.Hamiltonian(t, t+maxtstep, wr), -complexi*maxtstep*PI2)*Ut1;
	}

	//Utemp=Ut1;
	curro=roeq;
	for(int k=0;k<ToEle.size();k+=2)
		fid[0]+=roeq(ToEle[k], ToEle[k+1])*PowderAve.weight(curp);

	for(int j=1;j<npts;++j){
		curro.prop(Ut1);
		for(int k=0;k<ToEle.size();k+=2)
			fid[j]+=curro(ToEle[k], ToEle[k+1])*PowderAve.weight(curp);

		//Ut1=Utemp*Ut1;
	}
}

/************** MMQMAS TRANSITION CALCUATION ************/
/** Master MPI function for Caclualtion of MMQMAS Spectra **/
void mmqmas::calcMMQMAS()
{
	fid.resize(npts);
	fid.fill(0.0);
	int pos=0, done=-1;
	powit=powder::iterator(PowderAve);
	double t1=1.0/wr;
	maxtstep=(t1/double(myReduce.base))/std::ceil((t1/double(myReduce.base))/maxtstep);
	if(maxtstep>1.0/wr/myReduce.base ){
		std::cerr<<std::endl<<"Your 'maxtstep' is too big"<<std::endl;
		std::cerr<<" for the spinning speed and Base you have choosen"<<std::endl;
		std::cerr<<" Please remedy this situation."<<std::endl;
		exit(1);
	}
	if(MPIworld.master()) std::cout<<"real Delta T="<<maxtstep<<std::endl;
	if(MPIworld.master()) cout<<"---(powder pt)/(total)---"<<endl;
	RealTime=t1+double(myReduce.factor)/double(myReduce.base)*t1;


	if(MPIworld.master()){
		pos=powit.curpos();
		for(int i=1;i<MPIworld.size();++i){
			pos=powit.curpos();
			MPIworld.put(pos, i); ++powit;
			if(!powit) break;
		}
		int get, dummy;
		while(powit){
			if(MPIworld.serial()){
				calcMMQMAS(powit.curpos());
			}else{
				get=MPIworld.getAny(dummy);
				pos=powit.curpos();
				MPIworld.put(pos, get);
			}
			cout<<"   "<<itost_form("%10d",powit.curpos())<<"\\"<<powit.size()<<"\r";
			cout.flush();
			++powit;
		}
		for(int i=1;i<MPIworld.size();++i)	MPIworld.put(done, i);

	}else if(MPIworld.parallel()){
		while(1)
		{
			MPIworld.get(pos, 0);
			if(pos==done) break;
			powit.begin(pos);
			calcMMQMAS(pos);
			MPIworld.put(pos, 0);
		}
	}
	MPIworld.reduce(fid, Reduce::Add);
}

void mmqmas::calcMMQMAS(int pos)
{
	powit.begin(pos);
	static  matrix Ut1=mysys.Fe();
	static 	matrix Ut2=mysys.Fe();
	static 	matrix curro=roeq;
	static 	matrix Utemp;
	static 	Vector<matrix> Foward(myReduce.base-1);
	static 	Vector<matrix> Backs(myReduce.maxBackReduce());
	static Vector<matrix> Usubs(myReduce.base);
	static Vector<matrix> Ut2s(myReduce.base);

//so as to fit within our sub rotor chunk
	double t1=1.0/wr;
	double t2step=double(myReduce.factor)/double(myReduce.base)*t1;

	Usubs.fill(mysys.Fe());
	Ut2s.fill(mysys.Fe());
	Ut2=mysys.Fe();

	//set the powder angles of the spin system
	mysys.setPowderAngles(powit.theta(), powit.phi(), powit.gamma());

	//being calculating the sub-chunks of propogators

	double tstop=0.0, tbegin=0.0,t;
	for(int i=0;i<myReduce.base; ++i){
		tbegin=double(i)/double(myReduce.base)*t1;
		tstop=double(i+1)/double(myReduce.base)*t1;
		int j=0;
		for(t=tbegin+maxtstep/2.0;t<tstop;t+=maxtstep)
		{
			mysys.setRotorAngles(PI2*(t)*wr,rotor);
			if(j==0){
				Usubs[i]=Mexp(mysys.H(), complexi*maxtstep*PI2); ++j;
			}else{
				Usubs[i]=Mexp(mysys.H(), complexi*maxtstep*PI2)*Usubs[i];
			}
		}
		if(i==1 && myReduce.base>2) Foward[i-1]=Usubs[i]*Usubs[i-1];
		if(i>1)		Foward[i-1]=Usubs[i]*Foward[i-2];
	}
	int ct=0, ct2=0;

	//calculate the 'back propogators'
	if(Usubs.size()>=2 && Backs.size()>0){
		Backs[0]=Usubs[myReduce.base-1]*Usubs[myReduce.base-2];
	}
	int jB=myReduce.base-3;
	for(int i=1;i<Backs.size();++i){
		Backs[i]=Backs[i-1]*Usubs[jB];
		--jB;
	}
	myReduce.generateProps(Usubs, Foward, Backs, Ut2s);

	Utemp=Foward[myReduce.base-2];
	for(int k=0;k<FromEle.size();k+=2)
		fid[0]+= roeq(FromEle[k], FromEle[k+1])*powit.weight();

	for(int j=1;j<npts;++j)
	{
		//rotor evolution
		curro=adjprop(Utemp, roeq);

		//Coherence transfer transfer (perfect)
		for(int k=0;k<FromEle.size();k+=2)
			curro(ToEle[k], ToEle[k+1])=curro(FromEle[k], FromEle[k+1]);

		//advanace propogators in T2
		Ut2=Ut2s[(j-1)%myReduce.base]*Ut2;

		//the t2 evolution and collection
		curro=adjprop(Ut2, curro);

		//we can only detect the
		for(int k=0;k<ToEle.size();k+=2)
			fid[j]+=curro(ToEle[k], ToEle[k+1])*powit.weight();

		//now=stopwatch();
		if(j!=npts-1) Utemp=Foward[myReduce.base-2]*Utemp;

	}
}

void mmqmas::writeFID(std::string fidout)
{
	if(MPIworld.master()) plotterFID(fid, fidout, RealTime);
}

int main(int argc,char* argv[]){

	MPIworld.start(argc, argv);
	double startT=stopwatch(), endT;
	std::string fname;
	if(MPIworld.master()) query_parameter(argc, argv, 1, "Enter in filie name to parse: ", fname);
	MPIworld.scatter(fname);

	Parameters pset(fname);
	pset.addSection("spins");

	mmqmas myRunner;

	myRunner.mysys=SolidSys(pset.section("spins"));
	pset.addSection("params");
	myRunner.mysys.setBfield(pset.getParamD("Bfield", "params")*1e6);

	myRunner.npts=pset.getParamI("npts", "params");
	myRunner.sw=pset.getParamD("sw", "params", false, 50000.0);
	myRunner.wr=pset.getParamD("wr", "params");
	myRunner.rotor=pset.getParamD("rotor", "params")*PI/180.0;

	myRunner.PowderAve=powder(pset.getParamS("avetype", "params"),
	                          pset.getParamI("thetasteps", "params", false),
	                          pset.getParamI("phisteps", "params", false),
	                          pset.getParamI("gammasteps", "params",false));

	std::string initS=pset.getParamS("initro", "params", false, "Iz");
	std::string detS=pset.getParamS("detect", "params", false, "Ip");

	HamiltonianGen hamgen;
	myRunner.roeq=hamgen.Hamiltonian(myRunner.mysys, initS);
	myRunner.roeq/=sqrt(trace(myRunner.roeq,adjoint(myRunner.roeq)));

	if(MPIworld.master())	cout<<myRunner.mysys<<endl;

	myRunner.maxtstep=pset.getParamD("maxtstep", "params");
	myRunner.setElements(pset.getParamVectorI("FromElement", "params"),
				         pset.getParamVectorI("ToElement", "params"));


	if(MPIworld.master()) cout<<"INITIALIZATION TIME: "<<stopwatch()-startT<<endl;

	std::string expType=pset.getParamS("expType", "params");

	if(expType=="centralTrans"){
		myRunner.calcCentral();
	}else if(expType=="mqmas" || expType=="mmqmas"){
		int Base=pset.getParamI("baseSym", "params");
		int factorSym=pset.getParamI("extent", "params");
		ofstream oo("Reduce.out");
		myRunner.myReduce=PropReduce(Base, factorSym, &oo);
		myRunner.myReduce.reduce();
		myRunner.calcMMQMAS();
	}else{
		if(MPIworld.master()) cout<<"Error: "<<expType<<" not yet implimented "<<endl;
		MPIworld.end(); return 0;
	}
	myRunner.writeFID(pset.getParamS("dataout", "params", false, "fid"));

	MPIworld.end();
}

