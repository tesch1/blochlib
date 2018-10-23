


/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08-06-01
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*
	blochsolver.h--> The main integrator class for the bloch parameters setups
*/

#ifndef _solver_h_
#define _solver_h_


#include "bloch/bloch.h"
#include "timetrain/timetrain.h"
#include "driver/bs.h"
#include "driver/ckrk.h"
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"
#include "container/complex.h"
#include "utils/plotter.h"
#include "bloch/lyapunov.h"

BEGIN_BL_NAMESPACE


class SolverOps{
	public:
/* **Write Policies
SolverOps::Hold-->do not write
SolverOps::HoldUntilEnd-->when finished dump the data
SolverOps::Continous--> output data as we go
*/

enum writePol{Hold, HoldUntilEnd, Continous};	//hold all data till the end, or dump it continously

/* **Progress Bar options
SolverOps::On-->a little bar the travels across the screen
SolverOps::Off--> no progress output
*/

enum progress{On, Off};			//progress bar output

/* **data saving policies :: the 'data_' container below can get VERY large if
 there is alot of times steps and spins...resulting in lots of memmory and
 a slowing of the program to keep adding new data to the list

 to turn this off means 1) you are dumping data to file as you go
 						2) you only care about the final point in time and not the progression

 SolverOps::All--> collect all points (SolverOps::FIDonly and SolverOps::FinalPoint included)
 SolverOps::FIDonly--> colect just the fid data	(SolverOps::FinalPoint Included) ***FASTER****
 SolverOps::FinalPoint-> collect the final point						  ***FASTEST***
*/

enum collection{All=1,Magnetization=2, MagAndFID=3, FIDonly=4, FinalPoint=5};


/* Lyapunov flags
   When calculating the Variation Equations, it helps to have a reason to spend
   the large amount of time calculating them...this is why...to calculate the Lyaponuv exponents..

   SolverOps::HoldUntilEnd --> cacluate the lyapon but hold all the data until the end
   SolverOps::Continous--> calculate the lyapon and dump them out as we go
   NoCalc --> do not calculate them

*/

enum lyapunovcalc{LypHold, LypContinous, LypNoCalc};

};


template<class BlochEngine_t, class Solve_t=bs<BlochEngine_t, coord<> > >
class BlochSolver {

	private:
		BlochEngine_t *myBloch_;					//ptr to the function driver
		Vector<coord<> > *out_;					//pointer to the out from the ODESolver
		Vector<coord<> > InitialCon_;			//initial condition
		Solve_t driv_;	//ode solver default is the 'bs' driver...can be ckrk
		Vector< Vector< coord<> > > data_;		//total collected data

		Vector<complex> fid_;					//fid of the data
		Vector<coord<> > M_;		//(Mx, My,Mz) of detected spins...
		Vector<int> fidpos_;					//the positions in the data_ to grab the fid from
		double dt_;								//the time step for the FID collection (if any)
		double begint_;
		std::string detect_;					//detion spin ("1H", "13C", etc)

		Lyapunov<coord<> > myLyps_;				//Lypaunov data sets

		void AddBitsToFID(Vector<coord<> > *datin, int place)	//adds a complex pt to the fid data
		{
			if(fid_.size()<place){	fid_.resizeAndPreserve(place);	}
			if(gotfidpos_){
				complex tm(0.,0.);
				for(int j=0;j<fidpos_.size();j++){
					if(j<(*datin).size())	tm+=complex((*datin)(fidpos_(j)).x(), (*datin)(fidpos_(j)).y());
				}
				fid_(place)=(tm);
			}else{
				complex tm(0.,0.);
				if(LyapunovPol != SolverOps::LypNoCalc){
					for(int j=0;j<myBloch_->parameters()->size();j++){
						tm+=complex((*datin)(j).x(), (*datin)(j).y());
					}
				}else{
					for(int j=0;j<(*datin).size();j++){
						tm+=complex((*datin)(j).x(), (*datin)(j).y());
					}
				}
				fid_(place)=(tm);
			}
		}


		void AddBitsToM(Vector<coord<> > *datin, int place)	//adds a complex pt to the fid data
		{
			if(M_.size()<place){	M_.resizeAndPreserve(place);	}
			if(gotfidpos_){
				coord<> tm=0;
				for(int j=0;j<fidpos_.size();j++){
					if(j<(*datin).size())	tm+=(*datin)(fidpos_(j));
				}
				M_(place)=(tm);
			}else{
				if(LyapunovPol != SolverOps::LypNoCalc){
					M_(place)=sum((*datin)(Range(0,myBloch_->parameters()->size()-1)));
				}else{
					M_(place)=sum(*datin);
				}
			}
		}


		bool gottime_;
		bool gotinitial_;
		bool gotfidpos_;
		bool createdFID_;
		bool createdMag_;

		const void RunErr(){
			std::cout<<std::endl<<"Error : BlochSolver::solve(times)"<<std::endl;
			std::cout<<" No Bloch System Specified Cannot run simulation"<<std::endl;
		}

		const void FIDRunErr(){
			std::cout<<std::endl<<"Error : BlochSolver::setDetect(symbol)"<<std::endl;
			std::cout<<" No Bloch System Specified create fid tags"<<std::endl;
			std::cout<<" will proceed collecting ALL the spins"<<std::endl;
		}

		const void InitCondErr(){
			std::cout<<std::endl<<"Error : BlochSolver::solve(times)"<<std::endl;
			std::cout<<" No initial condition specified"<<std::endl;
			std::cout<<" Using the Bloch Engines default"<<std::endl;

		}

		const void NumErr(){
			std::cout<<std::endl<<std::endl<<"Error : BlochSolver::solve(times)"<<std::endl;
			std::cout<<" Number Overflow (is-nan detected).."<<std::endl;
			std::cout<<" Shrinking the time steps size and restarting section...."<<std::endl;

		}

		const void FIDsizeErr(){
			std::cout<<std::endl<<"Error : BlochSolver::setDetect(symbol)"<<std::endl;
			std::cout<<" No spin type in your parameters list..."<<std::endl;
			std::cout<<" ALL spins will be collected"<<std::endl;

		}

		const void CollectionErr(){
			std::cout<<std::endl<<"Error : BlochSolver::MakeFIDFromData()"<<std::endl;
			std::cout<<" Your data saving policy did not save the data"<<std::endl;
			std::cout<<" so i can create no FID for you..."<<std::endl;

		}

		const void WriteFIDerr(){
			std::cout<<std::endl<<"Error : BlochSolver::writeSpectrum()"<<std::endl;
			std::cout<<" Cannot Write FID..."<<std::endl;

		}

		const void WriteMerr(){
			std::cout<<std::endl<<"Error : BlochSolver::writeMag()"<<std::endl;
			std::cout<<" Cannot Write SolverOps::Magnetization..."<<std::endl;

		}

		std::string progBar;							//the progress bar
		void PrintProg(int onstep, int totsteps, double ontime, double finaltime){
			int x=(int)onstep*50/totsteps;
			int staroff=11;
			for(int i=staroff;i<x+staroff;i++)	progBar[i]='*';
			printf("time: %2.2f/%2.2f | ode steps:%d |", ontime,finaltime, driv_.ngood+driv_.nbad);
			printf(progBar.c_str());
			//flush(std::cout);
		}



	public:

		std::string rawdatafilename_;				//file to dump all info from 'data_'
		std::string LypFilename_;					//file name to dump lypaunov eqs to

		SolverOps::writePol writePol;					//out file writing policy
		SolverOps::progress progressBar;					//progress bar printing
		SolverOps::collection collectionPol;
		SolverOps::lyapunovcalc LyapunovPol; 		//lypunoc calulation policy

		std::ios::openmode rawoutios;	//raw data writing policy...either 'std::ios::out' or 'std::ios::append'

		BlochSolver():
			myBloch_(NULL), out_(NULL), dt_(0), detect_(""),
			gottime_(false), gotinitial_(false),gotfidpos_(false),
			createdFID_(false),createdMag_(false), rawdatafilename_("BlochSolver.out"),LypFilename_("Lyapunov.out"),
			writePol(SolverOps::Hold), progressBar(SolverOps::On), collectionPol(SolverOps::MagAndFID ),LyapunovPol(SolverOps::LypNoCalc),
			rawoutios(std::ios::out)

		{
			//data_.RecordSep="\n"; 	//sets the print policy for the out put to be newline ended
			progBar="progress: [..................................................]\r";
		}

		BlochSolver(
			BlochEngine_t &inBloch,
			const Vector<coord<> > &init,
			std::string rawdataname_="BlochSolver.out",
			std::string lypfile_="Lyapunov.out",
			SolverOps::collection cl=SolverOps::MagAndFID ,
			SolverOps::writePol wr=SolverOps::Hold,
			SolverOps::progress ts=SolverOps::On,
			SolverOps::lyapunovcalc lyp=SolverOps::LypNoCalc):

				myBloch_(&inBloch),
				InitialCon_(init),
				driv_(init,0,.1,.1, *myBloch_),
				dt_(0), detect_(""),
				gottime_(false), gotinitial_(true),gotfidpos_(false),
				createdFID_(false),createdMag_(false),rawdatafilename_(rawdataname_),LypFilename_(lypfile_),
				writePol(wr), progressBar(ts), collectionPol(cl),LyapunovPol(SolverOps::LypNoCalc),
				rawoutios(std::ios::out)

		{
			//data_.RecordSep="\n"; 	//sets the print policy for the out put to be newline ended
			out_=(&InitialCon_);
			driv_.setInitialCondition(InitialCon_);
			progBar="progress: [..................................................]\r";
		}

		BlochSolver(BlochSolver &copy):
				myBloch_(copy.myBloch),
				out_(copy.out_),
				InitialCon_(copy.InitialCon_),
				driv_(copy.InitialCon_,0,.1,.1, myBloch_),
				dt_(copy.dt_),
				detect_(copy.detect_),
				gottime_(copy.gottime_), gotinitial_(copy.gotinitial_),
				gotfidpos_(copy.gotfidpos_),
				rawdatafilename_(copy.rawdatafilename_),
				LypFilename_(copy.LypFilename_),
				writePol(copy.writePol),
				progressBar(copy.progressBar),
				collectionPol(copy.collectionPol),
				LyapunovPol(copy.calcLyapunov),
				createdFID_(copy.createdFID_),
				createdMag_(copy.createdFID_),
				rawoutios(cp.rawoutios)

		{
			//data_.RecordSep="\n"; 	//sets the print policy for the out put to be newline ended
			progBar="progress: [..................................................]\r";
		}

		~BlochSolver(){	myBloch_=NULL; out_=NULL;	}

//set the raw data output filename
		void setRawOut(std::string inname, std::ios::openmode iosf=std::ios::out){	rawdatafilename_=inname; rawoutios=iosf;	}

//setting policies after initialization
		void setWritePolicy(SolverOps::writePol in)			{	writePol=in;			}
		void setProgressBar(SolverOps::progress in)			{	progressBar=in;			}
		void setCollectionPolicy(SolverOps::collection in)	{	collectionPol=in;		}


/*****setting the bits for Lyspunov calculations */
	//setting policies
		void setLyapunovPolicy(SolverOps::lyapunovcalc in)
		{
			if(myBloch_->VariationalPolicy()==BlochOps::Variational && in!=SolverOps::LypNoCalc){
				LyapunovPol=in;
			}else if(in != SolverOps::LypNoCalc){
				std::cout<<std::endl<<"Warning: BlochSovler::setLyapunovPolicy()"<<std::endl;
				std::cout<<" setting the 'Bloch' to calculate the BlochOps::Variational Equations"<<std::endl;
				std::cout<<" Otherwise, i could not calculate the Lyapunov's"<<std::endl;
				myBloch_->calcVariational();
				LyapunovPol=in;
			}else{
				LyapunovPol=in;
			}
		}

	//set a lyapunov from an external decalred class
		void setLyapunov(Lyapunov<coord<> > &in){	myLyps_=in;	}
	//set the 'calculation step' for the Lyapunov class
		void setLyapunovStep(int in){	setCalcStep(in);	}
	//set output file name
		void setLypDataFile(const std::string &in){	LypFilename_=in;	}
/****************/

//set the initial condition...
		void setInitialCondition(const Vector<coord<> > &InitialCon) 
		{
			InitialCon_=InitialCon;
			driv_.setInitialCondition(InitialCon_);
			out_=driv_.solvedData();
		}

//set the initial condition for the variational euqations...
		void setVariationalInitCond(const Vector<coord<> > &InitialCon) 
		{
			if(InitialCon_.size()==InitialCon.size()/3){
				int pres=myBloch_->parameters()->size();
				InitialCon_.resizeAndPreserve(InitialCon_.size()*4);
				InitialCon_(Range(pres, InitialCon_.size()-1))=InitialCon;
			}else if(InitialCon_.size()==InitialCon.size()){
				InitialCon_=InitialCon;
			}else{
				std::string _mess = 
				" Size mismatch of input"+
				"\n -input should be 4 * the size of the system"+
				"\n if setting both the Magnitiziation AND BlochOps::Variational"+
				"\n -input should be 3 * the size of the system"+
				"\n if setting only BlochOps::Variational";
				BLEXCEPTION(_mess)
			}
			driv_.setInitialCondition(InitialCon_);
			out_=driv_.solvedData();
		}

//set the Bloch data set...
		void setBloch(BlochEngine_t &in){	myBloch_=&in;	}

//this fixes the TimeTrain to a smaller 'dt' if it is able to inorder
// to compensate for any 'isnan' found in the solver

		template<class TimeEng>
		TimeTrain<TimeEng> FixTimer(TimeTrain<TimeEng> &in)
		{

			return TimeTrain<TimeEng>(TimeEng(
									in.beginTime(),
									in.endTime(),
									in.size()*2,
									in.step(1)));
		}


//Main driver
		template<class TimerEng>
		bool solve(TimeTrain<TimerEng> &therun) 
		{
			progBar="progress: [..................................................]\r";
			if(!myBloch_){
				RunErr();
				return false;
			}

			std::ofstream oo;
			std::ofstream lyoo;
			if(writePol==SolverOps::Continous ||
				writePol==SolverOps::HoldUntilEnd)
			{
				oo.open(rawdatafilename_.c_str(), rawoutios);
			}
			if(LyapunovPol==SolverOps::LypContinous ||
			   LyapunovPol==SolverOps::LypHold)
				{	lyoo.open(LypFilename_.c_str());	}

			gottime_=true;
			int ct=1;
			out_=driv_.get_out();			//set ouput ptr...

			if(LyapunovPol != SolverOps::LypNoCalc)
			{
				myLyps_.setData(out_);			//set ptr insoed the Lyp class
				myLyps_.setSize(myBloch_->size());
				myLyps_.reset();
			}

			typename TimeTrain<TimerEng>::iterator Itrun(therun);

			int numsteps=therun.size();
			double ontime=Itrun.beginStepTime();
			begint_=therun.beginTime();
			double finaltime=therun.endTime();

			if(collectionPol==SolverOps::All)
			{
				data_.resize(numsteps);
				data_(0)= (*out_);
			}

			if(collectionPol==SolverOps::All ||
				collectionPol==SolverOps::FIDonly ||
				collectionPol==SolverOps::MagAndFID)
			{
				fid_.resize(numsteps,0);
				AddBitsToFID(out_, 0);
			}

			if(collectionPol==SolverOps::All ||
				collectionPol==SolverOps::Magnetization ||
				collectionPol==SolverOps::MagAndFID)
			{
				M_.resize(numsteps,0);
				AddBitsToM(out_, 0);
			}
			//dump out the first point...
			if(writePol==SolverOps::Continous)
				oo<<(*out_)(Range(0, myBloch_->parameters()->size()-1))
				  <<std::endl;

			dt_=Itrun.endStepTime()-Itrun.beginStepTime();
			//therun.reset();
			while(Itrun){
				if(progressBar==SolverOps::On)
				   PrintProg(ct, numsteps, ontime, finaltime);

				//cout<<"t: "<<Itrun.beginStepTime()<<(*out_)<<endl;
				driv_.set_time(Itrun.beginStepTime(), Itrun.endStepTime());
				driv_.odeint();

			#ifdef HAVE_ISNAN
				if(hasnan(*out_))
				{
					std::cout<<std::endl
					<<std::endl<<"***On Time:::"
					<<Itrun.beginStepTime()<<std::endl;

					NumErr();
					therun=TimeTrain<TimerEng>(
						therun.beginTime(),
						therun.endTime(),therun.size(),
						therun.substep(0)*2);

					numsteps=therun.size();
					typename TimeTrain<TimerEng>::iterator tmit(therun);
					Itrun=tmit;
					Itrun.reset();
					ct=1;
					progBar="progress: [..................................................]\r";
					driv_.set_y(InitialCon_);
					fid_.resize(numsteps,0);
					if(collectionPol==SolverOps::All){
						data_.resize(numsteps);
						data_(0)= (*out_);
					}else if(collectionPol==SolverOps::FIDonly){
						AddBitsToFID(out_, 0);
					}
				}
			#endif

				ontime=Itrun.beginStepTime();

				if(collectionPol==SolverOps::All)
					{	data_(ct)= (*out_);	}

				if(LyapunovPol==SolverOps::LypContinous
				|| LyapunovPol==SolverOps::LypHold)
					{	myLyps_.calcLyapunov(ontime, Itrun.dt());}

				if(LyapunovPol==SolverOps::LypContinous)
					{	lyoo<<myLyps_; }

				if(collectionPol==SolverOps::All ||
					collectionPol==SolverOps::FIDonly ||
					collectionPol==SolverOps::MagAndFID)
					{	AddBitsToFID(out_, ct); }

				if(collectionPol==SolverOps::All ||
					collectionPol==SolverOps::Magnetization ||
					collectionPol==SolverOps::MagAndFID )
					{	AddBitsToM(out_, ct);	}

				if(writePol==SolverOps::Continous)
				{	oo<<(*out_)(Range(0, myBloch_->parameters()->size()-1))<<std::endl;}

				//std::cout<<std::endl<<sum(*out_)<<std::endl;
				ct++;
				++Itrun;

			}
			if(writePol==SolverOps::HoldUntilEnd &&
				collectionPol==SolverOps::All)
				{	writeData(oo);	}

			if(collectionPol==SolverOps::All ||
				collectionPol==SolverOps::FIDonly ||
				collectionPol==SolverOps::MagAndFID)
				{ createdFID_=true;	}

			if(collectionPol==SolverOps::All ||
				collectionPol==SolverOps::Magnetization ||
				collectionPol==SolverOps::MagAndFID)
				{ createdMag_=true;	}

			if(LyapunovPol==SolverOps::LypHold)
				{lyoo<<myLyps_<<std::endl;	}
			return true;
		}


	//set the fidpos_...this is a list of the specific spin types posistions in the data lists
		void setDetect(std::string in)  
		{
			detect_=in;
			if(myBloch_){

				typename BlochEngine_t::iterator
					myIt(*(myBloch_->parameters()));

				int i=0;
				while(myIt)
				{
					if(myIt.symbol()==in){
						fidpos_.push_back(i);
					}
					++i;
					++myIt;
				}
				if(fidpos_.size()==0){
					FIDsizeErr();
					gotfidpos_=false;
				}else{
					gotfidpos_=true;
				}
			}else{
				FIDRunErr();
			}
		}

		void setDetect(Vector<int> &in) 
		{
			fidpos_=in;
			if(fidpos_.size()==0
				|| max(fidpos_>=myBloch_.size()))
			{
				FIDsizeErr();
				gotfidpos_=false;
			}else{
				gotfidpos_=true;
			}
		}
		
		void setDetect(Range &in)
		{
			fidpos_=in;
			if(fidpos_.size()==0
				|| max(fidpos_>=myBloch_.size()))
			{
				FIDsizeErr();
				gotfidpos_=false;
			}else{
				gotfidpos_=true;
			}
		}

		bool createFIDfromData() 
		{

			if(collectionPol==SolverOps::FIDonly ||
				collectionPol==SolverOps::All ) return true;

			if(collectionPol==SolverOps::FinalPoint)
				{ CollectionErr(); return false;	}

			if(!createdFID_){
				fid_.resize(data_.size(),0);
				if(gotfidpos_){
					for(int i=0;i<data_.size();i++){
						for(int j=0;j<fidpos_.size();j++){
							if(j<data_.size()){
								fid_(i)+=complex(data_(i)(j).x(), data_(i)(j).y());
							}
						}
					}
				}else{
					for(int i=0;i<data_.size();i++){
						for(int j=0;j<data_(i).size();j++){
							fid_(i)+=complex(data_(i)(j).x(), data_(i)(j).y());
						}
					}
				}
			}
			createdFID_=true;
			return true;
		}

		bool createMfromData() 
		{

			if(collectionPol==SolverOps::Magnetization ||
				collectionPol==SolverOps::All ||
				collectionPol==SolverOps::MagAndFID)
					return true;

			if(collectionPol==SolverOps::FinalPoint)
				{ CollectionErr(); return false;	}

			if(!createdFID_){
				M_.resize(data_.size(),0);
				if(gotfidpos_){
					for(int i=0;i<data_.size();i++){
						for(int j=0;j<fidpos_.size();j++){
							if(j<data_.size())	M_(i)+=data_(i)(j);
						}
					}
				}else{
					for(int i=0;i<data_.size();i++){
						for(int j=0;j<data_(i).size();j++){
							M_(i)+=data_(i)(j);
						}
					}
				}
			}
			createdMag_=true;
			return true;
		}

		Vector<complex> FID()	{	return fid_;	}
		Vector<coord<> > M()	{	return M_;	}

//the last point integrated...
		Vector<coord<> > lastPoint(){	return *out_;	}

		bool writeSpectrum(std::string fname) 
		{
			if(createFIDfromData()){
				if(dt_!=0){
					std::ofstream oo(fname.c_str());
					Vector<complex> fft(fid_.copy());
					if(fid_.size()%2==0){	fft=FFT(fft);	}

					double dw=1./dt_/fid_.size();
					double min_w=-1./dt_/2.;
					double t=begint_, w=min_w;
					for(int i=0;i<fid_.size();i++){
						oo<<t<<" "<<w<<" "<<Re(fft(i))<<" "
						<<Im(fft(i))<<" "<<Re(fid_(i))<<" "
						<<Im(fid_(i))<<" "<<norm(fft(i))<<std::endl;

						t+=dt_;
						w+=dw;
					}
				}else{
					std::ofstream oo(fname.c_str());
					for(int i=0;i<fid_.size();i++){
						oo<<Re(fid_(i))<<" "<<Im(fid_(i))<<std::endl;
					}
				}
			}else{
				WriteFIDerr();
				return false;
			}

			return true;
		}

		bool writeSpectrum(std::ofstream &oo) 
		{
			if(createFIDfromData()){
				if(dt_!=0){
					Vector<complex> fft=fid_;
					//if(fid_.size()%2==0){	fft=FFT(fid_);	}
					double dw=1./dt_/fid_.size();
					double min_w=-1./dt_/2.;
					double t=begint_, w=min_w;
					for(int i=0;i<fid_.size();i++){
						oo<<t<<" "<<w<<" "<<Re(fft(i))<<" "
						<<Im(fft(i))<<" "<<Re(fid_(i))<<" "
						<<Im(fid_(i))<<" "<<norm(fft(i))<<std::endl;

						t+=dt_;
						w+=dw;
					}
				}else{
					for(int i=0;i<fid_.size();i++){
						oo<<Re(fid_(i))<<" "<<Im(fid_(i))<<std::endl;
					}
				}
			}else{
				WriteFIDerr();
				return false;
			}

			return true;
		}

		bool writeMag(std::string fname) 
		{
			if(createMfromData()){
				double t=begint_;
				if(dt_!=0){
					std::ofstream oo(fname.c_str());
					for(int i=0;i<M_.size();i++){
						oo<<t<<" "<<M_[i]<<std::endl;
						t+=dt_;
					}
				}else{
					std::ofstream oo(fname.c_str());
					for(int i=0;i<M_.size();i++){
						oo<<M_[i]<<std::endl;
					}
				}
			}else{
				WriteMerr();
				return false;
			}

			return true;
		}

		bool writeMag(std::ofstream &oo) 
		{
			if(createMfromData()){
				double t=begint_;
				if(dt_!=0){
					for(int i=0;i<M_.size();i++){
						oo<<t<<" "<<M_[i]<<std::endl;
						t+=dt_;
					}
				}else{
					for(int i=0;i<M_.size();i++){
						oo<<M_[i]<<std::endl;
					}
				}
			}else{
				WriteMerr();
				return false;
			}

			return true;
		}


		void writeData(std::ostream &oo) 
		{
			for(int i=0;i<data_.size();++i){
				oo<<data_(i)<<std::endl;
			}
		}



};



END_BL_NAMESPACE


#endif
