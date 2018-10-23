

#include "propogation.h"
#include "sequenceparse.h"
#include "paramset.h"
#include "sequencerun.h"

#include <signal.h>
#include <stdio.h>

using namespace BlochLib;
using namespace std;

//this catchs any sigkills s to terminbate MPI
int Solid_gotsig=0;
extern "C"{
void sigdie(int);
void sigdie(int)
{
/*	signal(SIGHUP, SIG_DFL);
	signal(SIGQUIT ,SIG_DFL);
	signal(SIGKILL  , SIG_DFL);
	signal(SIGINT  , SIG_DFL);
	signal(SIGSEGV  , SIG_DFL);
	signal(SIGABRT , SIG_DFL);
*/
	if(MPIworld.master()){
		printf("\n\n ***Caught Termination signal\n");
		printf(" ***exiting Solid...\n");
	}
	if(MPIworld.size()>1)
	{
		if(MPIworld.master()){
			printf(" ***NOTE: Aborting MPI programs can leave zombies...\n");
			printf(" ***NOTE: Make sure to KILL THEM ALL...\n");
			printf(" ***NOTE: look for them via a 'ps -ex' (on linux & cygwin) command...\n");
		}
		MPIworld.abort();
	}
	Solid_gotsig=1;
//	int err=0;
//	MPIworld.abort();
//	exit(1);
	signal(SIGINT, SIG_DFL);
	signal(SIGQUIT, SIG_DFL);
	raise(SIGINT);
//	if(!MPIworld.end()){
//		printf(" ERROR: could not kill MPI...you may need to  \n");
//		printf(" kill them yourself via 'killall solid' command  ");
//	}
}
}

//gets me the output from a system command
std::string getComand(std::string cmd)
{
    /* execute the command */
std::FILE* cmdp =  popen(cmd.c_str(), "r");
    if (!cmdp) {
std::perror("popen");
        return("");
    }
std::string out="";
    char result[256];
    /* now read the "grep" command outputs */
    while (fgets(result, sizeof(result), cmdp))
        out+=std::string(result);   // add the string to the end
    pclose(cmdp);                // close the stream
    return (out);
}
    
int main(int argc,char* argv[])
{

	/*int pid=0;
	if ((pid = fork()) < 0) {
		perror("fork");
		exit(1);
	}
	std::cout<<pid<<std::endl;
	if(pid==0){
*/
		//set the singal handler functions
		signal(SIGHUP, sigdie);
		signal(SIGQUIT ,sigdie);
		signal(SIGKILL  , sigdie);
		signal(SIGINT  , sigdie);
		signal(SIGSEGV  , sigdie);
		signal(SIGABRT , sigdie);

//		for(;;);
//	}else{
	//while(Solid_gotsig==0){
		MPIworld.start(argc, argv);


		if(MPIworld.master())
		{
			std::cout<<std::endl<<" Solid 2.0"<<std::endl;
			std::cout<<"..a program to perform solid state NMR simulations "<<std::endl;
			std::cout<<" for documentation and licence see http://waugh.cchem.berkeley.edu/solid/ "<<std::endl;
			std::cout<<"  Last compiled on "<<__DATE__<<" "<<__TIME__<<std::endl;
			std::cout<<" Current Working directory:"<<getComand("pwd")<<std::endl;
            
		}

		std::string fname;
		if(MPIworld.master())
			query_parameter(argc,argv,1, "Enter the input file to parse: ", fname);

		MPIworld.scatter(fname);

		int Error=0, SpinError=0, ParamsErr=0, PulsesErr=0;

		Parameters pset(fname);
		if(!pset.addSection("spins"))
		{
			std::cerr<<"Error: Solid..."<<std::endl;
			std::cerr<<" must define a 'spins{...}' section"<<std::endl;
			SpinError=1;Error=1;

		}
		SequenceRun myRunner;
		
		Parameters SpinPset;
		if(SpinError!=1)
		{
			try{
				SpinPset=(pset.section("spins"));
				myRunner.setSystems(SpinPset);
			}catch(BL_exception e){
				Error=1;
				e.print(std::cerr);
			}
		}

		if(!pset.addSection("parameters"))
		{
			std::cerr<<"Error: Solid..."<<std::endl;
			std::cerr<<" must define a 'parameters{...}' section"<<std::endl;
			ParamsErr=1; Error=1;

		}

		Parameters ParamsPset;

		if(SpinError!=1)
		{
			ParamsPset=(pset.section("parameters"));
			try{
				myRunner.setParams(ParamsPset);

			}catch(BL_exception e){
				Error=1;
				e.print(std::cerr);
			}

			try{
				ParamsPset=(pset.section("parameters"));
				myRunner.setPowders(ParamsPset);
			}catch(BL_exception e){
				Error=1;
				e.print(std::cerr);
			}
		}


	//get the spin sections
		if(!pset.addSection("pulses"))
		{
			std::cerr<<"Error: Solid..."<<std::endl;
			std::cerr<<" must define a 'pulses{...}' section"<<std::endl;
			PulsesErr=1; Error=1;

		}

		Parameters subPulse;

		if(PulsesErr!=1 && Error!=1)
		{
			try	{
				signal(SIGHUP, sigdie);
				signal(SIGQUIT , sigdie);
				signal(SIGKILL  , sigdie);
				signal(SIGINT  , sigdie);
				signal(SIGSEGV  , sigdie);
				signal(SIGABRT , sigdie);
				signal(SIGTERM , sigdie);
			

				subPulse=(pset.section("pulses"));
				myRunner.parsePulse(subPulse);
				myRunner.toRun=true;
				myRunner.run();
			}catch(BL_exception e){
				e.print(std::cerr);
				myRunner.dumpState();
				Error=1;
			}
		}

		//if(Error!=1) std::cout<<myRunner<<std::endl;
		if(MPIworld.master() ) std::cout<<std::endl;

		return MPIworld.end();
	//}

}
