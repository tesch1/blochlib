
/* mpi_config.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 04-20-02
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
 	this holds the simple constants and global functions
 	for use with MPI and BLochlib...this is only a VERY
 	basic interface to the MPI pieces...and is at its infancy
 	USE WITH CAUTION
 */

#ifndef _BL_MPI_CONTROLLER_CC_
#define _BL_MPI_CONTROLLER_CC_ 1

#include <signal.h>

#ifndef ON_WINDOWS
#include "blochconfig.h"
#endif


#ifdef HAVE_CLIMITS
 #include <climits>
#else
 #include <limits.h>
#endif

#ifdef HAVE_MPI
	#include "mpi.h"
#endif

#include <string>
#include "mpi/mpi_config.h"
#include "mpi/mpi_controller.h"
#include "mpi/mpi_tools.h"

BEGIN_BL_NAMESPACE


//Our ever present master of MPI
MPIcontroller MPIworld;
int StartTag=1000;
int MasterControlRef=0;

MPIcontroller::MPIcontroller(int &argv, char *argc[])
{		start(argv, argc);	}

MPIcontroller::MPIcontroller(const MPI_Comm &commin)
{		start(commin);	}

MPIcontroller::MPIcontroller()
{
	communicator_=0;
	rank_=0;
	size_=1;
	tag_=StartTag;
	hasended_=true;
}

MPIcontroller::MPIcontroller(const MPIcontroller &dup)
{
/*	rank_=dup.rank_;
	size_=dup.size_;
	tag_=dup.tag_;
	hasended_=dup.hasended_;
#ifdef HAVE_MPI
	MPI_Comm_dup(communicator_, const_cast<MPI_Comm *>(&dup.communicator_));
#endif*/
	if(dup.communicator_){
		start(dup.communicator_);
	}else{
		rank_=0;
		size_=1;
		tag_=StartTag;
		hasended_=true;
	}
}

int MPIcontroller::start(const MPI_Comm &commin)
{
#ifdef HAVE_MPI
	int err=MPI_Comm_rank(commin, &rank_);   /* Get my rank   */
	err=MPI_Comm_size(commin, &size_);   /* Get the total */
	MPI_Comm_dup(commin, &communicator_);
	if(err>0){
		std::cerr<<std::endl<<"Error:: MPIcontroller::start()"<<std::endl;
		std::cerr<<" unknonwn error on start up, setting mode for serial running"<<std::endl;
		rank_=0;
		size_=1;
	}
	tag_=StartTag;
	hasended_=false;
	++MasterControlRef;
	//cout<<"REFED: "<<MasterControlRef<<" rank: "<<rank()<<endl;
	return err;
#else
	tag_=StartTag;
	rank_=0;
	size_=1;
	hasended_=true;
	MasterControlRef=1;
    return 1;
#endif
}

int MPIcontroller::start(int &argv, char *argc[])
{
#ifdef HAVE_MPI
	int err=MPIstart(argv, argc);
	start(MPI_COMM_WORLD);
	/*err=MPI_Comm_rank(MPI_COMM_WORLD, &rank_);   // Get my rank
	err=MPI_Comm_size(MPI_COMM_WORLD, &size_);   // Get the total
	MPI_Comm_dup(MPI_COMM_WORLD, &communicator_);
	//communicator_=MPI_COMM_WORLD;
	if(err>0){
		std::cerr<<std::endl<<"Error:: MPIcontroller::start()"<<std::endl;
		std::cerr<<" unknonwn error on start up, setting mode for serial running"<<std::endl;
		rank_=0;
		size_=1;
	}
	hasended_=false;
	tag_=StartTag;*/
	return err;
#else
	tag_=StartTag;
	rank_=0;
	hasended_=true;
	size_=1;
	MasterControlRef=1;
    return 1;
#endif
}

//processor name
std::string MPIcontroller::name()
{
#ifdef HAVE_MPI
	int  namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Get_processor_name(processor_name,&namelen);
	return std::string(processor_name, namelen);
#else
	return "";
#endif

}

//Loop splitter
Range MPIcontroller::splitLoop(int &begin, int &end, int &div, int numproc)
{
	return MPIsplitLoop(begin, end, div, numproc);
}

int MPIcontroller::abort()
{
#ifdef HAVE_MPI
	int flag=0;
	MPI_Finalized(&flag);
		//cout<<endl<<"NUMBER of REfs: "<<MasterControlRef<<" rank: "<<rank()<<" finalized: "<<flag<<endl;
	if(!flag){
		MPI_Abort(communicator(), SIGINT);
	}
	return 1;
#else
	return 1;
#endif
}

int MPIcontroller::end()
{
#ifdef HAVE_MPI
	--MasterControlRef;
	int flag=1;
	if(!hasended_) MPI_Finalized(&flag);
	//cout<<endl<<"NUMBER of REfs: "<<MasterControlRef<<" rank: "<<rank()<<" finalized: "<<flag<<endl;
	if(!flag){
//		cout<<endl<<"NUMBER of REfs: "<<MasterControlRef<<" rank: "<<rank()<<" finalized: "<<flag<<endl;
		hasended_=true;
		return MPI_Finalize();		//kill everything MPI realted
	}else{
		return 1;
	}

#else
	return 1;
#endif
}

int MPIcontroller::stop(){	return end();	}

int MPIcontroller::barrier()
{
#ifdef HAVE_MPI
	if(size_!=1){
		return MPI_Barrier(communicator_); //make sure everything gets to this point
	}
	return 0;
#else
	return 1;
#endif
}

MPIcontroller::~MPIcontroller()
{
#ifdef HAVE_MPI
	/*if(!hasended_ && &communicator_)
	{
		std::cout<<communicator_<<" "<<&communicator_<<" "<<hasended_<<std::endl;
		--MasterControlRef;
		MPI_Comm_free(&communicator_);
		std::cout<<communicator_<<" "<<&communicator_<<" "<<hasended_<<std::endl;

	}*/
#endif
}//	end();	}


END_BL_NAMESPACE


#endif


