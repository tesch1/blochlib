
/* mpi_config.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 01-20-02
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

#ifndef _BL_MPI_CONTROLLER_H_
#define _BL_MPI_CONTROLLER_H_ 1


#include <string>
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include "mpi/mpi_config.h"
#include "mpi/mpi_tools.h"

BEGIN_BL_NAMESPACE



//need to give this some sort of type for utility without MPI
//#ifndef HAVE_MPI
//	typedef MPI_Comm int;
//#endif

class MPIcontroller;

//This is used as our 'global' controller
extern MPIcontroller MPIworld;
extern int MasterControlRef;

class MPIcontroller{
	private:

		MPI_Comm communicator_;
		int rank_;
		int size_;
		int tag_;
		bool hasended_;


	public:
		MPIcontroller();
		MPIcontroller(const MPI_Comm &communin);
		MPIcontroller(int &argv, char *argc[]);
		MPIcontroller(const MPIcontroller &communin);

	//start up from a new input set
		int start(int &argv, char *argc[]);

	//start up from an MPI Communicator
		int start(const MPI_Comm &communin );

	//end the controlers and clean up mpi
		int end();

	//end the controlers and clean up mpi
		int stop();

	//aborts the MPI communcator...
		int abort();

	//basic MPI barrier command
		int barrier();

	//the destuctor is the same as 'end()'
		~MPIcontroller();

	//get the inerds of the class
		inline int rank()const {	return rank_;	}
		inline int size()const	{	return size_;	}
		inline int tag()const	{	return tag_;	}

		inline void advanceTag(){	(tag_>INT_MAX)?tag_=0:tag_++;	}

	//processor name
		std::string name();

		//inline MPI_Comm communicator(){	return communicator_;	}

	//returns true if the rank==0;
		inline bool master(){	return rank_==0;	}
		inline bool slave(){	return rank_!=0;	}

	//returns 'true' not running using MPI
		inline bool serial(){	return size_==1 && rank_==0;	}
		inline bool parallel(){	return size_>1;	}

	//get the communicator
		inline MPI_Comm communicator(){	return communicator_;	}

	//Loop splitter
		Range splitLoop(int &begin, int &end, int &div, int numproc=MPIsize);

	//Reductions
		template<class T>
		void reduce(T &thing, int ReduceType, int MasterProc=0)
		{		if(parallel()) MPIreduce(thing, ReduceType, MasterProc, communicator_);		}

	//Reductions
		template<class T>
		void reduce(T &thing, int ReduceType, int MasterProc, MPI_Comm comm)
		{		if(parallel()) MPIreduce(thing, ReduceType, MasterProc, comm);		}

	//reconstructions
		template<class T>
		void reconstruct(T &data, int begin, int end, int MasterProc=0)
		{		if(parallel()) MPIreconstruct(data, begin, end, MasterProc, communicator_);		}

	//reconstructions
		template<class T>
		void reconstruct(T &data, int begin, int end, int MasterProc, MPI_Comm comm)
		{		if(parallel()) MPIreconstruct(data, begin, end, MasterProc, comm);		}

	//scatter
		template<class T>
		void scatter(T &toscat, int MasterProc=0)
		{		if(parallel()) MPIscatter(toscat, MasterProc, communicator_);		}

	//scatter
		template<class T>
		void scatter(T &toscat, int MasterProc, MPI_Comm comm)
		{		if(parallel()) MPIscatter(toscat, MasterProc, comm);		}


	//the put functions to another proc
		template<class T>
		int put( T &out, int To)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIput(out, To, tag_, communicator_);
#else
			return 0;
#endif
		}

	//the put functions to another proc
		template<class T>
		int put( T &out, int To, int tagi)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIput(out, To, tagi, communicator_);
#else
			return 0;
#endif
		}

	//the put functions to another proc
		template<class T>
		int put( T &out, int To, int tagi, MPI_Comm comm)
		{
#ifdef HAVE_MPI
			RunTimeAssert(!hasended_);
			return MPIput(out, To, tagi, comm);
#else
			return 0;
#endif
		}

//Get from ANY source
// returns wehre it came from...
		template<class T>
		int getAny( T &out)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIgetAny(out, MPI_ANY_TAG, communicator_);
#else
			return 1;
#endif
		}

//Get from ANY source
// returns wehre it came from...
		template<class T>
		int getAny( T &out, int tag)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIgetAny(out, tag, communicator_);
#else
			return 1;
#endif
		}


	//the get functions
		template<class T>
		int get( T &out, int From)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIget(out, From, tag_, communicator_);
#else
			return 1;
#endif
		}

	//the get functions
		template<class T>
		int get( T &out, int From, int tagi)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIget(out, From, tagi, communicator_);
#else
			return 1;
#endif
		}

	//the get functions
		template<class T>
		int get( T &out, int From, int tagi, MPI_Comm comm)
		{
#ifdef HAVE_MPI
			RunTimeAssert(!hasended_);
			return MPIget(out, From, tagi, comm);
#else
			return 1;
#endif
		}

/****** NON BLOCKING GETS AND PUTS ****/
	//the put functions to another proc
		template<class T>
		int send( T &out, int To)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIsend(out, To, tag_, communicator_);
#else
			return 0;
#endif
		}

	//the put functions to another proc
		template<class T>
		int send( T &out, int To, int tagi)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIsend(out, To, tagi, communicator_);
#else
			return 0;
#endif
		}

	//the put functions to another proc
		template<class T>
		int send( T &out, int To, int tagi, MPI_Comm comm)
		{
#ifdef HAVE_MPI
			RunTimeAssert(!hasended_);
			return MPIsend(out, To, tagi, comm);
#else
			return 0;
#endif
		}

	//the get functions
		template<class T>
		int recv( T &out, int From)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIrecv(out, From, tag_, communicator_);
#else
			return 1;
#endif
		}

	//the get functions
		template<class T>
		int recv( T &out, int From, int tagi)
		{
#ifdef HAVE_MPI
			RunTimeAssert(communicator_!=0);
			RunTimeAssert(!hasended_);
			return MPIrecv(out, From, tagi, communicator_);
#else
			return 1;
#endif
		}

	//the get functions
		template<class T>
		int recv( T &out, int From, int tagi, MPI_Comm comm)
		{
#ifdef HAVE_MPI
			RunTimeAssert(!hasended_);
			return MPIrecv(out, From, tagi, comm);
#else
			return 1;
#endif
	}
};


END_BL_NAMESPACE


#endif


