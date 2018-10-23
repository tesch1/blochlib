
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

#ifndef _BL_MPI_CONFIG_H_
#define _BL_MPI_CONFIG_H_ 1


extern int MPIsize;
extern int MPIrank;


#include <new>
#include <string>
#include "container/containers.h"

BEGIN_BL_NAMESPACE


extern int MPItag;

#ifdef HAVE_MPI
	#include "mpi.h"
#else
/*	typedef int MPI_Status;
	typedef int MPI_COMM_WORLD;
	inline int MPI_Init(int *, char ***){ return 0; }
	inline int MPI_Comm_rank(int, int &){ return 0;	}
	inline int MPI_Comm_size(int, int &){ return 0;	}
	template<class T>
	inline int MPI_Recv(T &, int, int, int, int, int, MPI_Status){	return 0;	}
	template<class T>
	inline int MPI_Send(T &, int, int, int, int, MPI_Status){	return 0;	}
	inline int MPI_Finalize(){	return 0;	}
	inline int MPI_Barrier(int){	return 0;	}
	void MPI_Get_processor_name(char *,int *);
*/
	#define MPI_Comm int
	#define MPI_COMM_WORLD 1
	#define MPI_ANY_TAG 1
	#define MPI_ANY_SOURCE 1
#endif



void MPIsync();

//defined OUTside the HVE_MPI becuase we want to include it regardless of
//MPI presesnce or not
//a function that will return the 'Range' splitter
// begin=the start value to split
// end=the end value to split
// div=(end-begin)/MPIsize
// numproc=MPIsize (default)
//NOTE** the initial values of begin and 'end' MUST be those
//as if there was NO loop splitting!!
Range MPIsplitLoop(int &begin, int &end, int &div, int numproc=MPIsize);

class Reduce{
	public:
		enum{
			Add		=0x00001,
			Multiply=0x00002,
			Subtract=0x00004,
			Divide	=0x00008,
			Max		=0x00010,
			Min		=0x00020
		} ops;
};

//defined OUTside the HVE_MPI becuase we want to include it regardless of
//MPI presesnce or not
//a function that will 'reduce' a vector from all the procs by a type
// assumes that all the incoming 'T's are going to the main Proc 'MainProc'
// whose default is '0' and that they come from all
// begin=the start value to split
// end=the end value to split
// div=(end-begin)/MPIsize
// numproc=MPIsize (default)

/*//special case for matrix as the 'mat/=mat' is not defined but mat=mat/mat is
void MPIreduce(matrix &thing, int ReduceType, int MasterProc);
void MPIreduce(rmatrix &thing, int ReduceType, int MasterProc);
void MPIreduce(rdmatrix &thing, int ReduceType, int MasterProc);
void MPIreduce(hmatrix &thing, int ReduceType, int MasterProc);
void MPIreduce(smatrix &thing, int ReduceType, int MasterProc);
void MPIreduce(dmatrix &thing, int ReduceType, int MasterProc);
*/
template<class T>
void MPIreduce(T &thing, int ReduceType, int MasterProc=0, MPI_Comm comm=MPI_COMM_WORLD);

//special one for Coord<T,N> where 'min and max' are not defined
template<class T, int N>
void MPIreduce(coord<T,N> &thing, int ReduceType, int MasterProc=0, MPI_Comm comm=MPI_COMM_WORLD);

//special one for VEctor<Coord<T,N>> where 'min and max' are not defined
template<class T, int N>
void MPIreduce(Vector<coord<T,N> > &thing, int ReduceType, int MasterProc=0, MPI_Comm comm=MPI_COMM_WORLD);

template<class T, class struct_t>
void MPIreduce(_matrix<T, struct_t> &thing, int ReduceType, int MasterProc=0, MPI_Comm comm=MPI_COMM_WORLD);

//these function reduce the diesired variable then scatter it to
//every proc
template<class T>
void MPIreduceAndScatter(T &thing, int ReduceType, int MasterProc=0);

//This function simply scatters the variable to every proc
//The 'MasterProc' is the one with the 'correct' data
template<class T>
void MPIscatter(T &toscat, int MasterProc=0, MPI_Comm comm=MPI_COMM_WORLD);

//this function first sends to the 'master cpu' the appropriate chunk of a vector
//as well its 'range' to recollect
//then reconstructs the master vector from each of these peices
//on the master proc, then sends out the reconstructed vector BACK to
// all of the procs, so all the procs have the full vector
// it will take in the vector to reconstruct, and the 'begin' and 'end'
//as generated from 'MPIsplitLoop()'
template<class T>
void MPIreconstruct(Vector<T> &data, int begin, int end, int MasterProc=0, MPI_Comm comm=MPI_COMM_WORLD);

//will also scatter the result
template<class T>
void MPIreconstructAndScatter(Vector<T> &data, int begin, int end, int MasterProc=0);


/********** implimentations of Templated functions *************/

//specail function for matrix where 'Max and Min' are
// not really defined
template<class T, class struct_t>
void MPIreduce(_matrix<T, struct_t> &thing, int ReduceType, int MasterProc, MPI_Comm comm)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc, comm);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				_matrix<T, struct_t> tmpVC;
				MPIget(tmpVC, i,comm);
				switch(ReduceType){
					case Reduce::Add:
						thing+=tmpVC; break;
					case Reduce::Multiply:
						thing*=tmpVC; break;
					case Reduce::Divide:
						thing=thing/tmpVC; break;
					case Reduce::Subtract:
						thing-=tmpVC; break;
					default: break;
				}
			}
		}
	}
#endif
}


template<class T>
void MPIreduce(T &thing, int ReduceType, int MasterProc, MPI_Comm comm)
{

#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc, comm);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				T tmpVC;
				MPIget(tmpVC, i, comm);
				switch(ReduceType){
					case Reduce::Add:
						thing+=tmpVC; break;
					case Reduce::Multiply:
						thing*=tmpVC; break;
					case Reduce::Divide:
						thing/=tmpVC; break;
					case Reduce::Subtract:
						thing-=tmpVC; break;
					case Reduce::Max:
						thing=max(tmpVC, thing); break;
					case Reduce::Min:
						thing=min(tmpVC, thing); break;
					default: break;
				}
			}
		}
	}

#endif
}

//specail function for coord<> and VEctor<coord<> > where 'Max and Min' are
// not really defined
template<class T, int N>
void MPIreduce(coord<T,N> &thing, int ReduceType, int MasterProc, MPI_Comm comm)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc,MPItag, comm);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				coord<T,N> tmpVC;
				MPIget(tmpVC, i, MPItag, comm);
				switch(ReduceType){
					case Reduce::Add:
						thing+=tmpVC; break;
					case Reduce::Multiply:
						thing*=tmpVC; break;
					case Reduce::Divide:
						thing/=tmpVC; break;
					case Reduce::Subtract:
						thing-=tmpVC; break;
					default: break;
				}
			}
		}
	}
#endif
}

//specail function for coord<> and Vector<coord<> > where 'Max and Min' are
// not really defined
template<class T, int N>
void MPIreduce(Vector<coord<T,N> > &thing, int ReduceType, int MasterProc, MPI_Comm comm)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc, MPItag, comm);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				Vector<coord<T,N> > tmpVC;
				MPIget(tmpVC, i, MPItag, comm);
				switch(ReduceType){
					case Reduce::Add:
						thing+=tmpVC; break;
					case Reduce::Multiply:
						thing*=tmpVC; break;
					case Reduce::Divide:
						thing/=tmpVC; break;
					case Reduce::Subtract:
						thing-=tmpVC; break;
					default: break;
				}
			}
		}
	}
#endif
}



// The scatter function
template<class T>
void MPIscatter(T &toscat, int MasterProc, MPI_Comm comm)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank==MasterProc){
		for(int i=0;i<MPIsize;++i){
			if(i!=MasterProc)	MPIput(toscat, i, MPItag, comm);
		}
	}else{
		MPIget(toscat, MasterProc, MPItag, comm);
	}
	MPI_Barrier(comm);
#endif
}

//these function reduce the diesired variable then scatter it to
//every proc
template<class T>
void MPIreduceAndScatter(T &thing, int ReduceType, int MasterProc)
{
	MPIreduce(thing, ReduceType, MasterProc);
	MPIscatter(thing, MasterProc);
}

//this function first sends to the 'master cpu' the appropriate chunk of a vector
//as well its 'range' to recollect
//then reconstructs the master vector from each of these peices
//on the master proc, then sends out the reconstructed vector BACK to
// all of the procs, so all the procs have the full vector
// it will take in the vector to reconstruct, and the 'begin' and 'end'
//as generated from 'MPIsplitLoop()'
template<class T>
void MPIreconstruct(Vector<T> &data, int begin, int end, int MasterProc, MPI_Comm comm)
{
#ifdef HAVE_MPI
	RunTimeAssert(end>begin && data.size()>=end);
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(begin, MasterProc, MPItag, comm);
		MPIput(end, MasterProc, MPItag+1, comm);
		Vector<T> tosend(data(Range(begin,end-1)));
		MPIput(tosend, MasterProc, MPItag+2, comm);
	}else{
		for(int i=0;i<MPIsize;++i){
			if(i!=MasterProc){
				Vector<T> tmpVC;
				int rB=0, rE=0;
				MPIget(rB, i, MPItag, comm);
				MPIget(rE, i, MPItag+1, comm);
				MPIget(tmpVC, i,MPItag+2, comm);
				data(Range(rB,rE-1))=tmpVC;
			}
		}
	}

#endif
}

template<class T>
void MPIreconstruct(Vector<T> &data, Range R, int MasterProc, MPI_Comm comm)
{
#ifdef HAVE_MPI
	RunTimeAssert(end>begin && data.size()>=end);
	advanceMPItag();
	if(MPIrank!=MasterProc){
		int b=R.first(0), endd=R.last(data.size());
		MPIput(b, MasterProc, MPItag, comm);
		MPIput(endd, MasterProc, MPItag+1, comm);
		Vector<T> tosend(data(Range(begin,end-1)));
		MPIput(tosend, MasterProc, MPItag+2, comm);
	}else{
		for(int i=0;i<MPIsize;++i){
			if(i!=MasterProc){
				Vector<T> tmpVC;
				int rB=0, rE=0;
				MPIget(rB, i, MPItag, comm);
				MPIget(rE, i, MPItag+1, comm);
				MPIget(tmpVC, i,MPItag+2, comm);
				data(Range(rB,rE-1))=tmpVC;
			}
		}
	}

#endif
}

template<class T>
void MPIreconstructAndScatter(Vector<T> &data, int begin, int end, int MasterProc)
{
	MPIreconstruct(data, begin, end, MasterProc);
	MPIscatter(data, MasterProc);
}


//acts as the 'MPI_Init(..)' function
//that automatically sets MPIsize and MPIrank...
bool MPIstart(int &argv, char *argc[]);

//acts as the 'MPI_Init(..)' function
//that automatically sets MPIsize and MPIrank...
void MPIend();

//returns the processor Name
std::string MPIpname();

//advances the MPI tag to the next value
void advanceMPItag();


END_BL_NAMESPACE


#endif


