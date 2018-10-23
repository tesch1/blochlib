
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

#ifndef _BL_MPI_CONFIG_CC_
#define _BL_MPI_CONFIG_CC_ 1

int MPIsize=1;
int MPIrank=0;



#ifndef ON_WINDOWS
#include "blochconfig.h"
#endif


#include <new>
#include <string>

#ifdef HAVE_MPI
	#include "mpi.h"

	#ifdef HAVE_CLIMITS
	 #include <climits>
	#else
	 #include <limits.h>
	#endif

#endif

#include "mpi/mpi_config.h"
#include "mpi/mpi_tools.h"

BEGIN_BL_NAMESPACE



/*
void MPIreduce(matrix &thing, int ReduceType, int MasterProc)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				matrix tmpVC;
				MPIget(tmpVC, i);
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


void MPIreduce(rmatrix &thing, int ReduceType, int MasterProc)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				rmatrix tmpVC;
				MPIget(tmpVC, i);
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

void MPIreduce(hmatrix &thing, int ReduceType, int MasterProc)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				hmatrix tmpVC;
				MPIget(tmpVC, i);
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

void MPIreduce(smatrix &thing, int ReduceType, int MasterProc)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				smatrix tmpVC;
				MPIget(tmpVC, i);
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

void MPIreduce(dmatrix &thing, int ReduceType, int MasterProc)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				dmatrix tmpVC;
				MPIget(tmpVC, i);
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

void MPIreduce(rdmatrix &thing, int ReduceType, int MasterProc)
{
#ifdef HAVE_MPI
	advanceMPItag();
	if(MPIrank!=MasterProc){
		MPIput(thing, MasterProc);
	}else{
		for(int i=0;i<MPIsize;++i)
		{
			if(i!=MasterProc){
				rdmatrix tmpVC;
				MPIget(tmpVC, i);
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
}*/

void MPIsync()
{
#ifdef HAVE_MPI
		char ch;
		MPI_Status status;
		if (MPIrank==0)
		{
			MPI_Send(&ch, 1, MPI_BYTE, 0, 1, MPI_COMM_WORLD);
			MPI_Recv(&ch, 1, MPI_BYTE, 0, 1, MPI_COMM_WORLD, &status);
			MPI_Send(&ch, 1, MPI_BYTE, 0, 1, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Recv(&ch, 1, MPI_BYTE, MPIrank, 1, MPI_COMM_WORLD, &status);
			MPI_Send(&ch, 1, MPI_BYTE, MPIrank, 1, MPI_COMM_WORLD);
			MPI_Recv(&ch, 1, MPI_BYTE, MPIrank, 1, MPI_COMM_WORLD, &status);
		}
#endif
}

Range MPIsplitLoop(int &begin, int &end, int &divs, int numproc)
{
	divs=1;
#ifdef HAVE_MPI
	if(MPIsize<=1) return Range(Range::Start, Range::End);
	RunTimeAssert(end>begin);
	int inend=end;
	divs=int(double(end-begin)/double(numproc));
	if(divs<=1){
	}else{
		begin=MPIrank*divs;
		end=(MPIrank+1)*divs;
	}
	while(MPIrank==MPIsize-1 && end < inend) end++;
#endif
	return Range(begin, end);
}

int MPItag=0;

//acts as the 'MPI_Init(..)' function
//that automatically sets MPIsize and MPIrank...
//returns false if something fails
bool MPIstart(int &argc, char *argv[])
{
#ifdef HAVE_MPI
	int err=MPI_Init(&argc, &argv);
  	err=MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);   /* Get my rank   */
  	err=MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);   /* Get the total */
	if(err>0){
		std::cerr<<std::endl<<"Error:: MPIstart()"<<std::endl;
		std::cerr<<" unknonwn error on start up "<<std::endl;
		return false;
	}
	return err;
#else
	return 0;
#endif
}

//acts as the 'MPI_Init(..)' function
//that automatically sets MPIsize and MPIrank...
//returns false if something fails
void MPIend()
{
#ifdef HAVE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
}

std::string MPIpname()
{
#ifdef HAVE_MPI
	int  namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Get_processor_name(processor_name,&namelen);
	return std::string(processor_name, namelen);
#else
	return "Serial Proc";
#endif
}

void advanceMPItag()
{
#ifdef HAVE_MPI
	if(MPItag>INT_MAX) MPItag=0;
	MPItag++;
#endif
}


END_BL_NAMESPACE



#endif


