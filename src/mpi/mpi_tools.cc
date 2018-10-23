
/* mpi_tools.h ********/


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
 	Several Send an Get functions with little data structures to make sending
 	and reciving messages a bit easier
 	This is at its infancy
 	USE WITH CAUTION
 */

#ifndef _BL_MPI_TOOLS_CC_
#define _BL_MPI_TOOLS_CC_ 1

#ifdef  HAVE_MPI
	#include "mpi.h"
#endif

#include <new>
#include <string>
#include "container/containers.h"
#include "mpi/mpi_tools.h"
#include "mpi/mpi_config.h"
#include "mpi/mpi_packer.h"

BEGIN_BL_NAMESPACE


/**** CHAR ****/
//puts an char....
int MPIput(const  char &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	int err=MPI_Send(const_cast<char *>(&tosend), 1, MPI_CHAR, To, tag, comm);
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'From'
int MPIget(char &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Status sts;
	int err=MPI_Recv(&toget, 1, MPI_CHAR, From, tag, comm, &sts);
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'Any Proc'
//returns where is came from
int MPIgetAny(char &toget, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Status sts;
	int err=MPI_Recv(&toget, 1, MPI_CHAR, MPI_ANY_SOURCE, tag, comm, &sts);
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}


//puts an char.... NONBLOCKING
int MPIsend(const  char &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Request sts;
	int err=MPI_Isend(const_cast<char *>(&tosend), 1, MPI_CHAR, To, tag, comm,&sts);
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'From' NONBLOCKING
int MPIrecv(char &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Request sts;
	int err=MPI_Irecv(&toget, 1, MPI_CHAR, From, tag, comm, &sts);
	return err;
#else
	return 0;
#endif
}


/***** INTEGERS ****/
//gets an int from processor 'From'
int MPIget(int &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Status sts;
	int err=MPI_Recv(&toget, 1, MPI_INT, From, tag, comm, &sts);
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'Any Proc'
// returns where is came from
int MPIgetAny(int &toget, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	MPI_Status sts;
	MPI_Recv(&toget, 1, MPI_INT, MPI_ANY_SOURCE, tag, comm, &sts);
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}

//puts an int.... NONBLOCKING
int MPIput(const  int &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	int err=MPI_Send(const_cast<int *>(&tosend), 1, MPI_INT, To, tag, comm);
	return err;
#else
	return 0;
#endif
}

//puts an int.... NONBLOCKING
int MPIsend(const  int &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Request sts;
	int err=MPI_Isend(const_cast<int *>(&tosend), 1, MPI_INT, To, tag, comm,&sts);
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'From' NONBLOCKING
int MPIrecv(int &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Request sts;
	int err=MPI_Irecv(&toget, 1, MPI_INT, From, tag, comm, &sts);
	return err;
#else
	return 0;
#endif
}

/******* FLOAT *****/

//puts an float....
int MPIput(const  float &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	int err=MPI_Send(const_cast<float *>(&tosend), 1, MPI_FLOAT, To, tag, comm);
	return err;
#else
	return 0;
#endif
}


//get an float....
int MPIget( float &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Status sts;
	int err=MPI_Recv(&tosend, 1, MPI_FLOAT, To, tag, comm, &sts);
	return err;
#else
	return 0;
#endif
}

//gets an float from processor 'Any Proc'
// returns where is came from
int MPIgetAny(float &toget, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	MPI_Status sts;
	MPI_Recv(&toget, 1, MPI_FLOAT, MPI_ANY_SOURCE, tag, comm, &sts);
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}

//puts an int.... NONBLOCKING
int MPIsend(const  float &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Request sts;
	int err=MPI_Isend(const_cast<float *>(&tosend), 1, MPI_FLOAT, To, tag, comm,&sts);
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'From' NONBLOCKING
int MPIrecv(float &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Request sts;
	int err=MPI_Irecv(&toget, 1, MPI_FLOAT, From, tag, comm, &sts);
	return err;
#else
	return 0;
#endif
}

/******* DOUBLE *****/
//gets an double from processor 'From'
int MPIget(double &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Status sts;
	int err=MPI_Recv(&toget, 1, MPI_DOUBLE, From, tag, comm, &sts);
	return err;
#else
	return 0;
#endif
}

//puts an double....
int MPIput(const  double &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	int err=MPI_Send(const_cast<double *>(&tosend), 1, MPI_DOUBLE, To, tag, comm);
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'Any Proc'
// returns where is came from
int MPIgetAny(double &toget, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	MPI_Status sts;
	MPI_Recv(&toget, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, comm, &sts);
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}

//puts an int.... NONBLOCKING
int MPIsend(const  double &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Request sts;
	int err=MPI_Isend(const_cast<double *>(&tosend), 1, MPI_DOUBLE, To, tag, comm,&sts);
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'From' NONBLOCKING
int MPIrecv(double &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	static MPI_Request sts;
	int err=MPI_Irecv(&toget, 1, MPI_DOUBLE, From, tag, comm, &sts);
	return err;
#else
	return 0;
#endif
}

/*** COMPLEX *****/

//puts an complex....

int MPIput(const Complex<double> &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	buffer=new char[Packer<DELEGATE, Complex<double> >::size(tosend) ];
	int pos=Packer<DELEGATE, Complex<double> >::pack(tosend, buffer);
	int err=MPI_Send(buffer, pos, MPI_CHAR, To, tag, comm);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets an double from processor 'From'
int MPIget(Complex<double> &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	static MPI_Status sts;
	int ss=Packer<DELEGATE, Complex<double> >::size(toget);
	buffer=new char[ss];
	int err=MPI_Recv(buffer, ss, MPI_CHAR, From, tag, comm, &sts);
	Packer<DELEGATE, Complex<double> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'Any Proc'
// returns where is came from
int MPIgetAny(Complex<double> &toget, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Status sts;
	int ss=Packer<DELEGATE, Complex<double> >::size(toget);
	buffer=new char[ss];
	MPI_Recv(buffer, ss, MPI_CHAR, MPI_ANY_SOURCE, tag, comm, &sts);
	Packer<DELEGATE, Complex<double> >::unpack(toget, buffer);
	delete [] buffer;
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}

//puts an complex....
int MPIsend(const Complex<double> &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Request sts;
	buffer=new char[Packer<DELEGATE, Complex<double> >::size(tosend) ];
	int pos=Packer<DELEGATE, Complex<double> >::pack(tosend, buffer);
	int err=MPI_Isend(buffer, pos, MPI_CHAR, To, tag, comm,&sts);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets an double from processor 'From'
int MPIrecv(Complex<double> &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Request sts;
	int ss=Packer<DELEGATE, Complex<double> >::size(toget);
	buffer=new char[ss];
	int err=MPI_Irecv(buffer, ss, MPI_CHAR, From, tag, comm, &sts);
	Packer<DELEGATE, Complex<double> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

/**** Complex<float> ***/

int MPIput(const Complex<float> &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	buffer=new char[Packer<DELEGATE, Complex<float> >::size(tosend) ];
	int pos=Packer<DELEGATE, Complex<float> >::pack(tosend, buffer);
	int err=MPI_Send(buffer, pos, MPI_CHAR, To, tag, comm);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets an double from processor 'From'
int MPIget(Complex<float> &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	static MPI_Status sts;
	int ss=Packer<DELEGATE, Complex<float> >::size(toget);
	buffer=new char[ss];
	int err=MPI_Recv(buffer, ss, MPI_CHAR, From, tag, comm, &sts);
	Packer<DELEGATE, Complex<float> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets an int from processor 'Any Proc'
// returns where is came from
int MPIgetAny(Complex<float> &toget, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Status sts;
	int ss=Packer<DELEGATE, Complex<float> >::size(toget);
	buffer=new char[ss];
	MPI_Recv(buffer, ss, MPI_CHAR, MPI_ANY_SOURCE, tag, comm, &sts);
	Packer<DELEGATE, Complex<float> >::unpack(toget, buffer);
	delete [] buffer;
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}

//puts an complex....
int MPIsend(const Complex<float> &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Request sts;
	buffer=new char[Packer<DELEGATE, Complex<float> >::size(tosend) ];
	int pos=Packer<DELEGATE, Complex<float> >::pack(tosend, buffer);
	int err=MPI_Isend(buffer, pos, MPI_CHAR, To, tag, comm,&sts);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets an double from processor 'From'
int MPIrecv(Complex<float> &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Request sts;
	int ss=Packer<DELEGATE, Complex<float> >::size(toget);
	buffer=new char[ss];
	int err=MPI_Irecv(buffer, ss, MPI_CHAR, From, tag, comm, &sts);
	Packer<DELEGATE, Complex<float> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

/*** STRING ***/

//puts a string....
int MPIput(const std::string &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	int ss=Packer<DELEGATE, std::string >::size(tosend);
	char *buffer;
	buffer=new char[ss];
	int pos=Packer<DELEGATE, std::string >::pack(tosend, buffer);
	int err=MPI_Send(&ss, 1, MPI_INT, To, tag, comm);
	err=MPI_Send(buffer, pos, MPI_CHAR, To, tag+1, comm);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}


//gets a sting from processor 'From'
int MPIget(std::string &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	static MPI_Status sts;
	int err=MPI_Recv(&ss, 1, MPI_INT, From, tag, comm, &sts);
	buffer=new char[ss];
	err=MPI_Recv(buffer, ss, MPI_CHAR, From, tag+1, comm, &sts);
	toget.resize(ss);
	Packer<DELEGATE, std::string >::unpack(toget, buffer);

	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets a sting from processor 'ANY PROC'
//retuns where it came from
int MPIgetAny(std::string &toget,  int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	MPI_Status sts;
	MPI_Recv(&ss, 1, MPI_INT, MPI_ANY_SOURCE, tag, comm, &sts);
	buffer=new char[ss];
	MPI_Recv(buffer, ss, MPI_CHAR, MPI_ANY_SOURCE, tag, comm, &sts);
	toget.resize(ss);
	Packer<DELEGATE, std::string >::unpack(toget, buffer);

	delete [] buffer;
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}


//puts a string....
int MPIsend(const std::string &tosend, int To, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	int ss=Packer<DELEGATE, std::string >::size(tosend);
	char *buffer;
	buffer=new char[ss];
	static MPI_Request sts;
	int pos=Packer<DELEGATE, std::string >::pack(tosend, buffer);
	int err=MPI_Isend(&ss, 1, MPI_INT, To, tag, comm,&sts);
	err=MPI_Isend(buffer, pos, MPI_CHAR, To, tag+1, comm,&sts);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}


//gets a sting from processor 'From'
int MPIrecv(std::string &toget, int From, int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	static MPI_Request sts;
	int err=MPI_Irecv(&ss, 1, MPI_INT, From, tag, comm, &sts);
	buffer=new char[ss];
	err=MPI_Irecv(buffer, ss, MPI_CHAR, From, tag+1, comm, &sts);
	toget.resize(ss);
	Packer<DELEGATE, std::string >::unpack(toget, buffer);

	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

END_BL_NAMESPACE



#endif


