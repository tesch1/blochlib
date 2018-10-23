
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
 	This is at its infancy...	USE WITH CAUTION
 	the functions here are a simple case of massive function overloading

 */

#ifndef _BL_MPI_TOOLS_H_
#define _BL_MPI_TOOLS_H_ 1

#ifdef HAVE_MPI
	#include "mpi.h"
#endif

#include <new>
#include <string>
#include "container/containers.h"
#include "mpi/mpi_packer.h"
#include "mpi/mpi_config.h"

BEGIN_BL_NAMESPACE


//puts a char....
int MPIput(const char &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a char from processor 'From'
int MPIget(char &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a char from processor ANy Proc
//returns where it came from
int MPIgetAny(char &toget, int Tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD );
//puts a char.... NONBLOCKING
int MPIsend(const char &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a char from processor 'From' NONBLOCKING
int MPIrecv(char &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );


//puts a int....
int MPIput(const  int &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a int from processor 'From'
int MPIget(int &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a int from processor ANy Proc
//returns where it came from
int MPIgetAny(int &toget, int Tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD );
//puts a int....NONBLOCKING
int MPIsend(const  int &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a int from processor 'From' NONBLOCKING
int MPIrecv(int &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );

//puts a double....
int MPIput(const double &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a double from processor 'From'
int MPIget(double &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a double from processor ANy Proc
//returns where it came from
int MPIgetAny(double &toget, int Tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD );
//puts a double....NONBLOCKING
int MPIsend(const  double &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a double from processor 'From'...NONBLOCKING
int MPIrecv(double &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );

//puts a complex....
int MPIput(const  Complex<double> &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a complex from processor 'From'
int MPIget(Complex<double> &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a complex from processor ANy Proc
//returns where it came from
int MPIgetAny(Complex<double> &toget, int Tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD );
//puts a complex........NONBLOCKING
int MPIsend(const  Complex<double> &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a complex from processor 'From'....NONBLOCKING
int MPIrecv(Complex<double> &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );

int MPIput(const  Complex<float> &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a complex from processor 'From'
int MPIget(Complex<float> &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a complex from processor ANy Proc
//returns where it came from
int MPIgetAny(Complex<float> &toget, int Tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD );
//puts a complex........NONBLOCKING
int MPIsend(const  Complex<float> &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a complex from processor 'From'....NONBLOCKING
int MPIrecv(Complex<float> &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );


//puts a float....
int MPIput(const  float &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a float from processor 'From'
int MPIget(float &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a float from processor ANy Proc
//returns where it came from
int MPIgetAny(float &toget, int Tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD );
//puts a float.........NONBLOCKING
int MPIsend(const  float &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a float from processor 'From'.....NONBLOCKING
int MPIrecv(float &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );

//puts a string....
int MPIput(const std::string &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a sting from processor 'From'
int MPIget(std::string &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a string from processor ANy Proc
//returns where it came from
int MPIgetAny(std::string &toget, int Tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD );
//puts a string........NONBLOCKING
int MPIsend(const std::string &tosend, int To, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );
//gets a sting from processor 'From'....NONBLOCKING
int MPIrecv(std::string &toget, int From, int Tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD );

/**** VECTOR ****/
//puts a Vector<T>....
template<class T>
int MPIput(const Vector<T> &tosend, int To, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	int ss=Packer<DELEGATE, Vector<T> >::size(tosend);
	char *buffer;
	buffer=new char[ss];
	int position=Packer<DELEGATE, Vector<T> >::pack(tosend, buffer);
	int err=MPI_Send(&ss, 1, MPI_INT, To, tag, comm);
	err=MPI_Send(buffer, position, MPI_CHAR, To, tag+1, comm);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}


//gets a Vector<T>
template<class T>
int MPIget(Vector<T> &toget, int From, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	MPI_Status sts;
	int err=MPI_Recv(&ss, 1, MPI_INT, From, tag, comm, &sts);
	buffer=new char[ss];
	err=MPI_Recv(buffer, ss, MPI_CHAR, From, tag+1, comm, &sts);
	Packer<DELEGATE, Vector<T> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets a Vector<T> from ANY source
//returns where it came from
template<class T>
int MPIgetAny(Vector<T> &toget,  int tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	MPI_Status sts;
	int err=MPI_Recv(&ss, 1, MPI_INT, MPI_ANY_SOURCE, tag, comm, &sts);
	buffer=new char[ss];
	err=MPI_Recv(buffer, ss, MPI_CHAR, MPI_ANY_SOURCE, tag, comm, &sts);
	Packer<DELEGATE, Vector<T> >::unpack(toget, buffer);
	delete [] buffer;
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}

//send a vector NONBLOCKING
template<class T>
int MPIsend(const Vector<T> &tosend, int To, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	int ss=Packer<DELEGATE, Vector<T> >::size(tosend);
	char *buffer;
	buffer=new char[ss];
	MPI_Request req;
	int position=Packer<DELEGATE, Vector<T> >::pack(tosend, buffer);
	int err=MPI_Isend(&ss, 1, MPI_INT, To, tag, comm, &req);
	err=MPI_Isend(buffer, position, MPI_CHAR, To, tag+1, comm, &req);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}


//gets a Vector<T> NONBLOCKING
template<class T>
int MPIrecv( Vector<T> &toget, int From, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	MPI_Request sts;
	int err=MPI_Irecv(&ss, 1, MPI_INT, From, tag, comm, &sts);
	buffer=new char[ss];
	err=MPI_Irecv(buffer, ss, MPI_CHAR, From, tag+1, comm, &sts);
	Packer<DELEGATE, Vector<T> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

/*** COORD ***/

//put a coord<T,N>
template<class T, int N>
int MPIput(const coord<T,N> &tosend, int To, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	int ss=Packer<DELEGATE, coord<T,N> >::size(tosend);
	char *buffer;
	buffer=new char[ss];
	int position=Packer<DELEGATE, coord<T,N> >::pack(tosend, buffer);
	int err=MPI_Send(buffer, position, MPI_CHAR, To, tag, comm);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}


//gets a coord<T,N>
template<class T, int N>
int MPIget(coord<T,N> &toget, int From, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Status sts;
	int ss=Packer<DELEGATE, coord<T,N> >::size(toget);
	buffer=new char[ss];
	int err=MPI_Recv(buffer, ss, MPI_CHAR, From, tag, comm, &sts);
	Packer<DELEGATE, coord<T,N> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//gets a coord<T,N> from NAY source
//returns where it came from
template<class T, int N>
int MPIgetAny(coord<T,N> &toget, int tag=MPI_ANY_TAG, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Status sts;
	int ss=Packer<DELEGATE, coord<T,N> >::size(toget);
	buffer=new char[ss];
	int err=MPI_Recv(buffer, ss, MPI_CHAR, MPI_ANY_SOURCE, tag, comm, &sts);
	Packer<DELEGATE, coord<T,N> >::unpack(toget, buffer);
	delete [] buffer;
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}



//put a coord<T,N> NONBLOCKING
template<class T, int N>
int MPIsend(const coord<T,N> &tosend, int To, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	int ss=Packer<DELEGATE, coord<T,N> >::size(tosend);
	char *buffer;
	MPI_Request req;
	buffer=new char[ss];
	int position=Packer<DELEGATE, coord<T,N> >::pack(tosend, buffer);
	int err=MPI_Isend(buffer, position, MPI_CHAR, To, tag, comm, &req);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}


//gets a coord<T,N> NONBLOCKING
template<class T, int N>
int MPIrecv(coord<T,N> &toget, int From, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	MPI_Request sts;
	int ss=Packer<DELEGATE, coord<T,N> >::size(toget);
	buffer=new char[ss];
	int err=MPI_Irecv(buffer, ss, MPI_CHAR, From, tag, comm, &sts);
	Packer<DELEGATE, coord<T,N> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

/**** MATRIX ****/

//put a _matrix<T,Struct>
template<class Num_t, class Struct_T>
int MPIput(const _matrix<Num_t,Struct_T> &tosend, int To, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	int ss=Packer<DELEGATE, _matrix<Num_t,Struct_T> >::size(tosend);
	char *buffer;
	buffer=new char[ss];
	int position=Packer<DELEGATE, _matrix<Num_t,Struct_T> >::pack(tosend, buffer);
	int err=MPI_Send(&ss, 1, MPI_INT, To, tag, comm);
	err=MPI_Send(buffer, position, MPI_CHAR, To, tag+1, comm);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}


//get a _matrix<T,Struct>
template<class Num_t, class Struct_T>
int MPIget(_matrix<Num_t,Struct_T>  &toget, int From, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	MPI_Status sts;
	int err=MPI_Recv(&ss, 1, MPI_INT, From, tag, comm, &sts);
	buffer=new char[ss];
	err=MPI_Recv(buffer, ss, MPI_CHAR, From, tag+1, comm, &sts);
	Packer<DELEGATE, _matrix<Num_t,Struct_T> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

//get a _matrix<T,Struct> from ANY source
//returns where it came from
template<class Num_t, class Struct_T>
int MPIgetAny(_matrix<Num_t,Struct_T>  &toget, int tag=MPI_ANY_SOURCE, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	MPI_Status sts;
	int err=MPI_Recv(&ss, 1, MPI_INT, MPI_ANY_SOURCE, tag, comm, &sts);
	buffer=new char[ss];
	err=MPI_Recv(buffer, ss, MPI_CHAR, MPI_ANY_SOURCE, tag, comm, &sts);
	Packer<DELEGATE, _matrix<Num_t,Struct_T> >::unpack(toget, buffer);
	delete [] buffer;
	return sts.MPI_SOURCE;
#else
	return 0;
#endif
}

//put a _matrix<T,Struct> NONBLOCKING
template<class Num_t, class Struct_T>
int MPIsend(const _matrix<Num_t,Struct_T> &tosend, int To, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	int ss=Packer<DELEGATE, _matrix<Num_t,Struct_T> >::size(tosend);
	char *buffer;
	MPI_Request sts;
	buffer=new char[ss];
	int position=Packer<DELEGATE, _matrix<Num_t,Struct_T> >::pack(tosend, buffer);
	int err=MPI_Isend(&ss, 1, MPI_INT, To, tag, comm);
	err=MPI_Isend(buffer, position, MPI_CHAR, To, tag+1, comm, &sts);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}


//put a _matrix<T,Struct> NONBLOCKING
template<class Num_t, class Struct_T>
int MPIrecv(_matrix<Num_t,Struct_T>  &toget, int From, int tag=MPItag, MPI_Comm comm=MPI_COMM_WORLD )
{
#ifdef HAVE_MPI
	char *buffer;
	int ss;
	MPI_Request sts;
	int err=MPI_Irecv(&ss, 1, MPI_INT, From, tag, comm, &sts);
	buffer=new char[ss];
	err=MPI_Irecv(buffer, ss, MPI_CHAR, From, tag+1, comm, &sts);
	Packer<DELEGATE, _matrix<Num_t,Struct_T> >::unpack(toget, buffer);
	delete [] buffer;
	return err;
#else
	return 0;
#endif
}

END_BL_NAMESPACE


#endif


