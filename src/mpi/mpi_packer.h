
/* mpi_packer.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 01-24-02
 * -
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
 	mpi_packer.h-->this is essentially adapted from the "Serilaize.h" file of
 	"Cheeta (v1.1.1)"...becuase 'cheeta' kept giving me odd errors during runtime,
 	i have decided to right/rewrite some of its classes here and hopefully
 	debug some of its code...

 	much thanks to 'Cheeta' (http://www.acl.lanl.gov/cheeta/) for some helpful
 	hints here....

 	--> below you will find the "**CHEETA**" tag whenever a function
 	was simply coppied from Cheeta (most of this is...)
 */

#ifndef _MPI_PACKER_H_
#define _MPI_PACKER_H_ 1

#ifdef HAVE_MPI

#include <new>
#include <string>
#include "container/containers.h"

BEGIN_BL_NAMESPACE


/*----------------------------------------------------------------------

 Packer attempts to make the passing of complex C++ objects through MPI
 a bit easier to construct...such that eventually you can use a command
 like "put(MPI_COMM, Object, ToProc)" and not worry about the packing of
 that object...and then of course perform a "get(MPI_COMM, Object, FromProc)"
*/


/*
 To write your own versions of 'Packer' you need to be able to write
 these 4 functions for that object

    Return the storage needed to pack the item of type T
   static int size(const T &item);

    Pack an item of type T into the given buffer.  Return space used.
   static int pack(const T &item, char *buffer);

    Unpack an item of type T from the given buffer.  Set the given
    pointer to point at this item.  Return bytes unpacked.
   static int unpack(T* &p, char *buffer);

    Delete the item pointed to by the given pointer, that was
    unpacked with a previous call to unpack().
   static void cleanup(T *p);

----------------------------------------------------------------------
*/

//**CHEETA**----------------------------------------------------------------------
// Returns padding necessary for word alignment.
//----------------------------------------------------------------------

static inline int padding(int size);

static inline int padding(int size)
{
  int extra = size % sizeof(void*);
  return (extra == 0) ? 0 : sizeof(void*) - extra;
}



// The general tag type used to specialize the template "Packer" Below
// this is meant for 'single' objects (Arraies are below)

struct SINGLE
{
  inline SINGLE() { }
  inline ~SINGLE() { }
};


template<class Tag, class T>
class Packer { };



/*
	This class Assumes that there exsits a "sizeof" function for
	'class T' basically this is meant for simple objects

	that do not have complex storage mechanism...

	this is a **CHEETA** class with some extra unessesary things taken out
*/

template<class T>
class Packer<SINGLE,T>
{
public:
  //**CHEETA** Return the storage needed to pack the item of type T.
  // For the default case, this is just sizeof(T), but perhaps rounded
  // up to be pointer-word-size aligned.

  static inline int size(const T &)
  {

    return sizeof(double) * ((sizeof(T) + sizeof(double) - 1) / sizeof(double));
  }

  // Pack an item of type T into the given buffer.  Return space used.
  // By default, this just does a placement-new into the buffer,
  // assuming the storage required is sizeof(T).

  static inline int pack(const T &item, char *buffer)
  {
    new ((void*)buffer) T(item);
    return size(item);
  }

  // Unpack an item of type T from the given buffer.  Set the given
  // pointer to point at this item.  Return bytes unpacked.
  // By default, this just recasts the current buffer pointer.

  static inline int unpack(T* &p, char *buffer)
  {
    p = reinterpret_cast<T *>(buffer);
    return size(*p);
  }

  // Delete the item pointed to by the given pointer, that was
  // unpacked with a previous call to unpack().
  // By default, this just runs the destructor on the data, which for
  // many things will do nothing.

  static inline void cleanup(T *p)
  {
    p->~T();
  }
};


//----------------------------------------------------------------------
// ARRAY serialize specialization
//----------------------------------------------------------------------

struct ARRAY
{
  inline ARRAY() { }
  inline ~ARRAY() { }
};


// A specialization arraies of objects

template<class T>
class Packer< ARRAY, T>
{
public:

  // Return the storage needed to pack count items of type T,
  // This includes the bytes needed to store the size of the array.

  static inline int size(const T* items, const int& count)
  {
    int arraySize = count*sizeof(T);
    return ( Packer<SINGLE, int>::size(count)
            + arraySize + padding(arraySize) );
  }

  // Pack an item of type T into the given buffer.  Return space used.
  // By default, this just does a placement-new into the buffer,
  // assuming the storage required is sizeof(T).

  static inline int pack(const T* items, char* buffer, const int& count)
  {
     int n = Packer<SINGLE, int>::pack(count, buffer);
     memcpy(n+buffer, items, count*sizeof(T));
     return size(items, count);
  }

  // Unpack an item of type T from the given buffer.  Set the given
  // pointer to point at this item.  Return bytes unpacked.

  static inline int unpack(T* &p, char *buffer, int& count)
  {
     int* iPtr;
     int n = Packer<SINGLE, int>::unpack(iPtr, buffer);
     count = *iPtr;
     p = reinterpret_cast<T *>(n+buffer);
     return size(p, count);
  }

  // Delete the item pointed to by the given pointer, that was unpacked with a
  //  previous call to unpack(). By default, this just runs the destructor on
  // the data, which for many things will do nothing.  Memory has been
  // allocated from the provided buffer so no freeing of memory need be done
  // here.

  static inline void cleanup(T *p)
  {
    p->~T();
  }
};


//----------------------------------------------------------------------
// **CHEETA** DELEGATE serialize specialization
//
// A specialization tag which provides for delegation of serialization
// responsibilities to the item being packed.  Constructs objects
// from the heap rather than in place from the provided buffer.
//----------------------------------------------------------------------

struct DELEGATE
{
  inline  DELEGATE() { }
  inline ~DELEGATE() { }
};

//
// Delegation is used when DelegateType::delegate == true (the default).
// Otherwise, CHEETAH serialization is used for fundamental types that
// have been specialized with delegate == false (see below).
//

template<class T> class DelegateType
{
  public: enum { delegate = true };
};


//
// This class is used so that serialization routines can be specialized
// for either delegation (WrappedBool<true>) or CHEETAH
// (WrappedBool<false>).
//

template<bool flag> class WrappedBool
{
public:
  WrappedBool()  {}
  ~WrappedBool() {}
};


template<class T>
class Packer< DELEGATE, T>
{
public:

  static inline int size(const T& item)
  {
    return size(item, WrappedBool<DelegateType<T>::delegate>());
  }

  static inline int pack(const T& item, char* buffer)
  {
    return pack(item, buffer, WrappedBool<DelegateType<T>::delegate>());
  }

  static inline int unpack(T* &p, char* buffer)
  {
    return unpack(p, buffer, WrappedBool<DelegateType<T>::delegate>());
  }

  static inline void cleanup(T* p)
  {
    cleanup(p, WrappedBool<DelegateType<T>::delegate>());
  }

private:

  //
  // These private methods, templated on a boolean flag, are
  // used so that fundamental types can be used with DELEGATE
  // serialization.  DelegateType<> must be specialized for
  // all types for which delegation methods are not provided.
  //

  static inline int size(const T& item, const WrappedBool<true>& f)
  {
    //
    // Delegation specialization doesn't construct objects in place
    // (in the provided buffer), so size shouldn't have to worry about word
    // alignment.  However, DELEGATE may be used in conjunction with
    // CHEETAH serialization, so word alignment must be conserved.
    //
    int itemSize = item.size();
    return itemSize + padding(itemSize);
  }

  static inline int size(const T& item, const WrappedBool<false>& f)
  {
    return Packer<SINGLE, T>::size(item);
  }

  static inline int pack(const T& item, char* buffer, const WrappedBool<true>& f)
  {
    int itemSize = item.pack(buffer);
    return itemSize + padding(itemSize);
  }

  static inline int pack(const T& item, char* buffer, const WrappedBool<false>& f)
  {
    return Packer<SINGLE, T>::pack(item, buffer);
  }

  static int unpack(T* &p, char* buffer, const WrappedBool<true>& f)
  {
    p = new T();
    int itemSize = p->unpack(buffer);
    return itemSize + padding(itemSize);
  }

  static int unpack(T* &p, char* buffer, const WrappedBool<false>& f)
  {
    return Packer<SINGLE, T>::unpack(p, buffer);
  }

  static void cleanup(T* p, const WrappedBool<true>& f)
  {
    p->cleanup();
    delete p;
  }

  static void cleanup(T* p, const WrappedBool<false>& f)
  {
    Packer<SINGLE, T>::cleanup(p);
  }
};


//
// DelegateType specializations for common types
//

template<> class DelegateType<bool>	      { public: enum { delegate = false }; };
template<> class DelegateType<char>	      { public: enum { delegate = false }; };
template<> class DelegateType<unsigned char>  { public: enum { delegate = false }; };
template<> class DelegateType<short>	      { public: enum { delegate = false }; };
template<> class DelegateType<unsigned short> { public: enum { delegate = false }; };
template<> class DelegateType<int>	      { public: enum { delegate = false }; };
template<> class DelegateType<unsigned int>   { public: enum { delegate = false }; };
template<> class DelegateType<long>	      { public: enum { delegate = false }; };
template<> class DelegateType<unsigned long>  { public: enum { delegate = false }; };
template<> class DelegateType<float>	      { public: enum { delegate = false }; };
template<> class DelegateType<double>	      { public: enum { delegate = false }; };


//
// DELEGATE specializations for STL strings
//

template<>
class Packer< DELEGATE, std::string>
{
public:

  static inline int size(const std::string& str)
  {
    return Packer<ARRAY, char>::size(0, str.length());
  }

  static int pack(const std::string &str, char* buffer)
  {
    return Packer<ARRAY, char>::pack(str.data(), buffer, str.length());
  }

  static int unpack(std::string &str, char* buffer)
  {
    char* ptr;
    int size;

    int n = Packer<ARRAY, char>::unpack(ptr, buffer, size);
    str = std::string(ptr, size);
    return n;
  }

  static void cleanup(std::string* str) { delete str; }

};

//
// DELEGATE specializations for complex values
//

template<>
class Packer< DELEGATE, Complex<float> >
{
public:

  static inline int size(const Complex<float>& str)
  {
    return Packer<ARRAY, float>::size(0, 2);
  }

  static int pack(const Complex<float> &str, char* buffer)
  {
    float moo[2]={str.Re(), str.Im()};
  //  buffer=new char[size(str)];
    return Packer<ARRAY, float>::pack(moo, buffer, 2);
  }

  static int unpack(Complex<float> &str, char* buffer)
  {
    float* ptr;
    int siz;
    int n = Packer<ARRAY, float>::unpack(ptr, buffer, siz);
   	str = Complex<float>(ptr[0], ptr[1]);
    return n;
  }

  static void cleanup(Complex<float>* str) { delete str; }

};

//
// DELEGATE specializations for complex values
//

template<>
class Packer< DELEGATE, Complex<double> >
{
public:

  static inline int size(const Complex<double>& str)
  {
    return Packer<ARRAY, double>::size(0, 2);
  }

  static int pack(const Complex<double> &str, char* buffer)
  {
    double moo[2]={str.Re(), str.Im()};
  //  buffer=new char[size(str)];
    return Packer<ARRAY, double>::pack(moo, buffer, 2);
  }

  static int unpack(Complex<double> &str, char* buffer)
  {
    double* ptr;
    int siz;
    int n = Packer<ARRAY, double>::unpack(ptr, buffer, siz);
   	str = Complex<double>(ptr[0], ptr[1]);
    return n;
  }

  static void cleanup(Complex<double>* str) { delete str; }

};

//
// DELEGATE specializations for coord values
//

template<class T, int N>
class Packer< DELEGATE, coord<T,N> >
{
public:

  static inline int size(const coord<T,N>& str)
  {
    return Packer<ARRAY, T>::size(0, N);
  }

  static int pack(const coord<T,N> &str, char* buffer)
  {
    T *moo=const_cast<T *>(str.data());
    return Packer<ARRAY, T>::pack(moo, buffer, N);
  }

  static int unpack(coord<T,N> &str, char* buffer)
  {
    T* ptr;
    int siz;
    int n = Packer<ARRAY, T>::unpack(ptr, buffer, siz);
   	str = coord<T,N>(ptr);
    return n;
  }

  static void cleanup(coord<T,N>* str) { delete str; }

};

//
// DELEGATE specializations for Vector arraies
//

template<class T>
class Packer< DELEGATE, Vector<T> >
{
public:

  static inline int size(const Vector<T>& str)
  {
    return Packer<ARRAY, T>::size(0, str.size());
  }

  static int pack(const Vector<T> &str, char* buffer)
  {
    return Packer<ARRAY, T>::pack(str.data(), buffer, str.size());
  }

  static int unpack(Vector<T> &str, char* buffer)
  {
    T* ptr;
    int siz;
    int n = Packer<ARRAY, T>::unpack(ptr, buffer, siz);
   	str = Vector<T>(ptr, siz);
    return n;
  }

  static void cleanup(Vector<T>* str) { delete str; }

};

//
// DELEGATE specializations for matrix
//

template<class T, class Struct_t>
class Packer< DELEGATE, _matrix<T, Struct_t> >
{
public:

  static inline int size(const _matrix<T, Struct_t>& str)
  {
	return Packer<SINGLE, int>::size(str.rows())+
		Packer<SINGLE, int>::size(str.cols())+
		Packer<ARRAY, T>::size(0, str.numElements());
  }

  static int pack(const _matrix<T, Struct_t> &str, char* buffer)
  {
   		int ss = Packer<SINGLE, int>::pack(str.rows(), buffer);
    	ss += Packer<SINGLE, int>::pack(str.cols(), buffer+ss);
    	return 	ss+Packer<ARRAY, T>::pack(str.data(), buffer+ss, str.numElements());
  }

  static int unpack( _matrix<T, Struct_t>  &str, char* buffer)
  {
    T* ptr;
    int siz,*r,*c,n;
    n = Packer<SINGLE, int>::unpack(r, buffer);
    n += Packer<SINGLE, int>::unpack(c, buffer+n);
    n += Packer<ARRAY, T>::unpack(ptr, buffer+n, siz);
    str =  _matrix<T, Struct_t>(*r,*c,ptr);
    return n;
  }

  static void cleanup(_matrix<T, Struct_t> * str) { delete str; }

};

//
// DELEGATE specializations for identity matrix(a weird 'one element' array)
// here we do not pass anything except the rows and columns size

template<class T>
class Packer< DELEGATE, _matrix<T, IdentityMatrix> >
{
public:

  static inline int size(const _matrix<T, IdentityMatrix>& str)
  {
	return Packer<SINGLE, int>::size(str.rows())+
		Packer<SINGLE, int>::size(str.cols());
  }

  static int pack(const _matrix<T, IdentityMatrix> &str, char* buffer)
  {
   		int ss = Packer<SINGLE, int>::pack(str.rows(), buffer);
    	return ss+Packer<SINGLE, int>::pack(str.cols(), buffer+ss);

  }

  static int unpack( _matrix<T, IdentityMatrix>  &str, char* buffer)
  {
    int *r,*c,n;
    n = Packer<SINGLE, int>::unpack(r, buffer);
    n += Packer<SINGLE, int>::unpack(c, buffer+n);
    str =  _matrix<T, IdentityMatrix>(*r,*c);
    return n;
  }

  static void cleanup(_matrix<T, IdentityMatrix> * str) { delete str; }

};

END_BL_NAMESPACE
#endif //HAVE MPI

#endif


