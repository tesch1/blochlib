/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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

/****************************************************************************
 * Vector.h
 	A really nice expression templated class for vectors
 	this is probably as fast as C++ can get for a general vector of objects
 	uses template expressions to expand all operations, removes all
 	temporary objects, and excutes ONE loop, for each big long
 	line of operations

 	based roughly on the excelent teachings of blitz++
 	(http://oonumerics.org/blitz)

 	if you never seen expression templates, this code will
 	look like a pile of raging letters with little meaning
 	the coding is not intuitive.
 */

//#define VBoundCheck 1

#ifndef _Vector_h_
#define _Vector_h_ 1

#ifndef ON_WINDOWS
#include "blochconfig.h"
#endif

//#include <math.h>
#include "container/rankType.h"
#include <string>
#include "container/complex.h"
#include "container/MemChunk.h"
#include<iostream>
#include<iomanip>
#include "container/operations.h"  		//contains the applicative template operations like Add, Mul, Div, etc
#include "container/range.h"


BEGIN_BL_NAMESPACE


template<class T>
class Vector;

// Vector
template<class NumType_t>
struct ObjectTrait<Vector<NumType_t> > {
    typedef NumType_t numtype;
};

template<class T, int N>
class coord;

//class complex;
template<class T>	class VExpr;
template<class T>	class VExprNonConst;
template<class T>	class VIter;
template<class T>	class VIterConst;

template<class expr, class Op>	class VBinUnaryOp;
template<class T>	class VExprConst;
template<class T1, class T2, class Op> class VBinExprOp;

template<class T>
class Vector : public MemChunkRef<T>{

private:
  using MemChunkRef<T>::data_;
    int length_;
    void lenError(int i)const throw()
    {
		if(i<0){
			BLEXCEPTION(" reqeust for index too small (negative) in vector")
		}
		if(i>=length_){
			BLEXCEPTION(" reqeust for index too large in vector")
		}
	}

 	void lenError(int i){
		if(i<0){
			BLEXCEPTION(" reqeust for index too small (negative) in vector")
		}
		if(i>=length_){
			BLEXCEPTION(" reqeust for index too large in vector")
		}
	}

        void RowErr() const {
		BLEXCEPTION(" request for non exsistant element....")
	}



public:
   typedef SumType(T) sumType;
   typedef T numtype;
   typedef VIter<T> iterator;
   typedef VIterConst<T> const_iterator;

   /* this was supposed to act like flag to put at the end of
      the print BUT it takes up ALOT of CPU to carry around and initialize
      SO it was killed
   */
   //std::string RecordSep;
	~Vector() {}

    explicit Vector(int n) : MemChunkRef<T>(n), length_(n){} //,RecordSep(" ") {}

    Vector(): MemChunkRef<T>(0), length_(0){}//,RecordSep(" ") {}

	Vector(T *in, int N) : MemChunkRef<T>(N, in), length_(N){}// ,RecordSep(" ") { }

	Vector(Vector<T> &in) : MemChunkRef<T>(in), length_(in.length()){}// ,RecordSep(in.RecordSep) {}
	Vector(const Vector<T> &in) : MemChunkRef<T>(in), length_(in.length()){} //,RecordSep(in.RecordSep){}

	template<class expr>
	Vector(const VExpr<expr> &in): MemChunkRef<T>(in.guessLen()){//, RecordSep(" "){
		length_=in.guessLen();
		(*this)=in;
	}

	template<class T2>
	Vector(const Vector<T2> &in): MemChunkRef<T>(in.length()){//, RecordSep(" "){
		length_=in.length();
		(*this)=in;
	}

	Vector(int le, const T &in): MemChunkRef<T>(le){//, RecordSep(" "){
		length_=le;
		(*this)=in;
	}

	template<class expr>
	Vector(int le, const VExpr<expr> &in): MemChunkRef<T>(le){//, RecordSep(" "){
		length_=le;
		(*this)=in;
	}

	template<class T1>
	Vector(int le, const T1 &in): MemChunkRef<T>(le){//, RecordSep(" "){
		length_=le;
		(*this)=T(in);
	}

//A Vector from a range
	Vector(Range r)
		: MemChunkRef<T>(r.guessLen())
	{
		length_ = r.guessLen();
		(*this) = VExpr<Range>(r);
    }

//A Vector from a spread
	template<class numt>
	Vector(Spread<numt> r)
		: MemChunkRef<T>(r.guessLen())
	{
		length_ = r.guessLen();
		(*this) = VExpr<Spread<numt> >(r);
    }

//a vector froma range and an old vector
    Vector(Vector<T>& vec, Range r)
        : MemChunkRef<T>(vec, r.first(0))
    {

        if(r.first(0) <0 || r.first(0) > vec.length())
        {
        	std::cerr<<std::endl<<" Error: Vector(Vector, range) "<<std::endl;
        	std::cerr<<" Cannot create vector, first element in Range"<<std::endl;
        	std::cerr<<" is too large (or small) "<<std::endl;
        }
        if(r.last(vec.length()-1) <0 || r.last(vec.length()-1) > vec.length())
        {
        	std::cerr<<std::endl<<" Error: Vector(Vector, range) "<<std::endl;
        	std::cerr<<" Cannot create vector, last element in Range"<<std::endl;
        	std::cerr<<" is too large (or small) "<<std::endl;
        }

        length_ = (r.last(vec.length()-1) - r.first(0)) / r.stride() + 1;
    }

    inline T& operator[](int i) {
#ifdef VBoundCheck
		lenError(i);
#endif
		return data_[i];
	}

	inline T operator[](int i)const {
#ifdef VBoundCheck
		lenError(i);
#endif
		return data_[i];
	}

	inline T operator()(int i)const {
#ifdef VBoundCheck
		lenError(i);
#endif
		return data_[i];
	}

	inline T &operator()(int i) {
#ifdef VBoundCheck
		lenError(i);
#endif
		return data_[i];
	}

	Vector      operator()(Range r)
	{
		return Vector(*this, r);
	}

	Vector      operator[](Range r)
	{
		return Vector(*this, r);
	}


    const_iterator const_iter() const { return const_iterator(*this); }

    int begin() const {	return 0;	}
    int end()	const	{	return length_;	}

    int begin(int glen) const {	return 0;	}
    int end(int gend)	const	{	return length_;	}

    iterator iter() { return iterator(*this);	}

    int length() const { return length_; }
    int length(int slen) const {	return length_;	}
    int guessLen() const { return length_;	}
	bool empty()const { return length_>0?false:true;	}
    int size() const { return length_; }

    const T* data()const{ return data_;	}
    T* data(){ return data_;	}

	template<class T2> inline Vector<T>& operator+=(const Vector<T2> &);
    template<class T2> inline Vector<T>& operator-=(const Vector<T2> &);
    template<class T2> inline Vector<T>& operator*=(const Vector<T2> &);
    template<class T2> inline Vector<T>& operator/=(const Vector<T2> &);

	template<class expr> inline Vector<T>& operator+=(const VExpr<expr> &);
    template<class expr> inline Vector<T>& operator-=(const VExpr<expr> &);
    template<class expr> inline Vector<T>& operator/=(const VExpr<expr> &);
    template<class expr> inline Vector<T>& operator*=(const VExpr<expr> &);

 	template<class T1,int N>
 	inline Vector<T>& operator+=(const coord<T1,N> &);
    template<class T1,int N>
    inline Vector<T>& operator-=(const coord<T1,N> &);
    template<class T1,int N>
    inline Vector<T>& operator/=(const coord<T1,N> &);
    template<class T1,int N>
    inline Vector<T>& operator*=(const coord<T1,N> &);


 	inline Vector<T>& operator+=(const complex &);
    inline Vector<T>& operator-=(const complex &);
    inline Vector<T>& operator*=(const complex &);
    inline Vector<T>& operator/=(const complex &);

    inline Vector<T>& operator+=(const double &);
    inline Vector<T>& operator-=(const double &);
    inline Vector<T>& operator*=(const double &);
    inline Vector<T>& operator/=(const double &);

	inline Vector<T>& operator+=(const float &);
    inline Vector<T>& operator-=(const float &);
    inline Vector<T>& operator*=(const float &);
    inline Vector<T>& operator/=(const float &);

    inline Vector<T>& operator+=(const int &);
    inline Vector<T>& operator-=(const int &);
    inline Vector<T>& operator*=(const int &);
    inline Vector<T>& operator/=(const int &);

    inline Vector<T>& operator+=(const short &);
    inline Vector<T>& operator-=(const short &);
    inline Vector<T>& operator*=(const short &);
    inline Vector<T>& operator/=(const short &);

    inline Vector<T>& operator+=(const long &);
    inline Vector<T>& operator-=(const long &);
    inline Vector<T>& operator*=(const long &);
    inline Vector<T>& operator/=(const long &);

    inline Vector<T>& operator+=(const char &);
    inline Vector<T>& operator-=(const char &);
    inline Vector<T>& operator*=(const char &);
    inline Vector<T>& operator/=(const char &);

    inline Vector<T>& operator+=(const bool &);
    inline Vector<T>& operator-=(const bool &);
    inline Vector<T>& operator*=(const bool &);
    inline Vector<T>& operator/=(const bool &);

    template<class T1>
    Vector& operator=(const VExpr<T1>& x){
   		resize(x.length(length_));

#ifdef VBoundCheck
		if(x.length(length_) != length_){
			BLEXCEPTION(" vectors must be same length in Expression (i.e. s '=' 3*w) assignment")
		}
#endif
		for(int i=x.begin(0);i<x.end(length_);++i){
			data_[i]=(T)x(i);
		};
		return *this;
	}

	template<class T1>
    Vector &operator=(VExpr<T1>& x){
    	resize(x.length(length_));

#ifdef VBoundCheck
		if(x.length(length_) != length_){
			BLEXCEPTION(" vectors must be same length in Expression (i.e. s '=' 3*w) assignment")
		}
#endif

		for(int i=x.begin(0);i<x.end(length_);++i){
			data_[i]=(T)x(i);
		};
		return *this;
	}
    /*template<class T1>
    Vector& operator=(const Vector<T1>& x){
    	(*this)=VExpr<Vector<T1> >(x);
    	return *this;
    }*/


	template<class T1>
	void operator=(const Vector<T1> &x){
		resize(x.length());
		for(int i=0;i<length_;++i){
			data_[i]=T(x[i]);
		}
		//return *this;
	}

	void operator=(const Range &x){
		resize(x.length());
		for(int i=0;i<length_;++i){
			data_[i]=T(x(i));
		}
		//return *this;
	}


	void operator=(const Vector &x){
		if(this==&x) return;
		resize(x.length());
		for(int i=0;i<length_;++i){
			data_[i]=x[i];
		}
	}

    template<class T1>
    void operator=(const T1& x){
    	for(int i=0;i<length_;++i){
    		data_[i]=T(x);
    	}
    	//return *this;
    }

    void operator=(const T& x){
    	for(int i=0;i<length_;++i){
    		data_[i]=x;
    	}
    //	return *this;
    }

	Vector copy() const{
	   Vector<T> newCopy(length_);
	   memcpy(newCopy.data(), data(), length_ * sizeof(T));
	   return newCopy;
	}

	template<class T1>
	void copy(const Vector<T1> &in){
		*this=in;
	}

	void reference(Vector<T>& x){
	    MemChunkRef<T>::changeBlock(x, 0);
	    length_ = x.length_;
	}

	void resize(int length){
	    if (length != length_){
	        MemChunkRef<T>::newBlock(length);
	        length_ = length;
	    }
	}

	void resize(int length, T in){
		if (length != length_){
			MemChunkRef<T>::newBlock(length);
			length_ = length;
		}
		fill(in);
	}

	void fill(T in){	(*this)=in; 	}

//aplies a 'function' across an entire vector
//The input class MUST have the operator(T, int) defined

	template<class Func>
	void apply( Func &f)
	{
		for(int i=0;i<length();++i)
		{
			(*this)(i)=f((*this)(i), i);
		}
	}

//aplies a 'function' across an a Range in the vector
//The input class MUST have the operator(T, int) defined

	template<class Func>
	void apply( Func &f, const Range &r)
	{
		for(int i=0;i<r.length(length());++i)
		{
			(*this)(r(i))=f(get(r(i)), r(i));
		}
	}

	T get(int i){	return (*this)[i];	}

	template<class T1>
	void put(int i, T1 in){ (*this)[i]=T(in);	}

	void put(int i, T in){ (*this)[i]=in;	}

//put functions for Ranges
	template<class T1>
	void put(const Range &r,const Vector<T1> &in)
	{
#ifdef VBoundCheck
		if(r.length(length())>length() || r.length(length())>in.length()) RowErr();
#endif
		for(int i=0;i<r.length(length());i++)
		{
			put(r(i),T(in(i)));
		}
	}

//put functions for Ranges
	template<class expr>
	void put(const Range &r,const VExpr<expr> &in)
	{
#ifdef VBoundCheck
		if(r.length(length())>length() || r.length(length())>in.length()) RowErr();
#endif
		for(int i=0;i<r.length(length());i++)
		{
			put(r(i),T(in(i)));
		}
	}



	void resizeAndPreserve(int newLength){
	    Vector<T> newVector(newLength);
	    if (newLength > length_){
				for(int i=0;i<length_;++i){		newVector[i]=(*this)[i];	}
		}else{
	   		for(int i=0;i<newLength;++i){	newVector[i]=(*this)[i];	}
		}
		reference(newVector);
		// *this=newVector;

	}

	template<class nt1>
	void resizeAndPreserve(int newLength, const nt1 &in){
	    Vector<T> newVector(newLength);
	    if (newLength > length_)
				for(int i=0;i<length_;++i)		newVector[i]=(*this)[i];
	    else
	   		for(int i=0;i<newLength;++i)	newVector[i]=(*this)[i];

			for(int i=length_-1; i<newLength;i++) newVector[i]=T(in);
	    reference(newVector);
	}

	template<class nt1>
	void push_back(nt1 in){
		//MemChunkRef<T>::push_back(in);
			Vector<T> newVec(length_+1);
			for(int i=0;i<length_;++i)		newVec[i]=(*this)[i];
			newVec.put(length_,T(in));
			reference(newVec);
	}

	void push_back(T in){
		//MemChunkRef<T>::push_back(in);
			Vector<T> newVec(length_+1);
			for(int i=0;i<length_;i++)		newVec[i]=(*this)[i];
			newVec.put(length_,in);
			reference(newVec);
		//	*this=newVec;
	}

	void pop(){
		if(length_>0) resizeAndPreserve(length_-1);
	}

#if 0
	void drop(int i){
#ifdef VBoundCheck
		if(i>=length()) RowErr();
#endif
		if(length_>0){
			Vector<T> cp(size()-1);
			for(int k=0;k<i;++k)	cp[k]=get(k);
			for(int j=i+1, k=i;j<size();++k, ++j) cp[k]=get(j);
			*this=cp;
		}
	}
#endif

	void print(std::ostream &);
	template<class NextThing>
	void print(std::ostream &, NextThing out);


};


//*********************************************
//		itterator for the vector
//*********************************************
template<class type>
class VIter {
	private:
		VIter() { }
		 type* data_;
		int length_;
		int cur_;
		int stop_;
		void lenError(int i)const{
			if(i<0){
				BLEXCEPTION(" reqeust for index too small (negative) in vector")
			}
			if(i>=length_){
				BLEXCEPTION(" reqeust for index too large in vector")
			}
		}
	public:
		typedef type numtype;

		explicit VIter(Vector<type>& x)  : data_(x.data()){
			cur_=0;
			length_ = x.length();
			stop_=length_;
		}


		VIter(type* data,int length)    :
			data_(data), length_(length), stop_(length)
		{ }

		VIter(const VIter<type>& x) {
			cur_=0;
			data_ = x.data_;
			length_ = x.length_;
			stop_=length_;
		}

		~VIter(){	data_=NULL;	}

		type * data(){	return data_;	}

		type operator[](int i) const{
#ifdef VBoundCheck
			lenError(i);
#endif
			return data_[i];
		}

		type&  operator[](int i){
#ifdef VBoundCheck
			lenError(i);
#endif
			return data_[i];
		}

		type operator()(int i) const{
#ifdef VBoundCheck
			lenError(i);
#endif
			return data_[i];
		}

		type& operator()(int i) {
#ifdef VBoundCheck
			lenError(i);
#endif
			return data_[i];
		}

		type shift(int i) const
		{
#ifdef VBoundCheck
			lenError(cur_+i);
#endif
			return data_[cur_+i];
		}

		type &shift(int i)
		{
#ifdef VBoundCheck
			lenError(cur_+i);
#endif
			return data_[cur_+i];
		}

		type &operator()(void){	return (data_[cur_]);	}
		type operator()(void) const {	return data_[cur_];	}

		void moveTo(int i) const
		{
#ifdef VBoundCheck
			lenError(i);
#endif
			int &ref=const_cast<int &>( cur_);
			ref=i;
		}

		void begin(int i)
		{
#ifdef VBoundCheck
			lenError(i);
#endif
			cur_=i;
		}

		void end(int i)
		{
#ifdef VBoundCheck
			lenError(i);
#endif
			stop_=i;
		}



		inline int end() const {return length_; }
		inline int begin() const { 	return 0; }

		void reset(){	cur_=0; stop_=length_;	}
		inline int size()	const {	return length_;	}

		void operator++()  {
			if(cur_<stop_){ ++cur_;	}
		}

		void operator--()  {
			if(cur_>0){ --cur_;	}
		}

		operator bool() const {		return cur_<(stop_);		}

		inline int curpos() const {	return cur_;	}
		inline int curpos() {	return cur_;	}

		type operator*() const{ return *data_; }
		type& operator*(){ return *data_; }

		int length(int) const { return length_; }
		int guessLen() const { return length_;	}


};

//*********************************************
//		Constant itterator for the vector
//*********************************************

template<class type>
class VIterConst {
	private:
		const type * data_;
		int length_;
		int cur_;
		 void lenError(int i) const{
			if(i<0){
				BLEXCEPTION(" reqeust for index too small (negative) in vector")
			}
			if(i>=length_){
				BLEXCEPTION(" reqeust for index too large in vector")
			}
		}
	public:
		typedef type numtype;

		explicit VIterConst(const Vector<type>& x)  : data_(x.data()) , cur_(0){
			length_ = x.length();
		}

		VIterConst(const VIterConst<type>& x){
			data_ = x.data_;
			cur_=0;
			length_ = x.length_;
		}

		type operator[](int i) const{
#ifdef VBoundCheck
			lenError(i);
#endif
			return data_[i];
		}

		type operator()(int i) const{
#ifdef VBoundCheck
			lenError(i);
#endif
			return data_[i];
		}

		inline int end() const {return length_; }
		inline int begin() const { 	return 0; }

		inline int end(int ge) const {return length_; }
		inline int begin(int gb) const { 	return 0; }

		void operator++() {
			if(cur_<length_) cur_++;
		}
		operator bool() const {
			return cur_>(length-1);
		}

		int length(int slen) const  { return length_; }
		int guessLen() const { return length_;	}

};



//******************************************************************
// * VExpr -- a vector expression iterator
//******************************************************************


template<class expr>
class VExpr {

private:
    expr iter_;
    VExpr() {}

public:
    typedef typename expr::numtype numtype;

    VExpr(const expr& a)  : iter_(a) { }

    int length(int slen) const {	return iter_.length(slen); }
    int begin(int blen) const	{	return iter_.begin(blen);	}
    int end(int elen) const	{	return iter_.end(elen);	}

    int guessLen() const { return iter_.guessLen();	}
    int size() const {	return iter_.size();	}

   	numtype operator()(int i) const{ return iter_(i);	}

    numtype operator[](int i) const  { return iter_[i]; }

    void operator++()
    { ++iter_; }
};

template<class expr>
class VExprNonConst {

private:
    expr iter_;
    VExprNonConst() {}

public:
    typedef typename expr::numtype numtype;

    VExprNonConst(const expr& a)  : iter_(a) { }

    int length(int slen) const {	return iter_.length(slen); }
    int begin(int blen) const	{	return iter_.begin(blen);	}
    int end(int elen) const	{	return iter_.end(elen);	}

    int guessLen() const { return iter_.guessLen();	}
    int size() const {	return iter_.size();	}

   	numtype operator()(int i) { return iter_(i);	}

    numtype operator[](int i)   { return iter_[i]; }

    void operator++()
    { ++iter_; }
};



/********************************
// for constant number operations
*********************************/
template<class T>
class VExprConst{

private:
    T iter_;
    VExprConst(){}

public:
	typedef T numtype;

    VExprConst(const numtype& a)    : iter_(a) { }

	int length(int slen)const{ return slen;	} //keeps the propogation of lengths 'good'
	int guessLen() const { return 0;	}

 	int begin(int blen) const	{	return blen;	}
    int end(int elen) const	{	return elen;	}

    numtype operator() (int i) const
    { return iter_; }

    numtype operator[](int i) const
    { return iter_; }

    void operator++() { }
};

/****************************************************************************
 * the major work horse...the guy at the end of the entire 'unrollin'
 * of the expression (a+b*abs(c)...) that 'applies' the operation
 * to each element to make it look like one loop
 * ie a[i]+b[i]*abs(c[i])...the acctual loop is done in the Vector class
 * upon assignment.. Vector<T> = a+b*abs(c)
 */
template<class T1, class T2, class Op>
class VBinExprOp {

private:
    T1 iter1_;
    T2 iter2_;
    VBinExprOp(){}

public:

    VBinExprOp(const T1& a, const T2& b) : iter1_(a), iter2_(b)  { }

	typedef typename T1::numtype type1;
	typedef typename T2::numtype type2;
	typedef OutType(type1,type2) numtype;

	int length(int slen) const {

#ifdef VBoundCheck
		if(iter1_.length(slen) != iter2_.length(slen)){
			BLEXCEPTION(" lengths must be the same to perform operation")
		}
#endif

		return iter1_.length(slen);
	}
	int guessLen() const { return (iter1_.guessLen()>iter2_.guessLen())?iter1_.guessLen():iter2_.guessLen();	}

 	int begin(int blen) const	{	return (iter1_.begin(blen)>iter2_.begin(blen))?iter1_.begin(blen):iter2_.begin(blen) ;	}
    int end(int elen) const	{	return (iter1_.end(elen)>iter2_.end(elen))?iter1_.end(elen):iter2_.end(elen);	}

    void operator++()
    { ++iter1_; ++iter2_; }

    numtype operator()(int i) const
    { return Op::apply(iter1_(i),iter2_(i)); }

    numtype operator[](int i) const
    { return Op::apply(iter1_[i], iter2_[i]); }
};

//******************************************************************
// * VExprUnary -- a vector expression iterator for ONE op
//	like abs(), cos(), etc
//******************************************************************


template<class expr, class Op>
class VBinUnaryOp{
	private:
		VBinUnaryOp() { }
		expr iter_;
	public:
		typedef typename Op::numtype numtype;

		VBinUnaryOp(const expr& iter): iter_(iter) { }

		VBinUnaryOp(const VBinUnaryOp<expr, Op>& x) : iter_(x.iter_) { }

 		int begin(int blen) const	{	return iter_.begin(blen);	}
    	int end(int elen) const	{	return iter_.end(elen);	}


		numtype operator[](int i) const
		{ return Op::apply(iter_[i]); }

		numtype operator()(int i) const
		{ return Op::apply(iter_(i)); }

		int length(int slen) const	{ return iter_.length(slen); }
		int guessLen() const { return iter_.guessLen();	}


};


#define VecFromVec(op)	\
	template<class T> template<class T1>\
	inline Vector<T> & Vector<T>::operator op (const Vector<T1> &a){	\
		resize(a.length());	\
		for (int i=0; i < length_; ++i) {                              \
			(*this)[i] op a[i];                                      \
		}    								\
		return *this;	\
	}	\

VecFromVec(+=)
VecFromVec(-=)
VecFromVec(/=)
VecFromVec(*=)



#define VecFromConst(type, op)	\
	template<class T>	\
	inline Vector<T> & Vector<T>::operator op (const type &a){	\
		(*this) op VExpr<VExprConst<T> >(VExprConst<T>(a));		\
		return *this;	\
	} \

VecFromConst(complex,+=)
VecFromConst(complex, -=)
VecFromConst(complex, /=)
VecFromConst(complex,*=)

VecFromConst(double,+=)
VecFromConst(double, -=)
VecFromConst(double, /=)
VecFromConst(double,*=)

VecFromConst(float,+=)
VecFromConst(float, -=)
VecFromConst(float, /=)
VecFromConst(float,*=)

VecFromConst(int,+=)
VecFromConst(int, -=)
VecFromConst(int, /=)
VecFromConst(int,*=)

VecFromConst(short,+=)
VecFromConst(short, -=)
VecFromConst(short, /=)
VecFromConst(short,*=)

VecFromConst(long,+=)
VecFromConst(long, -=)
VecFromConst(long, /=)
VecFromConst(long,*=)

VecFromConst(char,+=)
VecFromConst(char, -=)
VecFromConst(char, /=)
VecFromConst(char,*=)

VecFromConst(bool,+=)
VecFromConst(bool, -=)
VecFromConst(bool, /=)
VecFromConst(bool,*=)



//****************************************************************************
//		allow for a Vector<T1> from VExp<T1>
//		 for use with unaray ops i.e (+=, *=)
//****************************************************************************
#define VecFromExpr(op)                                            \
	template<class T> template<class expr>                        \
	inline Vector<T>& Vector<T>::operator op (const VExpr<expr> &Texpr){     \
		for (int i=0; i < length_; ++i) {                              \
			(*this)[i] op Texpr[i];                                      \
		}                                                                   \
    	return *this;				\
	}

VecFromExpr(+=)
VecFromExpr(-=)
VecFromExpr(*=)
VecFromExpr(/=)


//****************************************************************************


//*************************************************************************************************
//**************OPERATION"S IN A CAN***************************************************************
//*************************************************************************************************

//
// VECTOR BINARY OPERTIONS OPERATORS
//  ordered like sooo
// **********Vec OP number*********
// **********number OP Vec*********
// *************Expr OP number**************
// *************number OP Expr**************
// *****************Vec OP Vec****************
// *************Expr OP Vec**************
// *************Vec OP Expr**************
// *************Expr OP Expr**************


//macro'ed for our typing pleasure....

#define VecBinOps(FNAME, OP) \
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VIterConst<T1>,VIterConst<T2>,FNAME<T1, T2> > >   operator OP (const Vector<T1>& a, const Vector<T2>& b)		\
{		\
    typedef VBinExprOp<VIterConst<T1>,VIterConst<T2>,FNAME<T1,T2> > ExprT;		\
    return VExpr<ExprT>(ExprT(a.const_iter(),b.const_iter()));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VExpr<T1>,VIterConst<T2>,FNAME<typename T1::numtype, T2> > >  operator OP (const VExpr<T1>& a, const Vector<T2>& b)		\
{		\
    typedef VBinExprOp<VExpr<T1>,VIterConst<T2>,FNAME<typename T1::numtype, T2> > ExprT;		\
    return VExpr<ExprT>(ExprT(a,b.const_iter()));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VIterConst<T1>,VExpr<T2>,FNAME<T1, typename T2::numtype> > >  operator OP (const Vector<T1>& a, const VExpr<T2>& b)		\
{		\
    typedef VBinExprOp<VIterConst<T1>,VExpr<T2>,FNAME<T1, typename T2::numtype> > ExprT;		\
    return VExpr<ExprT>(ExprT(a.const_iter(),b));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VExpr<T1>,VExpr<T2>,FNAME<typename T1::numtype, typename T2::numtype> > >		\
	operator OP (const VExpr<T1>& a, const VExpr<T2>& b)		\
{		\
    typedef VBinExprOp<VExpr<T1>,VExpr<T2>,FNAME<typename T1::numtype, typename T2::numtype> > ExprT;		\
    return VExpr<ExprT>(ExprT(a,b));		\
}		\


VecBinOps(ApAdd, +)
VecBinOps(ApSubtract, -)
VecBinOps(ApMultiply, *)
VecBinOps(ApDivide, /)

#define VecBinOpsConst(TYPE, FNAME, OP) \
\
template<class T1>		\
inline VExpr<VBinExprOp<VIterConst<T1>,VExprConst<TYPE>,FNAME<T1, TYPE > > >   operator OP (const Vector<T1>& a, const TYPE& b)		\
{		\
    typedef VBinExprOp<VIterConst<T1>,VExprConst<TYPE>,FNAME<T1,TYPE > > ExprT;		\
    return VExpr<ExprT>(ExprT(a.const_iter(),VExprConst<TYPE>(b)));		\
}		\
		\
template<class T2>		\
inline VExpr<VBinExprOp<VExprConst<TYPE>,VIterConst<T2>,FNAME<TYPE, T2> > >   operator OP (const TYPE& a, const Vector<T2>& b)		\
{	\
    typedef VBinExprOp<VExprConst<TYPE>,VIterConst<T2>,FNAME<TYPE,T2> > ExprT;		\
    return VExpr<ExprT>(ExprT(VExprConst<TYPE>(a),b.const_iter()));		\
}		\
		\
template<class T1>		\
inline VExpr<VBinExprOp<VExpr<T1>,VExprConst<TYPE>,FNAME<typename T1::numtype, TYPE> > >  operator OP (const VExpr<T1>& a, const TYPE& b)		\
{		\
    typedef VBinExprOp<VExpr<T1>,VExprConst<TYPE>,FNAME<typename T1::numtype, TYPE> > ExprT;		\
    return VExpr<ExprT>(ExprT(a,VExprConst<TYPE>(b)));		\
}		\
		\
template<class T2>		\
inline VExpr<VBinExprOp<VExprConst<TYPE>,VExpr<T2>,FNAME<TYPE,typename T2::numtype> > >  operator OP (const TYPE& a,const VExpr<T2>& b)		\
{		\
    typedef VBinExprOp<VExprConst<TYPE>,VExpr<T2>,FNAME<TYPE,typename T2::numtype> > ExprT;		\
    return VExpr<ExprT>(ExprT(VExprConst<TYPE>(a),b));		\
}		\


VecBinOpsConst(complex, ApAdd, +)
VecBinOpsConst(complex, ApSubtract, -)
VecBinOpsConst(complex, ApMultiply, *)
VecBinOpsConst(complex, ApDivide, /)

VecBinOpsConst(double, ApAdd, +)
VecBinOpsConst(double, ApSubtract, -)
VecBinOpsConst(double, ApMultiply, *)
VecBinOpsConst(double, ApDivide, /)

VecBinOpsConst(float, ApAdd, +)
VecBinOpsConst(float, ApSubtract, -)
VecBinOpsConst(float, ApMultiply, *)
VecBinOpsConst(float, ApDivide, /)

VecBinOpsConst(int, ApAdd, +)
VecBinOpsConst(int, ApSubtract, -)
VecBinOpsConst(int, ApMultiply, *)
VecBinOpsConst(int, ApDivide, /)

VecBinOpsConst(long, ApAdd, +)
VecBinOpsConst(long, ApSubtract, -)
VecBinOpsConst(long, ApMultiply, *)
VecBinOpsConst(long, ApDivide, /)

VecBinOpsConst(short, ApAdd, +)
VecBinOpsConst(short, ApSubtract, -)
VecBinOpsConst(short, ApMultiply, *)
VecBinOpsConst(short, ApDivide, /)

VecBinOpsConst(char, ApAdd, +)
VecBinOpsConst(char, ApSubtract, -)
VecBinOpsConst(char, ApMultiply, *)
VecBinOpsConst(char, ApDivide, /)

VecBinOpsConst(bool, ApAdd, +)
VecBinOpsConst(bool, ApSubtract, -)
VecBinOpsConst(bool, ApMultiply, *)
VecBinOpsConst(bool, ApDivide, /)

// ****************************************************
// SPECIAL CASE FOR "NON-Operators Binary Operations"
//*****************************************************

#define VecSpecialBinOps(FNAME, OP) \
\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VIterConst<T1>,VExprConst<T2>,FNAME<T1, T2> > >   OP (const Vector<T1>& a, const T2& b)		\
{		\
    typedef VBinExprOp<VIterConst<T1>,VExprConst<T2>,FNAME<T1,T2> > ExprT;		\
    return VExpr<ExprT>(ExprT(a.const_iter(),VExprConst<T2>(b)));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VExprConst<T1>,VIterConst<T2>,FNAME<T1, T2> > >   OP (const T1& a, const Vector<T2>& b)		\
{	\
    typedef VBinExprOp<VExprConst<T1>,VIterConst<T2>,FNAME<T1,T2> > ExprT;		\
    return VExpr<ExprT>(ExprT(VExprConst<T1>(a),b.begin()));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VExpr<T1>,VExprConst<T2>,FNAME<typename T1::numtype, T2> > >  OP (const VExpr<T1>& a, const T2& b)		\
{		\
    typedef VBinExprOp<VExpr<T1>,VExprConst<T2>,FNAME<typename T1::numtype, T2> > ExprT;		\
    return VExpr<ExprT>(ExprT(a,VExprConst<T2>(b)));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VExprConst<T1>,VExpr<T2>,FNAME<T1,typename T2::numtype> > >  OP (const T1& a,const VExpr<T2>& b)		\
{		\
    typedef VBinExprOp<VExprConst<T1>,VExpr<T2>,FNAME<T1,typename T2::numtype> > ExprT;		\
    return VExpr<ExprT>(ExprT(VExprConst<T1>(a),b));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VIterConst<T1>,VIterConst<T2>,FNAME<T1, T2> > >   OP (const Vector<T1>& a, const Vector<T2>& b)		\
{		\
    typedef VBinExprOp<VIterConst<T1>,VIterConst<T2>,FNAME<T1,T2> > ExprT;		\
    return VExpr<ExprT>(ExprT(a.const_iter(),b.const_iter()));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VExpr<T1>,VIterConst<T2>,FNAME<typename T1::numtype, T2> > >  OP (const VExpr<T1>& a, const Vector<T2>& b)		\
{		\
    typedef VBinExprOp<VExpr<T1>,VIterConst<T2>,FNAME<typename T1::numtype, T2> > ExprT;		\
    return VExpr<ExprT>(ExprT(a,b.const_iter()));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VIterConst<T1>,VExpr<T2>,FNAME<T1, typename T2::numtype> > >  OP (const Vector<T1>& a, const VExpr<T2>& b)		\
{		\
    typedef VBinExprOp<VIterConst<T1>,VExpr<T2>,FNAME<T1, typename T2::numtype> > ExprT;		\
    return VExpr<ExprT>(ExprT(a.const_iter(),b));		\
}		\
		\
template<class T1, class T2>		\
inline VExpr<VBinExprOp<VExpr<T1>,VExpr<T2>,FNAME<typename T1::numtype, typename T2::numtype> > >		\
	OP (const VExpr<T1>& a, const VExpr<T2>& b)		\
{		\
    typedef VBinExprOp<VExpr<T1>,VExpr<T2>,FNAME<typename T1::numtype, typename T2::numtype> > ExprT;		\
    return VExpr<ExprT>(ExprT(a,b));		\
}		\


VecSpecialBinOps(ApPow, pow)

//*********************************************************
//		unary expressions (i.e. *=) (Internal class members)
//*********************************************************


//
//
//***Single Operations like abs(vec) or abs(VExpr)
//
//
//THe main macro

#define MakeUnaryOpAlive(name, op)														\
	template<class T1>																	\
	inline VExpr<VBinUnaryOp<VIterConst<T1>, name<T1> > > 	op(const Vector<T1>& a)				\
	{																					\
		typedef VBinUnaryOp<VIterConst<T1>,name<T1> > ExprT;							\
		return VExpr<ExprT>(ExprT(a.const_iter()));											\
	}																					\
																						\
	template<class T1>																	\
	inline VExpr<VBinUnaryOp<VExpr<T1>,name<typename T1::numtype> > > op(const VExpr<T1>& a)	\
	{																					\
		typedef VBinUnaryOp<VExpr<T1>,name<typename T1::numtype> > ExprT;				\
		return VExpr<ExprT>(ExprT(a));													\
	}																					\
																						\

MakeUnaryOpAlive(ApAbs, abs)
MakeUnaryOpAlive(ApCos, cos)
MakeUnaryOpAlive(ApSin, sin)
MakeUnaryOpAlive(ApTan, tan)
MakeUnaryOpAlive(ApExp, exp)
MakeUnaryOpAlive(ApAcos, acos)
MakeUnaryOpAlive(ApAsin, asin)
MakeUnaryOpAlive(ApAtan, atan)
MakeUnaryOpAlive(ApLog, log)
MakeUnaryOpAlive(ApCosh, cosh)
MakeUnaryOpAlive(ApSinh, sinh)
MakeUnaryOpAlive(ApTanh, tanh)
MakeUnaryOpAlive(ApAcosh, acosh)
MakeUnaryOpAlive(ApAsinh, asinh)
MakeUnaryOpAlive(ApAtanh, atanh)
MakeUnaryOpAlive(ApRe, Re)
MakeUnaryOpAlive(ApIm, Im)
MakeUnaryOpAlive(ApRe, real)
MakeUnaryOpAlive(ApIm, imag)
MakeUnaryOpAlive(ApConj, conj)
MakeUnaryOpAlive(ApCiel, ceil)
MakeUnaryOpAlive(ApFloor, floor)
#ifdef HAVE_FINITE
MakeUnaryOpAlive(ApNan, isnan)
#elif HAVE_ISNAN
MakeUnaryOpAlive(ApNan, isnan)
#endif
MakeUnaryOpAlive(ApSqrt, sqrt)
MakeUnaryOpAlive(ApSqr, sqr)
MakeUnaryOpAlive(ApCube, cube)
MakeUnaryOpAlive(ApForthPow, quartic)
MakeUnaryOpAlive(ApFifthPow, quintic)
MakeUnaryOpAlive(ApChop, chop)


//
//
//***Special 'total' vector operations
//
//

//************Dot product*********************

template<class T1, class T2>
inline OutType(typename T1::numtype,typename T2::numtype) ApDot(const T1& a, const T2 &b){

#ifdef VBoundCheck
	if(a.guessLen()!=b.guessLen()){
		BLEXCEPTION(" must be of the same length to do a dot product")
	}
#endif
	OutType(typename T1::numtype,typename T2::numtype) ss=a[0]*b[0];
	for(int i=1;i<a.guessLen();i++){
		ss+=a[i]*b[i];
	}
	return ss;
}

template<class T1, class T2>
inline OutType(T1,T2) dot(const Vector<T1>& a, const Vector<T2> &b){																					\
	return ApDot(a,b);
}

template<class T1, class T2>
inline OutType(typename T1::numtype,T2) dot(const VExpr<T1>& a, const Vector<T2> &b){																					\
	return ApDot(a,b);
}

template<class T1, class T2>
inline OutType(T1,typename T2::numtype) dot(const Vector<T1>& a, const VExpr<T2> &b){																					\
	return ApDot(a,b);
}

template<class T1, class T2>
inline OutType(typename T1::numtype,typename T2::numtype) dot(const VExpr<T1>& a, const VExpr<T2> &b){																					\
	return ApDot(a,b);
}

//******************Equality type operations...VECTOR***********************



template<class T1, class T2>
inline bool ApEqualVec(const T1& a, const T2 &b){
	if(a.guessLen()!=b.guessLen()) return false;
	for(int i=0;i<a.guessLen();i++){
		if(b(i)!=a(i)) return false;
	}
	return true;
}

template<class T1, class T2>
inline bool ApNotEqualVec(const T1& a, const T2 &b){
	if(a.guessLen()!=b.guessLen()) return true;
	for(int i=0;i<a.guessLen();i++){
		if(b(i)==a(i)) return false;
	}
	return true;
}

template<class T1, class T2>
inline bool ApGreaterThenVec(const T1& a, const T2 &b){
	if(a.guessLen()!=b.guessLen()) return false;
	for(int i=0;i<a.guessLen();i++){
		if(a(i)<=b(i)) return false;
	}
	return true;
}

template<class T1, class T2>
inline bool ApGreaterThenEqualVec(const T1& a, const T2 &b){
	if(a.guessLen()!=b.guessLen()) return false;
	for(int i=0;i<a.guessLen();i++){
		if(a(i)<b(i)) return false;
	}
	return true;
}

template<class T1, class T2>
inline bool ApLessThenVec(const T1& a, const T2 &b){
	if(a.guessLen()!=b.guessLen()) return false;
	for(int i=0;i<a.guessLen();i++){
		if(a(i)>=b(i)) return false;
	}
	return true;
}

template<class T1, class T2>
inline bool ApLessThenEqualVec(const T1& a, const T2 &b){
	if(a.guessLen()!=b.guessLen()) return false;
	for(int i=0;i<a.guessLen();i++){
		if(a(i)>b(i)) return false;
	}
	return true;
}

#define EqTypeVecOps(FNAME,OP) \
template<class T1>		\
inline bool operator OP (const Vector<T1> & a , const Vector<T1> &b){	\
	return FNAME(a,b);	\
}	\
\
template<class T1>		\
inline bool operator OP (const Vector<T1> & a , const VExpr<T1> &b){	\
	return FNAME(a,b);	\
}	\
\
template<class T1>		\
inline bool operator OP (const VExpr<T1> & a , const Vector<T1> &b){	\
	return FNAME(a,b);	\
}	\
\
template<class T1>		\
inline bool operator OP (const VExpr<T1> & a , const VExpr<T1> &b){	\
	return FNAME(a,b);	\
}	\
\
template<class T1, class T2>		\
inline bool operator OP (const Vector<T1> & a , const Vector<T2> &b){	\
	return FNAME(a,b);	\
} \
\
template<class T1, class T2>		\
inline bool operator OP (const VExpr<T1> & a , const Vector<T2> &b){	\
	return FNAME(a,b);	\
} \
\
template<class T1, class T2>		\
inline bool operator OP (const Vector<T1> & a , const VExpr<T2> &b){	\
	return FNAME(a,b);	\
} \
\
template<class T1, class T2>		\
inline bool operator OP (const VExpr<T1> & a , const VExpr<T2> &b){	\
	return FNAME(a,b);	\
} \


EqTypeVecOps(ApEqualVec, ==)
EqTypeVecOps(ApNotEqualVec, !=)
EqTypeVecOps(ApGreaterThenVec, >)
EqTypeVecOps(ApGreaterThenEqualVec, >=)
EqTypeVecOps(ApLessThenVec, <)
EqTypeVecOps(ApLessThenEqualVec, <=)

//comparisons for number types...
#define VectCompNum(TYPE,NUM) \
	template<class T1>	\
	inline bool  operator == (const TYPE<T1> & a, const NUM & b)	\
	{	\
		for(int i=0;i<a.guessLen();++i)	\
			if(a[i] != b) return false; 	\
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator != (const TYPE<T1> & a, const NUM & b)	\
	{	\
		for(int i=0;i<a.guessLen();++i)	\
			if(a[i] == b) return false; \
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator >= (const TYPE<T1> & a, const NUM & b)	\
	{	\
		for(int i=0;i<a.guessLen();++i)	\
			if(a[i] < b) return false; 	\
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator <= (const TYPE<T1> & a, const NUM & b)	\
	{	\
		for(int i=0;i<a.guessLen();++i)	\
			if(a[i] > b) return false; \
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator < (const TYPE<T1> & a, const NUM & b)	\
	{	\
		for(int i=0;i<a.guessLen();++i)	\
			if(a[i] >= b) return false; 	\
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator > (const TYPE<T1> & a, const NUM & b)	\
	{	\
		for(int i=0;i<a.guessLen();++i)	\
			if(a[i] <= b) return false; 	\
		return true;	\
	}	\
	\
	template<class T1>	\
	inline bool  operator == (const NUM & a, const TYPE<T1> & b)	\
	{	\
		for(int i=0;i<b.guessLen();++i)	\
			if(a != b[i]) return false; 	\
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator != (const NUM & a, const TYPE<T1> & b)	\
	{	\
		for(int i=0;i<b.guessLen();++i)	\
			if(a == b[i]) return false; \
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator >= (const NUM & a, const TYPE<T1> & b)	\
	{	\
		for(int i=0;i<b.guessLen();++i)	\
			if(a < b[i]) return false; 	\
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator <= (const NUM & a, const TYPE<T1> & b)	\
	{	\
		for(int i=0;i<b.guessLen();++i)	\
			if(a > b[i]) return false; \
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator < (const NUM & a, const TYPE<T1> & b)	\
	{	\
		for(int i=0;i<b.guessLen();++i)	\
			if(a >= b[i]) return false; 	\
		return true;	\
	} \
	template<class T1>	\
	inline bool  operator > (const NUM & a, const TYPE<T1> & b)	\
	{	\
		for(int i=0;i<b.guessLen();++i)	\
			if(a <= b[i]) return false; 	\
		return true;	\
	}	\

#ifdef HAVE_ISNAN
	template<class T1>
	inline bool hasnan(const Vector<T1> &b)
	{
		for(int i=0;i<b.guessLen();++i)
			if(hasnan(b[i])) return true ;
		return false;
	}

	template<class T1>
	inline bool hasnan(const VExpr<T1> &b)
	{
		for(int i=0;i<b.guessLen();++i)
			if(hasnan(b[i])) return true ;
		return false;
	}

#endif



VectCompNum(Vector, double)
VectCompNum(Vector, float)
VectCompNum(Vector, complex)
VectCompNum(Vector, int)
VectCompNum(Vector, unsigned int)
VectCompNum(Vector, long)
VectCompNum(Vector, unsigned long)
VectCompNum(Vector, short)
VectCompNum(Vector, char)
VectCompNum(Vector, unsigned char)
VectCompNum(Vector, bool)

VectCompNum(VExpr, double)
VectCompNum(VExpr, float)
VectCompNum(VExpr, complex)
VectCompNum(VExpr, int)
VectCompNum(VExpr, unsigned int)
VectCompNum(VExpr, long)
VectCompNum(VExpr, unsigned long)
VectCompNum(VExpr, short)
VectCompNum(VExpr, char)
VectCompNum(VExpr, unsigned char)
VectCompNum(VExpr, bool)


//******************Min***********************

template<class T1>
inline typename T1::numtype ApMin(const T1& a){
	/*if(a.guessLen()<=0){
		std::cerr<<std::endl<<"Error: min(t1)"<<std::endl;
		std::cerr<<" length of object MUST be at least 1"<<std::endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}*/
	typename T1::numtype ss=a[0];
	for(int i=1;i<a.guessLen();i++){
		if(min(ss)>=min(a[i])) ss=a[i];
	}
	return ss;
}

template<class T1>
inline T1 min(const Vector<T1>& a){																					\
	return ApMin(a);
}

template<class T1>
inline typename T1::numtype min(const VExpr<T1>& a){																					\
	return ApMin(a);
}

//******************Max***********************

template<class T1>
inline typename T1::numtype ApMax(const T1& a){

#ifdef VBoundCheck
	if(a.guessLen()<=0){
		BLEXCEPTION(" length of object MUST be at least 1")
	}
#endif

	typename T1::numtype ss=a[0];
	for(int i=1;i<a.guessLen();i++){
		if(max(ss)<=max(a[i])) ss=a[i];
	}
	return ss;
}

template<class T1>
inline T1 max(const Vector<T1>& a){																					\
	return ApMax(a);
}

template<class T1>
inline typename T1::numtype max(const VExpr<T1>& a){																					\
	return ApMax(a);
}

//**************Norm************************
template<class T1>
inline FloatType(typename T1::numtype) ApNorm(const T1 &a){
    int length = a.guessLen();

#ifdef VBoundCheck
    if(length<=0){
		BLEXCEPTION(" length of object MUST be at least 1")
	}
#endif

	typedef typename T1::numtype numtype;
    typedef SumType(numtype)          Stype;
    typedef FloatType(numtype)        Ftype;

    Ftype vv=a[0]*a[0];
    Ftype ss =vv;

    for (int i=1; i < length; ++i){
		vv=Ftype(a[i]);
		ss += vv * vv;
    }
    return ApSqrt<Ftype>::apply(ss);
}

template<class T1>
inline FloatType(T1) norm(const Vector<T1>& a){																					\
	return ApNorm(a);
}

template<class T1>
inline FloatType(typename T1::numtype) norm(const VExpr<T1>& a){																					\
	return ApNorm(a);
}

//**************SqaureNorm************************
template<class T1>
inline FloatType(typename T1::numtype) ApSqNorm(const T1 &a){
    int length = a.guessLen();

#ifdef VBoundCheck
    if(length<=0){
		BLEXCEPTION(" length of object MUST be at least 1")
	}
#endif

	typedef typename T1::numtype numtype;
    typedef SumType(numtype)          Stype;
    typedef FloatType(numtype)        Ftype;

    Ftype vv=a[0]*a[0];
    Ftype ss =vv;

    for (int i=1; i < length; ++i){
		vv=Ftype(a[i]);
		ss += vv * vv;
    }
    return ss;
}

template<class T1>
inline FloatType(T1) square_norm(const Vector<T1>& a){																					\
	return ApSqNorm(a);
}

template<class T1>
inline FloatType(typename T1::numtype) square_norm(const VExpr<T1>& a){																					\
	return ApSqNorm(a);
}

//**************SUM************************
template<class T1>
inline SumType(typename T1::numtype) ApSum(const T1 &a){
    int length = a.guessLen();
	typedef typename T1::numtype numtype;
    typedef SumType(numtype)          Stype;
    Stype hh=Stype(0);
    for (int i=0; i < length; ++i){
			hh += a(i) ;
    }
    return hh;
}

template<class T1>
inline SumType(T1) sum(const Vector<T1>& a){																					\
	return ApSum(a);
}

template<class T1>
inline Vector<SumType(T1)> sum(const Vector<Vector<T1> >& a, int dir=1){																					\
    Vector<SumType(T1)> hh;
    if(dir==1){
		if(a.size()>0) hh.resize(a[0].size());
		hh.fill(SumType(T1)(0));
		for (int i=0; i < a.size(); ++i){
			hh += a(i) ;
		}
	}else{
		hh.resize(a.size());
		hh.fill(SumType(T1)(0));
		for (int i=0; i < a.size(); ++i){
			hh(i) += sum(a(i)) ;
		}
	}
    return hh;
}

template<class T1>
inline SumType(typename T1::numtype) sum(const VExpr<T1>& a){																					\
	return ApSum(a);
}

//**************Product************************
template<class T1>
inline FloatType(typename T1::numtype) ApProdv(const T1 &a){
    int length = a.guessLen();
	typedef typename T1::numtype numtype;
    typedef FloatType(numtype)          Stype;
    Stype hh=Stype(1);
    for (int i=0; i < length; ++i){
			hh *= a(i) ;
    }
    return hh;
}

template<class T1>
inline FloatType(T1) prod(const Vector<T1>& a){																					\
	return ApProdv(a);
}

template<class T1>
inline FloatType(typename T1::numtype) prod(const VExpr<T1>& a){																					\
	return ApProdv(a);
}

//***********Chop*************************
template<class T1>
inline void chop(Vector<T1> &a)
{
	for (int i=0; i < a.guessLen(); ++i)
	{
		chop(a(i), 1.e-12);
	}
}

template<class T1>
inline void chop(VExpr<T1> &a)
{
	for (int i=0; i < a.guessLen(); ++i)
	{
		chop(a(i), 1.e-12);
	}
}

template<class T1>
inline void chop(Vector<T1> &a, double lim)
{
	for (int i=0; i < a.guessLen(); ++i)
	{
		chop(a(i), lim);
	}
}

template<class T1>
inline void chop(VExpr<T1> &a, double lim)
{
	for (int i=0; i <a.guessLen(); ++i)
	{
		chop(a(i), lim);
	}
}
//**************IO*******************

template<class T>
void Vector<T>::print(std::ostream &out)
{	print(out, " ");	}

template<class T>
template<class NextThing>
void Vector<T>::print(std::ostream &otr, NextThing out)
{
	for (int i=0; i < length(); ++i)
		otr <<get(i)<<out;//x.RecordSep;
}


template<class T>
std::ostream& operator<<(std::ostream& otr, const Vector<T>& x){
	for (int i=0; i < x.length(); ++i)
		otr <<x[i]<<" ";//x.RecordSep;

    return otr;
}

template<class expr>
std::ostream& operator<<(std::ostream& otr, const VExpr<expr> &oo){
	for (int i=0; i < oo.guessLen(); ++i)	otr<<oo[i]<<" ";
//	otr<<tmp.RecordSep;
    return otr;
}

END_BL_NAMESPACE


#include "Vectoraux.h"

#endif

