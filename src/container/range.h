 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-16-01
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
	range.h--> allows for odd acsess to vectors that do not start at one
*/

//examples;;;
//Vector<double> a(6,0), b(3, 3), d(6,0),c;
//c=a(Range(3,6))+b;  //c is a length 3 vector


#ifndef _range_h_
#define _range_h_ 1

//#include "blochconfig.h"
#ifndef ON_WINDOWS
#include "blochconfig.h"
#endif

#ifdef HAVE_CLIMITS
 #include <climits>
#else
 #include <limits.h>
#endif

#include "container/rankType.h"
#include <iostream>

BEGIN_BL_NAMESPACE




class Range {
	private:
		int fir_; //first element
		int las_; //last element
		int str_; //stride

	public:
		typedef int numtype;
		enum {Start=INT_MIN, End=INT_MIN};

		static Range All();

		Range():
			fir_(Start), las_(End), str_(1)
		{}

		Range(int f, int l):
			fir_(f), las_(l), str_(1)
		{}

		Range(int f, int l, int s);

		inline int first() const
		{
			return fir_;
		}

		inline int last() const
		{
			return las_;
		}

		inline int first(int guesff) const
		{
			if(fir_==Start) return guesff;
			return fir_;
		}

		inline int begin(int gb) const {	return first(gb);	}

		inline int start(int guesff) const
		{
			if(fir_==Start) return guesff;
			return fir_;
		}

		inline int end(int guesll) const
		{
			if(las_==End) return guesll;
			return las_;
		}

		inline int last(int guesll) const
		{
			if(las_==End) return guesll;
			return las_;
		}

		inline int stride() const
		{
			return str_;
		}

		inline int length() const
		{
			return std::abs(fir_-las_)/str_ +1;
		}

		inline int length(int len) const
		{
			if(fir_==Start && las_==End) return len;
			if(fir_==Start) return std::abs(las_)/str_+1;
			return std::abs(fir_-las_)/str_ +1;
		}

		inline int guessLen() const
		{
			if(fir_==Start) return (las_)/str_ +1;
			return std::abs(las_-fir_)/str_ +1;
		}

		inline int guessLen(int len) const
		{
			if(fir_==Start && las_==End) return len;
			if(fir_==Start) return std::abs(las_)/str_+1;
			return abs(fir_-las_)/str_ +1;
		}

		inline int size() const
		{
			return abs(fir_-las_)/str_ +1;
		}

		// Operators
		Range operator-(int shift) const;
		Range operator+(int shift) const;
		void operator=(const Range &rhs);

		inline int operator[](unsigned i) const
		{		return first(0) + i * str_; 	}

		inline int operator()(unsigned i) const
		{		return first(0) + i * str_;		}


		friend std::ostream& operator<<(std::ostream& os, const Range& range);

};

template<class Num_t>
class Spread {
	private:
		Num_t fir_; //first element
		Num_t las_; //last element
		Num_t str_; //stride



	public:
		typedef Num_t numtype;

		Spread():
			fir_(0), las_(1), str_(1)
		{}

		Spread(Num_t f, Num_t l):
			fir_(f), las_(l), str_(1)
		{}

		Spread(Num_t f, Num_t l, Num_t s):
			fir_(f), las_(l), str_(s)
		{}

		Spread(const Spread &cp):
			fir_(cp.fir_), las_(cp.las_), str_(cp.str_)
		{}

		Spread(const Range &cp):
			fir_(cp.fir_), las_(cp.las_), str_(cp.str_)
		{}

		Spread &operator=(const Spread &rhs)
		{
			if(this==&rhs) return *this;
			fir_=rhs.fir_;
			las_=rhs.las_;
			str_=rhs.str_;
			return *this;
		}

		Spread &operator=(Range &rhs)
		{
			fir_=numtype(rhs.fir_);
			las_=numtype(rhs.las_);
			str_=numtype(rhs.str_);
			return *this;
		}

		template<class T>
		Spread &operator=(Spread<T> &rhs)
		{
			fir_=numtype(rhs.fir_);
			las_=numtype(rhs.las_);
			str_=numtype(rhs.str_);
			return *this;
		}

		inline Num_t first() const
		{
			return fir_;
		}

		inline Num_t last() const
		{
			return las_;
		}

		inline int begin(int gb) const {	return 0;	}

		inline Num_t first(int guesff) const
		{
			return fir_;
		}

		inline Num_t start(int guesff) const
		{
			return fir_;
		}

		inline int end(int guesll) const
		{
			return length(guesll);
		}

		inline Num_t last(int guesll) const
		{
			return las_;
		}

		inline Num_t stride() const
		{
			return str_;
		}

		inline int length() const
		{
			return int(std::abs(fir_-las_)/str_) +1;
		}

		inline int length(int len) const
		{
			return int(std::abs(fir_-las_)/str_) +1;
		}

		inline int guessLen() const
		{
			return int(std::abs(las_-fir_)/str_) +1;
		}

		inline int guessLen(int len) const
		{
			return int(std::abs(fir_-las_)/str_) +1;
		}

		inline int size() const
		{
			return int(std::abs(fir_-las_)/str_) +1;
		}

		// Operators
		Spread operator-(int shift) const
		{
			return Spread(fir_ - shift, las_ - shift, str_);
		}

		Spread operator+(int shift) const
		{
			return Spread(fir_ + shift, las_ + shift, str_);
		}

		numtype operator[](unsigned i) const
		{
			return fir_ + i * str_;
		}

		numtype operator()(unsigned i) const
		{
			return fir_ + i * str_;
		}


		friend inline std::ostream& operator<<(std::ostream& os, const Spread& range)
		{
			os << "Spread(" << range.first() << "," << range.last() << ","<< range.stride() << ")";

			return os;
    	}

};

END_BL_NAMESPACE


#endif



