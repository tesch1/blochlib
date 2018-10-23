/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-25-01
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
	powder.h--> generates the powder angles from various methods
	like, zcw (the best one), alderman, sophie, planar, spherical,
	and 'none' for liquid type sims
*/
/* Here we have the general powder avearge generators
it calculates a vector of theta, phi, and weight vectors
OR it reads a file of points and weights (if there are no
weights specified it assumes a '1')*/

#ifndef _powder_h_
#define _powder_h_ 1

#include "container/Vector/Vector.h"
#include "container/range.h"
#include "container/matrix/matrix.h"
#include<string>
#include<stdlib.h>
#include "utils/constants.h"
#include "utils/utils.h"
#include "utils/blassert.h"

BEGIN_BL_NAMESPACE


const std::string powtypes[]={"alderman", "sphere", "rect", "zcw", "sophie", "liquid"};

class powderIter;
class constpowderIter;

class powder {
	friend class powderIter;
	friend class constpowderIter;
	private:
		Vector<double> _theta;
		Vector<double> _phi;
		Vector<double> _weight;
		Vector<double> _gamma;

		int _tstep;
		int _pstep;
		int _gstep;
	public:
		enum pTypes{alderman, sphere, rect, zcw, sophie, liquid,None};

	private:
		pTypes _method;
		std::string _fname;

		int _ct;
		bool _underlimit;
		void __reset();


	public:

		typedef powderIter iterator;
		typedef constpowderIter const_iterator;

	//given a string will tell me if it is a valid powder algorithm..
		pTypes AssignM(std::string me);
		pTypes canCalculate(std::string me); //same function as above

		powder();
		powder(std::string method);
		powder(pTypes method);
		powder(std::string method, int tstep, int pstep, int gammast=1);
		powder(pTypes method, int tstep, int pstep, int gammast=1);
		powder(double oneth, double oneph, double oneg=0.0);
		powder(double oneth);
		powder(int size);

		void calcPow();
		void calcPow(std::string method, int tstep, int pstep, int gammast=1);
		void calcPow(pTypes method, int tstep, int pstep, int gammast=1);

		inline void reset(){	_ct=0; _underlimit=true;	}

		void read(std::string fname);
		void read(const char *fname);
		void read(std::ifstream &infile);

		inline int size()const{ return _theta.size();	}

		void operator=(const powder &rhs);

//		inline iterator begin(){	return iterator(*this);	}
//		inline const_iterator begin()const{	return const_iterator(*this);	}

		inline double getTheta(int i){
			if(this->size()>i && !_theta.empty())
				return _theta[i];
			else
				return 0.;
		}

		inline double getPhi(int i){
			if(this->size()>i && !_phi.empty())
				return _phi[i];
			else
				return 0.;
		}

		inline double getGamma(int i){
			if(this->size()>i && !_gamma.empty())
				return _gamma[i];
			else
				return 0.;
		}

		inline double getWeight(int i){
			if(this->size()>i && !_theta.empty())
				return _weight[i];
			else
				return 0.;
		}
		inline double theta(int i){	return getTheta(i);		}
		inline double phi(int i){	return getPhi(i);	}
		inline double gamma(int i){	return getGamma(i);	}
		inline double weight(int i){	return getWeight(i);	}

		operator bool() const { return _underlimit;	}
		void operator++(){
			if(!_underlimit){
				BLEXCEPTION(" itterating too far...")
			}
			_ct++;
			if(_ct>=size())	_underlimit=false;
		}

		inline double theta(){	return _theta[_ct];	}
		inline double phi(){		return _phi[_ct];	}
		inline double gamma(){		return _gamma[_ct];	}
		inline double weight(){	return _weight[_ct];	}

		void zero() {
			_ct=0;
			_underlimit=true;
		}

};

class constpowderIter{
	private:
		const powder *mp;
		int curpos_;
		int end_;
		int begin_;
		bool notended_;

	public:
		constpowderIter():
			curpos_(0), end_(0), begin_(0), notended_(false)
		{
			mp=NULL;
		}

		constpowderIter(const powder &in):
			curpos_(0), end_(in.size()), begin_(0), notended_(true)
		{
			mp=&in;
		}

		constpowderIter(const powder &in, Range r):
			curpos_(r.first(0)), end_(r.last(in.size())), begin_(r.first(0)), notended_(true)
		{
			mp=&in;
		}

		constpowderIter(const powder &in, int b, int e):
			curpos_(b), end_(e), begin_(b), notended_(true)
		{
			mp=&in;
		}


		constpowderIter &operator=(const constpowderIter &rhs)
		{
			if(&rhs==this)	return *this;
			mp=rhs.mp;
			begin_=rhs.begin_;
			curpos_=rhs.curpos_;
			end_=rhs.end_;
			return *this;
		}

		~constpowderIter(){
			mp=NULL;
		}

		operator bool() const { return notended_;	}
		void operator++(){
			if(!notended_){
				BLEXCEPTION(" itterating too far...")
			}
			curpos_++;
			if(curpos_>=end_)	notended_=false;
		}

		int curpos()const{	return curpos_;	}

		inline double theta() const {	return mp->_theta[curpos_];	}
		inline double phi() const {		return mp->_phi[curpos_];	}
		inline double gamma() const {		return mp->_gamma[curpos_];	}
		inline double weight() const {	return mp->_weight[curpos_];	}

		int end()	const	{	return end_;	}
		int begin()	const	{	return begin_;	}
		void begin(int in)	{	RunTimeAssert(in>0); begin_=in; curpos_=in;	}
		void end(int in)	{	RunTimeAssert(in<size()); end_=in;	}
		int size() const {	return mp->size();	}
};

class powderIter{
	private:
		powder *mp;
		int curpos_;
		int end_;
		int begin_;
		bool notended_;

	public:
		powderIter():
			curpos_(0), end_(0), begin_(0), notended_(false)
		{
			mp=NULL;
		}

		inline powderIter( powder &in):
			curpos_(0), end_(in.size()), begin_(0),notended_(true)
		{
			mp=&in;
		}

		inline powderIter(powder &in, Range r):
			curpos_(r.first(0)), end_(r.last(in.size())), begin_(r.first(0)), notended_(true)
		{
			mp=&in;
		}

		inline powderIter(powder &in, int b, int e):
			curpos_(b), end_(e), begin_(b), notended_(true)
		{
			mp=&in;
		}


		void operator=( powderIter rhs)
		{
			if(&rhs==this)	return;
			mp=rhs.mp;
			begin_=rhs.begin_;
			curpos_=rhs.curpos_;
			end_=rhs.end_;
			notended_=rhs.notended_;
		}



		~powderIter()
		{
			mp=NULL;
		}

		operator bool()  { return notended_;	}
		void operator++(){
			if(!notended_){
				BLEXCEPTION(" itterating too far...")
			}
			curpos_++;
			if(curpos_>=end_)	notended_=false;
		}

		inline int curpos()const{	return curpos_;	}

		inline double &theta()  {	return mp->_theta[curpos_];	}
		inline double &phi()  {		return mp->_phi[curpos_];	}
		inline double &gamma()  {		return mp->_gamma[curpos_];	}
		inline double &weight()  {	return mp->_weight[curpos_];	}

		inline int end()	const	{	return end_;	}
		inline int begin()	const	{	return begin_;	}

		inline int &end()		{	return end_;	}
		inline int &begin()		{	return begin_;	}
		void begin(int in)	{	RunTimeAssert(in>=0); begin_=in; curpos_=in;	}
		void end(int in)	{	RunTimeAssert(in<size()); end_=in;	}

		int size()  {	return mp->size();	}
};

//Display operator...
template<class Engine_t>
std::ostream &operator << (std::ostream &oo,powder &out);

//ASCII WRITE
template<class Engine_t>
std::ofstream &operator << (std::ofstream &OutFile,powder &ObjToWrite);


//BINARY WRITE
template<class Engine_t>
std::fstream &operator << (std::fstream &in, powder &out);


//BINARY READ
template<class Engine_t>
std::fstream &operator >> (std::fstream &in, powder &out);

//ASCII READ
template<class Engine_t>
std::ifstream &operator >> (std::ifstream &in, powder &out);


END_BL_NAMESPACE



#endif
