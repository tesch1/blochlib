/* listblochpars.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-10-01
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
 	listblochpars.h-->a large list of BlochParams<>  with Grids attached...
 */


#ifndef _listblochpars_h_
#define _listblochpars_h_ 1

#include "bloch/blochParams.h"
#include "bloch/blochParamsBasic.h"
#include "bloch/gradientgrid.h"
#include "container/grids/grids.h"

BEGIN_BL_NAMESPACE


/* maintains a list of bloch parameters
	there are two specializations,
		one with no grid (a 'NullGrid')
		and one with a GradientGrid (a grid snapped with a graident vector)
	all others inlcude simple grids.

*/


enum SolveUnits_{Dimensionless, Normal};

template<class GridEngine_t=NullGrid, int BPops=BPoptions::Density | BPoptions::HighField, class offset_T=double>
class ListBlochParams;

template<class GridEngine_t, int BPops, class offset_T>
class ListBlochParams<GradientGrid<GridEngine_t >, BPops, offset_T >;

//pre dec for typdef in bloch param...
template<class GridEngine_t=NullGrid, int BPops=BPoptions::Density | BPoptions::HighField, class offset_T=double>
class ListBlochParamsIterator;

template<int BPops, class offset_T>
class ListBlochParamsIterator<NullGrid, BPops, offset_T>;

template<class Engine_t,int BPops, class offset_T>
class ListBlochParamsIterator<GradientGrid<Engine_t>, BPops, offset_T >;


// The 'NullGrid' class...no other parameters other then a list BlochParams<>
template<int BPops, class offset_T >
class ListBlochParams<NullGrid, BPops, offset_T> {

	//friend class ListBlochParamsIterator<NullGrid, BPops, offset_T>;

	private:

		Vector<BlochParams<BPops, offset_T> > Vpars_;


		double totMo_;

		SolveUnits_ MySolveUnits_;


		//this finds the largest 'Mo' in Normal/real units of all the spins in the
		//list...this is used for setting a Dimensionless Mo where one could have
		//different concnetrations and gamma factors for each spin...
		double FindMaxMo();


		void reCalcMo();

	public:

		typedef ListBlochParamsIterator<NullGrid, BPops, offset_T> iterator;
		typedef BlochParams<BPops, offset_T> Parameter_T;
		typedef offset_T Offset_T;

		//the flag from 'InitialCondition'
		int InitCond;

		ListBlochParams():
			Vpars_(1, BlochParams<BPops, offset_T> ()),  totMo_(0),
			MySolveUnits_(Normal),
			InitCond(InitialCondition::AllUp)
		{
			calcTotalMo();
		}

		ListBlochParams(int nsp,int IntC=InitialCondition::AllUp, SolveUnits_ mysu=Normal):
			Vpars_(nsp, BlochParams<BPops, offset_T> ()), totMo_(0),
			MySolveUnits_(mysu),
			InitCond(IntC)
		{
			calcTotalMo();
		}

		ListBlochParams(int nsp, std::string inspin ,int IntC=InitialCondition::AllUp,SolveUnits_ mysu=Normal):
			Vpars_(nsp, BlochParams<BPops, offset_T> (inspin)),totMo_(0),
			MySolveUnits_(mysu),
			InitCond(IntC)
		{
			calcTotalMo();
		}


		ListBlochParams(const BlochParams<BPops, offset_T>  &one, int IntC=InitialCondition::AllUp,SolveUnits_ mysu=Normal):
			Vpars_(1, one),  totMo_(0),
			MySolveUnits_(mysu),
			InitCond(IntC)
		{
			calcTotalMo();
		}

		ListBlochParams(int nsp, const BlochParams<BPops, offset_T>  &dup, int IntC=InitialCondition::AllUp, SolveUnits_ mysu=Normal):
			Vpars_(nsp, dup), totMo_(0),
			MySolveUnits_(mysu),
			InitCond(IntC)
		{
			calcTotalMo();
		}

		ListBlochParams( const ListBlochParams &dup):
			Vpars_(dup.Vpars_),  totMo_(0),
			MySolveUnits_(dup.MySolveUnits_),
			InitCond(dup.InitCond)
		{
			calcTotalMo();
		}

		~ListBlochParams()	{}

		iterator begin();

		inline int size()	const	{	return Vpars_.size();	}

		void calcTotalMo();

		//uses the 'InitialCondition' class above
		void setInitialCondition(int IC){	InitCond=IC;	calcTotalMo(); }
		void calcInitialCondition(){	calcTotalMo();	}
		void calcInitialCondition(int IC){	InitCond=IC; calcTotalMo();	}
		inline int initialCondition(){	return InitCond;	}

		SolveUnits_ &SolveUnits();
		void SetSolveUnits(const SolveUnits_ in);


		inline offset_T Bo() const;
		inline offset_T &Bo();
		inline offset_T &Bo(int i);
		inline offset_T Bo(int i) const;
		inline void Bo(offset_T ne);

		inline offset_T &Mo(int i)	;
		inline offset_T Mo(int i)	const;

		inline double &TotMo() ;
		inline double TotMo() const;
		inline double &totalMo() ;
		inline double totalMo() const;

		inline double meanMo() const;

		inline double &temperature();
		inline double temperature() const;
		inline void temperature(double ne);
		inline double &temperature(int i);
		inline double temperature(int i) const ;

		inline offset_T &omega(int i) ;
		inline offset_T omega(int i) const;

		inline double &moles(int i);
		inline double moles(int i) const;

		inline double &offset(int i);
		inline double &offset(int i, double t); //time dependant grids

		//this 'offset' function returns the offset WITHOUT a gradient
		//so that one can 'set' the value of it (the BlochParameter offset)
		inline double &spin_offset(int i);

		std::string symbol(int i) const ;

		inline double gamma(int i) const;
		inline double gammaGauss(int i) const;

		inline double offset(int i) const;
		inline double offset(int i, double t) const; //time dependant grids

		//this 'offset' function returns the offset WITHOUT a gradient
		//so that one can 'set' the value of it (the BlochParameter offset)
		inline double spin_offset(int i) const;

		inline rmatrix Jacobian(int i);
		inline void Jacobian(rmatrix &out, int i);

		Parameter_T  &operator()(int i)	;

		Parameter_T  &operator[](int i)	;


		void operator=(const ListBlochParams &rhs);

		void operator=(const BlochParams<BPops, offset_T>  &rhs);

		ListBlochParams &operator+(ListBlochParams &rhs);

		ListBlochParams &operator+(const BlochParams<BPops, offset_T>  &rhs);

		void operator+=( ListBlochParams &rhs);
		void operator+=(const BlochParams<BPops, offset_T>  &rhs);

		void push_back(const BlochParams<BPops, offset_T>  &rhs);


		void print(std::ostream &oo) ;

		void print()  ;	//text print
		/* The entire list is saved like so...

		----------
		LISTBLOCHPARS   -->this is in ASCII
		totalbyte Bo temperature numberofspins
		pars1 pars2 ... parsN
		ENDLISTBLOCHPARS --->this is in ASCII
		--------------
		The first and last words allow us to find the parameters
		if nested in a larger file of other data (binary or otherwise sooo DONOT use those words for
		anything else...

		*/


		bool write(std::fstream &oo); //binary write

		bool read(std::fstream &in); //binary read
};

template<int BPops, class offset_T>
std::ostream& operator<<(std::ostream &oo,  ListBlochParams<NullGrid, BPops, offset_T> &out);

/*****
the class that uses a grid..but NOT gradients....
*****/

template<class GridEngine_t, int BPops, class offset_T>
class ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >; //pre dec for the construcutors and '=' for the class below

template<class GridEngine_t, int BPops, class offset_T>
class ListBlochParams : public ListBlochParams<NullGrid, BPops, offset_T> {

	//friend class ListBlochParamsIterator<GradientGrid<Engine_t> >;

	private:

		GridEngine_t *grads_;

	public:

		typedef ListBlochParamsIterator< GridEngine_t, BPops, offset_T > iterator;
		typedef typename GridEngine_t::iterator Griditerator;

		ListBlochParams(): ListBlochParams<NullGrid, BPops, offset_T>(),
			grads_(NULL)
		{}

		ListBlochParams(int nsp, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(nsp, IntC),
			grads_(NULL)
		{}

		ListBlochParams(int nsp, std::string inspin, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(nsp, inspin, IntC),
			grads_(NULL)
		{}

		ListBlochParams(int nsp, std::string inspin, GridEngine_t &gr, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(nsp, inspin, IntC)
		{
			grads_=&gr;
		}

		ListBlochParams(const BlochParams<BPops, offset_T>  &onein, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(onein, IntC),
			grads_(NULL)
		{}

		ListBlochParams(int nsp, const BlochParams<BPops, offset_T>  &dup, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(nsp,dup,IntC),
			grads_(NULL)
		{}

		ListBlochParams( const ListBlochParams<NullGrid, BPops, offset_T> &dup, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(dup),
			grads_(NULL)
		{}

		ListBlochParams( const ListBlochParams<NullGrid, BPops, offset_T> &dup,GridEngine_t &gr, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(dup),
			grads_(&gr)
		{}

		ListBlochParams( const ListBlochParams &dup):
			ListBlochParams<NullGrid, BPops, offset_T>(dup),
			grads_(dup.grads_)
		{}

	//defined in the methods
		ListBlochParams( const ListBlochParams<GradientGrid<GridEngine_t >, BPops, offset_T > &dup):
			ListBlochParams<NullGrid, BPops, offset_T>(dup),
			grads_(dup.Grid())
		{}


		~ListBlochParams(){ grads_=NULL;	}

		iterator begin(){	return iterator(*this);	}

		void SetGrid(GridEngine_t &in){	grads_=&in;	}

		inline GridEngine_t *Grid();
		inline const GridEngine_t *Grid()	 const;

		void operator=(const ListBlochParams &rhs);

		void operator=(const ListBlochParams<NullGrid, BPops, offset_T> &rhs);

		//template<class Engine_t>
		void operator=(const ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > &rhs);
};


template<class GridEngine_t, int BPops, class offset_T>
std::ostream& operator<<(std::ostream &oo,  ListBlochParams<GridEngine_t > &out);



/* the method for using the parameters with Gradients...
uses the above 'NullGrid' class as its base

*/

template<class Engine_t, int BPops, class offset_T>
class ListBlochParams< GradientGrid<Engine_t>, BPops, offset_T> : public ListBlochParams<NullGrid, BPops, offset_T> {

	//friend class ListBlochParamsIterator<GradientGrid<Engine_t> >;

	private:

		GradientGrid<Engine_t> *grads_;

		bool apply_;
	public:

		typedef ListBlochParamsIterator<GradientGrid<Engine_t> , BPops, offset_T> iterator;

		ListBlochParams(): ListBlochParams<NullGrid, BPops, offset_T>(),
			grads_(NULL), apply_(true)
		{}

		ListBlochParams(int nsp, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(nsp, IntC),
			grads_(NULL), apply_(true)
		{}

		ListBlochParams(int nsp, std::string inspin, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(nsp, inspin, IntC),
			grads_(NULL), apply_(true)
		{}

		ListBlochParams(int nsp, std::string inspin, GradientGrid<Engine_t> &gr, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(nsp, inspin, IntC), apply_(true)
		{
			grads_=&gr;
		}

		ListBlochParams(const BlochParams<BPops, offset_T>  &one, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(one, IntC),
			grads_(NULL), apply_(true)
		{}

		ListBlochParams(int nsp, const BlochParams<BPops, offset_T>  &dup, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(nsp,dup, IntC),
			grads_(NULL), apply_(true)
		{}

		ListBlochParams( const ListBlochParams<NullGrid, BPops, offset_T> &dup, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(dup, IntC),
			grads_(NULL), apply_(true)
		{}

		ListBlochParams( const ListBlochParams<NullGrid, BPops, offset_T> &dup,GradientGrid<Engine_t> &gr, int IntC=InitialCondition::AllUp):
			ListBlochParams<NullGrid, BPops, offset_T>(dup, IntC),
			grads_(&gr), apply_(true)
		{}

		ListBlochParams( const ListBlochParams<Engine_t, BPops, offset_T> &dup,GradientGrid<Engine_t> &gr):
			ListBlochParams<NullGrid, BPops, offset_T>(dup, IntC),
			grads_(&gr), apply_(true)
		{}

		ListBlochParams( const ListBlochParams<GradientGrid<Engine_t>, BPops, offset_T > &dup):
			ListBlochParams<NullGrid, BPops, offset_T>(dup),
			grads_(dup.grads_), apply_(true)
		{}

		~ListBlochParams(){ grads_=NULL;	}

		iterator begin(){	return iterator(*this);	}

		void SetGrads(GradientGrid<Engine_t> &in);

		void SetGrid(Engine_t &in);

		inline offset_T offset(int i);
		inline offset_T offset(int i, double t); //time dependant grids


		inline const GradientGrid<Engine_t> *Grads() const;
		inline GradientGrid<Engine_t> *Grads();

		inline Engine_t *Grid();
		inline const Engine_t *Grid() const;

		inline void off();
		inline void on();
		inline void GradOff();
		inline void GradOn();
		inline bool apply();

		inline rmatrix Jacobian(int i);
		inline void Jacobian(rmatrix &out, int i);


		void operator=(const ListBlochParams &rhs);

		void operator=(const ListBlochParams<NullGrid, BPops, offset_T> &rhs);

		void operator=(const ListBlochParams<Engine_t, BPops, offset_T> &rhs);
};


template<class GridEngine_t, int BPops, class offset_T>
std::ostream& operator<<(std::ostream &oo,  ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > &out);


/************ ITTERATORS ********************/
template<int BPops, class offset_T>
class ListBlochParamsIterator<NullGrid, BPops, offset_T>
{
	protected:
		BlochParams< BPops, offset_T>  *CurBlochP_;	//this acts like the current parameter...it speed things up alot
										//foreach itteration point, there are numerous acsess to the same part in
										//the list...this reduces it to one...it drops seconds of speed depending on the list size

		ListBlochParams<NullGrid, BPops, offset_T> *mylist_; //a ptr to the List we wish to iterate
		bool notended_;
		int curpos_;

		Range R_;

		const void IterErr();
	public:

		typedef BlochParams<BPops, offset_T> Parameter_T;

		ListBlochParamsIterator():
			CurBlochP_(NULL), mylist_(NULL), notended_(false), curpos_(0), R_(0,0)
		{}

		ListBlochParamsIterator(ListBlochParams<NullGrid, BPops, offset_T> &in):
			mylist_(&in),notended_(true), curpos_(0), R_(0,in.size(), 1)
		{
			CurBlochP_=&(*mylist_)(curpos_);
		}

	//RANGE ctors --> an important piece in the MPI ability
		ListBlochParamsIterator(ListBlochParams<NullGrid, BPops, offset_T> &in, Range R) :
			mylist_(&in),notended_(true), curpos_(R.first(0)), R_(R)
		{
			CurBlochP_=&(*mylist_)(curpos_);
		}

		ListBlochParamsIterator(const ListBlochParamsIterator &in) :
			CurBlochP_(in.CurBlochP_),
			mylist_(in.mylist_),
			notended_(in.notended_),
			curpos_(in.curpos_),
			R_(in.R_)
		{}

	//RANGE ctors --> an important piece in the MPI ability
		ListBlochParamsIterator(const ListBlochParamsIterator &in, Range R):
			CurBlochP_(in.CurBlochP_),
			mylist_(in.mylist_),
			notended_(in.notended_),
			curpos_(R.first(in.curpos_)),
			R_(R)
		{}


		~ListBlochParamsIterator()
		{
			CurBlochP_=NULL;
			mylist_=NULL;
		}

		inline void setRange(Range R){	R_=R;	}

		inline double gamma() const;
		inline double gamma(int i) const;
		inline double gammaGauss() const;
		inline double gammaGauss(int i) const;

		inline double &temperature();
		inline double temperature() const;
		inline double &temperature(int i);
		inline double temperature(int i) const;
		inline void temperature(double intemp);

		inline std::string symbol() const;

		inline offset_T &Mo();
		inline offset_T Mo() const;
		inline offset_T &Mo(int i);
		inline offset_T Mo(int i) const;
		inline void Mo(offset_T ne);

		inline offset_T &Bo();
		inline offset_T Bo() const;
		inline offset_T &Bo(int i);
		inline offset_T Bo(int i) const;
		inline void Bo(offset_T ne);

		inline offset_T &omega();
		inline offset_T omega() const;
		inline offset_T omega(int i) const;
		inline offset_T &omega(int i);
		//inline void omega(offset_T ne) const;

		void setInitialCondition(int IC, int Halfflag=0);

		inline SolveUnits_ SolveUnits() const;

		inline std::string symbol(int i) const;

		inline double &moles();
		inline double moles() const;
		inline double &moles(int i);
		inline double moles(int i) const;
		inline void moles(double ne);


		inline double &TotMo();
		inline double &totalMo();
		inline double meanMo() const;

		inline void offset(double ne);

		inline rmatrix Jacobian();
		inline void Jacobian(rmatrix &out);
		inline rmatrix Jacobian(int i);
		inline void Jacobian(rmatrix &out, int i);

		void calcMo();


		inline BlochParams<BPops, offset_T>  &operator()();

		inline void operator++();

		inline void operator++(int);


		void operator--();

		void reset();

		inline operator bool();

		inline int Position()const{	return curpos_;	}
		inline int curpos()const{	return curpos_;	}

		inline int size() const {	return mylist_->size();	}

		BlochParams<BPops, offset_T>  &CurParam();
};

/************ Simple GRID ITTERATORS ********************/

//
/// Itterators for the List of bloch params  WITH GRIDS attached
// --> Inherts all that from the Grid iterator and NullGRid List Bloch iterator (above0)
template<class GridEngine_t, int BPops, class offset_T>
class ListBlochParamsIterator  :
	public ListBlochParamsIterator<NullGrid, BPops, offset_T>,
	public GridEngine_t::iterator

{
	protected:
		ListBlochParams<GridEngine_t, BPops, offset_T > *mylist_;

	public:
		typedef typename GridEngine_t::iterator Griditerator;
	//ctors
		ListBlochParamsIterator():
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(),
			GridEngine_t::iterator(),
			mylist_(0)
		{}

		ListBlochParamsIterator(ListBlochParams<GridEngine_t, BPops, offset_T > &in):
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(in),
			GridEngine_t::iterator(in.Grid()),
			mylist_(&in)
		{}

		ListBlochParamsIterator(const ListBlochParamsIterator &in) :
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(in),
			GridEngine_t::iterator(in),
			mylist_(in.mylist_)
		{};

	//RANGE ctors--> An important element in the MPI ability of this lib
		ListBlochParamsIterator(ListBlochParams<GridEngine_t , BPops, offset_T> &in, Range R) :
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(in,R),
			GridEngine_t::iterator(in.Grid(), R),
			mylist_(&in)
		{};

		ListBlochParamsIterator( ListBlochParamsIterator< GridEngine_t, BPops, offset_T > &in, Range R):
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(in,R),
			GridEngine_t::iterator(in,R),
			mylist_(in.mylist_)
		{};


		~ListBlochParamsIterator()
		{
			mylist_=NULL;
		}

	//	inline coord<> &Point();

		inline void operator++();

		inline void operator++(int);

		inline operator bool();
		inline int curpos()const{	return curpos_;	}
		void reset();

};


/************ Gradient Grid ITTERATORS ********************/

//
/// Itterators for the List of bloch params  WITH GRIDS attached
//  specifically GRADIENTs (and Grids)
// --> Inherts all that from the Gradient Grid iterator and NullGRid List Bloch iterator (above0)
template<class GridEngine_t, int BPops, class offset_T>
class ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T > :
	public ListBlochParamsIterator<NullGrid, BPops, offset_T>,
	public GradientGridIter<GridEngine_t>

{
	protected:
		ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > *mylist_;

	public:

	//ctors
		ListBlochParamsIterator():
			mylist_(NULL)
		{}

		ListBlochParamsIterator(ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > &in) :
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(in),
			GradientGridIter<GridEngine_t>(in.Grads()),
			mylist_(&in)
		{
			CurBlochP_=&(*mylist_)(0);
		}

		ListBlochParamsIterator(const ListBlochParamsIterator &in) :
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(in),
			GradientGridIter<GridEngine_t>(in),
			mylist_(in.mylist_)
		{}

	//RANGE ctors--> An important element in the MPI ability of this lib
		ListBlochParamsIterator(ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > &in, Range R) :
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(in,R),
			GradientGridIter<GridEngine_t>(in, R)
		{}

		ListBlochParamsIterator( ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T > &in, Range R) :
			ListBlochParamsIterator<NullGrid, BPops, offset_T>(in,R),
			GradientGridIter<GridEngine_t>(in, R)
		{}



		~ListBlochParamsIterator()
		{
			mylist_=NULL;
		}

//		inline rmatrix Jacobian();
//		inline void Jacobian(rmatrix &out);
//		inline rmatrix Jacobian(int i);
//		inline void Jacobian(rmatrix &out, int i);

		inline offset_T offset();
		inline offset_T offset(int i);

		inline offset_T offset(double t); //time dependnat
		inline offset_T offset(int i, double t);

		inline int curpos() const {	return ListBlochParamsIterator<NullGrid, BPops, offset_T>::curpos();	}
		inline void operator++();

		inline void operator++(int);

		inline operator bool();
		void reset();

};

END_BL_NAMESPACE


#include "bloch/listblochpars_meth.h"

#endif

