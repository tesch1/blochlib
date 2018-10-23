/* listblochpars_meth.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-8-01
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
 	listblochpars_meth.h-->a large list of BlochParams<BPops, offset_T> <>  with Grids attached...
 */


#ifndef _listlochparas_meth_h_
#define _listlochparas_meth_h_ 1


/* maintains a list of bloch parameters
	there are two specializations,
		one with no grid (a 'NullGrid') SEE LISTBLOCHPARS_METH.CC
		and one with a GradientGrid (a grid snapped with a graident vector)
	all others inlcude simple grids.

*/

//using namespace std;


BEGIN_BL_NAMESPACE


template<int BPops, class offset_T >
double ListBlochParams<NullGrid, BPops, offset_T>::FindMaxMo()
{
	double tmp=0.;
	iterator myIt((*this));
	while(myIt){
		myIt.calcMo();
		if(abs(myIt.Mo())>abs(tmp)) tmp=myIt.Mo();
		++myIt;
	}
	return tmp;
}


template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::reCalcMo()
{
	iterator myIt((*this));
	while(myIt){
		myIt.calcMo();
		myIt.setInitialCondition(InitCond);
		++myIt;
	}
	//this is the chunk for dimmensionless integrations
	//we need to know the totalMo BEFORE normalization
	//becuase we want the total sum to be 1 so each individual Spin magnetization
	//would be its normal Mo/totMo

	if(MySolveUnits_!=Normal)
	{
		calcTotalMo();
		myIt.reset();
		double tmtotmo=0;
		while(myIt){
			myIt.Mo(myIt.Mo()/totMo_);
			tmtotmo+=norm(myIt.Mo());
			++myIt;
		}
		totMo_=tmtotmo;
	}
}


template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::calcTotalMo()
{
	ListBlochParams<NullGrid, BPops, offset_T>::iterator myIt(*this);
	totMo_=0;
	while(myIt){
		myIt.calcMo();
		myIt.setInitialCondition(InitCond);
		totMo_+=norm(myIt.Mo());
		++myIt;
	}
}


template<int BPops, class offset_T >
SolveUnits_ &ListBlochParams<NullGrid, BPops, offset_T>::SolveUnits(){	return MySolveUnits_;	}

template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::SetSolveUnits(const SolveUnits_ in)
{
	//do not do recalculate anything if we do not have to
	if(MySolveUnits_!=in){
		MySolveUnits_=in;
		reCalcMo();
	}
}

template<int BPops, class offset_T >
typename ListBlochParams<NullGrid, BPops, offset_T>::iterator  ListBlochParams<NullGrid, BPops, offset_T>::begin()
{	return iterator(*this);	}



template<int BPops, class offset_T >
BlochParams<BPops, offset_T>  &ListBlochParams<NullGrid, BPops, offset_T>::operator()(int i)
{ return Vpars_(i);	}

template<int BPops, class offset_T >
BlochParams<BPops, offset_T>  &ListBlochParams<NullGrid, BPops, offset_T>::operator[](int i)
{	return Vpars_(i);	}

template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::operator=(const ListBlochParams<NullGrid, BPops, offset_T> &rhs)
{
	if(&rhs==this) return;
	Vpars_=rhs.Vpars_;
}


template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::operator=(const BlochParams<BPops, offset_T>  &rhs)
{
	Vpars_.resize(1);
	Vpars_(0)=rhs;
}

template<int BPops, class offset_T >
ListBlochParams<NullGrid, BPops, offset_T> &ListBlochParams<NullGrid, BPops, offset_T>::operator+(ListBlochParams &rhs)
{
	for(int i=0;i<rhs.Vpars_.size();i++){
		Vpars_.push_back(rhs(i));
	}
	if(MySolveUnits_==Dimensionless) reCalcMo();
	else calcTotalMo();
	return *this;
}

template<int BPops, class offset_T >
ListBlochParams<NullGrid, BPops, offset_T> &ListBlochParams<NullGrid, BPops, offset_T>::operator+(const BlochParams<BPops, offset_T>  &rhs)
{
	Vpars_.push_back(rhs);
	if(MySolveUnits_==Dimensionless) reCalcMo();
	else calcTotalMo();

	return *this;
}

template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::operator+=( ListBlochParams &rhs)
{
	for(int i=0;i<rhs.Vpars_.size();i++){
		Vpars_.push_back(rhs(i));
	}
	if(MySolveUnits_==Dimensionless) reCalcMo();
	else calcTotalMo();
}

template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::operator+=(const BlochParams<BPops, offset_T> &rhs)
{
	Vpars_.push_back(rhs);
	if(MySolveUnits_==Dimensionless) reCalcMo();
	else calcTotalMo();
}


template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::push_back(const BlochParams<BPops, offset_T>  &rhs)
{
	Vpars_.push_back(rhs);
	if(MySolveUnits_==Dimensionless) reCalcMo();
	else calcTotalMo();
}


template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::print(std::ostream &oo)
{
	oo<<"Base Bfield="<<Bo()<<" (note::not a valid representation Bo() is different everywhere)"<<std::endl;
	oo<<"Base temperature="<<temperature()<<" (note::not a valid representation temp() is different everywhere)"<<std::endl;

	oo<<"Initial Condition:: "<<InitialCondition::name(InitCond)<<std::endl;

	if(MySolveUnits_!=Normal){
		oo<<"'Mo' for each spin parameter is calculated assuming dimmensionless units"<<std::endl;
		oo<<" so there grand total Mo will add up to 1"<<std::endl;
	}
//	out.Vpars_.RecordSep="\n";
	iterator myIt((*this));
	while(myIt){
		oo<<myIt.CurParam()<<std::endl;
		++myIt;
	}
}

template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::print(){ print(std::cout);	}	//text print
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


template<int BPops, class offset_T >
bool ListBlochParams<NullGrid, BPops, offset_T>::write(std::fstream &oo)	//binary write
{
	if(!oo){
		std::cerr<<std::endl<<"Error: ListBlochPars::write(fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return false;
	}
	oo<<std::endl<<"LISTBLOCHPARS\n";
	long totalsize=0;
	iterator myIt((*this));
	while(myIt){
		totalsize+=myIt.CurParam().binarySize();
		++myIt;
	}
	oo.write((char *)&totalsize, sizeof(long));
	oo.write((char *)&Vpars_(0).Bo_, sizeof(double));
	oo.write((char *)&Vpars_(0).temperature_, sizeof(double));
	if(MySolveUnits_==Normal){
		int normal=1;
		oo.write((char *)&normal, sizeof(int));
	}else{
		int normal=0;
		oo.write((char *)&normal, sizeof(int));
	}
	int nsp=Vpars_.size();
	oo.write((char *)&nsp, sizeof(int));
	myIt.reset();
	while(myIt){
		if(myIt.CurParam().write(oo)) return false;
		++(myIt);
	}
	oo<<std::endl<<"ENDLISTBLOCHPARS\n";
	return true;

}

//template<>
template<int BPops, class offset_T >
bool ListBlochParams<NullGrid, BPops, offset_T>::read(std::fstream &in)	//binart read
{
	if(!in){
		std::cerr<<std::endl<<"Error: ListBlochPars::read(fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return false;
	}
	std::string temp="";
	char tm[1000];
	while(temp!="LISTBLOCHPARS" && !in.eof()){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
	}
	if(in.eof())
	{
		std::cerr<<std::endl<<"Error: ListBlochPars::read(fstream"<<std::endl;
		std::cerr<<" No 'ListBlochParameters' found...."<<std::endl;
		return 0;
	}
	long totalsize=0;
	in.read((char *)&totalsize, sizeof(long));
	in.read((char *)&Vpars_(0).Bo_, sizeof(double));
	in.read((char *)&Vpars_(0).temperature_, sizeof(double));

	int normal=1;
	in.read((char *)&normal, sizeof(int));
	if(normal==1) MySolveUnits_=Normal;
	else	MySolveUnits_=Dimensionless;

	int nsp=1;
	in.read((char *)&nsp, sizeof(int));
	Vpars_.resize(nsp);
	iterator myIt((*this));
	while(myIt){
		myIt.CurParam().read(in);
		//stepbin+=Vpars_(curpos_).binarySize();
		++myIt;
	}
	while(temp!="ENDLISTBLOCHPARS"){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
	}
	return true;
}

template<int BPops, class offset_T >
std::ostream& operator<<(std::ostream &oo, ListBlochParams<NullGrid, BPops, offset_T> &out)
{
	out.print(oo);
	return oo;
}



/**************************** INLINE METHODS FOR "NULLGRID" ***/


/*** Symbol ***/
template<int BPops, class offset_T >
inline std::string ListBlochParams<NullGrid, BPops, offset_T>::symbol(int i) const
{	return Vpars_[i].symbol();	}

/*** temperature **/

template<int BPops, class offset_T >
inline double &ListBlochParams<NullGrid, BPops, offset_T>::temperature(int i)
{	return Vpars_[i].BasicBlochParams<offset_T>::temperature();	}

template<int BPops, class offset_T >
inline double &ListBlochParams<NullGrid, BPops, offset_T>::temperature()
{	return Vpars_[0].BasicBlochParams<offset_T>::temperature();	}

template<int BPops, class offset_T >
inline double ListBlochParams<NullGrid, BPops, offset_T>::temperature() const
{	return Vpars_[0].BasicBlochParams<offset_T>::temperature();	}

template<int BPops, class offset_T >
inline double ListBlochParams<NullGrid, BPops, offset_T>::temperature(int i) const
{	return Vpars_[i].BasicBlochParams<offset_T>::temperature();	}

template<int BPops, class offset_T >
inline void ListBlochParams<NullGrid, BPops, offset_T>::temperature(double Tempin)
{	Vpars_(0).temperature_=Tempin; reCalcMo();	}



/*** Bo ***/
template<int BPops, class offset_T >
inline offset_T &ListBlochParams<NullGrid, BPops, offset_T>::Bo(int i)
{	return Vpars_[i].BasicBlochParams<offset_T>::Bo();	}

template<int BPops, class offset_T >
inline offset_T ListBlochParams<NullGrid, BPops, offset_T>::Bo(int i) const
{	return Vpars_[i].BasicBlochParams<offset_T>::Bo();	}

template<int BPops, class offset_T >
inline offset_T &ListBlochParams<NullGrid, BPops, offset_T>::Bo()
{	return Vpars_[0].BasicBlochParams<offset_T>::Bo();	}

template<int BPops, class offset_T >
inline offset_T ListBlochParams<NullGrid, BPops, offset_T>::Bo() const
{	return Vpars_[0].BasicBlochParams<offset_T>::Bo();	}

template<int BPops, class offset_T >
void ListBlochParams<NullGrid, BPops, offset_T>::Bo(offset_T inBo)
{	Vpars_(0).Bo_=inBo; Vpars_(0).Ho_=Vpars_(0).Bo_*1/(permVac);	reCalcMo();	}



/*** Mo ***/

template<int BPops, class offset_T >
inline offset_T &ListBlochParams<NullGrid, BPops, offset_T>::Mo(int i)
{	return Vpars_[i].BasicBlochParams<offset_T>::Mo();	}

template<int BPops, class offset_T >
inline offset_T ListBlochParams<NullGrid, BPops, offset_T>::Mo(int i) const
{	return Vpars_[i].BasicBlochParams<offset_T>::Mo();	}

template<int BPops, class offset_T >
inline double &ListBlochParams<NullGrid, BPops, offset_T>::TotMo()
{	return totMo_;	}

template<int BPops, class offset_T >
inline double ListBlochParams<NullGrid, BPops, offset_T>::TotMo() const
{	return totMo_;	}

template<int BPops, class offset_T >
inline double &ListBlochParams<NullGrid, BPops, offset_T>::totalMo()
{	return totMo_;	}

template<int BPops, class offset_T >
inline double ListBlochParams<NullGrid, BPops, offset_T>::totalMo() const
{	return totMo_;	}

template<int BPops, class offset_T >
inline double ListBlochParams<NullGrid, BPops, offset_T>::meanMo() const
{	return totMo_/double(size());	}


/**** moles ***/
template<int BPops, class offset_T >
inline double &ListBlochParams<NullGrid, BPops, offset_T>::moles(int i)
{	return Vpars_[i].BasicBlochParams<offset_T>::moles();	}

template<int BPops, class offset_T >
inline double ListBlochParams<NullGrid, BPops, offset_T>::moles(int i) const
{	return Vpars_[i].BasicBlochParams<offset_T>::moles();	}

/** gamma ***/
template<int BPops, class offset_T >
inline double ListBlochParams<NullGrid, BPops, offset_T>::gamma(int i) const
{	return Vpars_[i].gamma();	}

template<int BPops, class offset_T >
inline double ListBlochParams<NullGrid, BPops, offset_T>::gammaGauss(int i) const
{	return Vpars_[i].gammaGauss();	}

/** omgea ***/
template<int BPops, class offset_T >
inline offset_T &ListBlochParams<NullGrid, BPops, offset_T>::omega(int i)
{	return Vpars_[i].omega();	}

template<int BPops, class offset_T >
inline offset_T ListBlochParams<NullGrid, BPops, offset_T>::omega(int i) const
{	return Vpars_[i].omega();	}


/******************** END ****************/


/****************** GRID PARAMETER BITS ********************/



// Gradient Grid compatability...
template<class GridEngine_t, int BPops, class offset_T>
void	ListBlochParams<GridEngine_t, BPops, offset_T>::operator=
	(const ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > &rhs)
{
	ListBlochParams<NullGrid, BPops, offset_T>::operator=(rhs);
	grads_=rhs.Grid();
}

//self equal operator
template<class GridEngine_t, int BPops, class offset_T>
void	ListBlochParams<GridEngine_t, BPops, offset_T>::operator=(const ListBlochParams<GridEngine_t, BPops, offset_T> &rhs)
{
	if(this==&rhs) return;
	ListBlochParams<NullGrid, BPops, offset_T>::operator=(rhs);
	grads_=rhs.grads_;
}

// equal operator with a "NullGrid"
template<class GridEngine_t, int BPops, class offset_T>
void	ListBlochParams<GridEngine_t, BPops, offset_T>::operator=(const ListBlochParams<NullGrid, BPops, offset_T> &rhs)
{
	ListBlochParams<NullGrid, BPops, offset_T>::operator=(rhs);
	grads_=0;
}


template<class GridEngine_t, int BPops, class offset_T>
inline const GridEngine_t *ListBlochParams<GridEngine_t, BPops, offset_T>::Grid()	 const
{
	return grads_;
}


template<class GridEngine_t, int BPops, class offset_T>
inline GridEngine_t *ListBlochParams<GridEngine_t, BPops, offset_T>::Grid()
{
	return grads_;
}

template<class GridEngine_t, int BPops, class offset_T>
std::ostream& operator<<(std::ostream &oo, ListBlochParams<GridEngine_t, BPops, offset_T > &out)
{
	out.print(oo);
	return oo;
}

/******************8 GRADIENT PARAMETER LISTS **********************/

/* the method for using the parameters with Gradients...
uses the above 'NullGrid' class as its base

*/



template<class GridEngine_t, int BPops, class offset_T>
void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::SetGrads(GradientGrid<GridEngine_t> &in) 
{
	if(in.size() != size()){
		BLEXCEPTION(" Gradient grid size must be the same as the number of spins...")
	}
	grads_=&in;
}

template<class GridEngine_t, int BPops, class offset_T>
void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::SetGrid(GridEngine_t &in)
{
	GradientGrid<GridEngine_t>::SetGrid(in);
}

template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::off()
{	apply_=false;	}

template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::on()
{		apply_=true;	}

template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::GradOff()
{	apply_=false;	}

template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::GradOn()
{		apply_=true;	}

template<class GridEngine_t, int BPops, class offset_T>
inline GradientGrid<GridEngine_t> *ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::Grads()
{ return grads_;	}

template<class GridEngine_t, int BPops, class offset_T>
inline const GradientGrid<GridEngine_t> *ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::Grads()const
{ return grads_;	}

template<class GridEngine_t, int BPops, class offset_T>
inline offset_T ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::offset(int i)
{ return gammaGauss(i)*grads_->offset(i);	}

template<class GridEngine_t, int BPops, class offset_T>
inline offset_T ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::offset(int i, double t)
{ return gammaGauss(i)*grads_->offset(i,t);	}

template<class GridEngine_t, int BPops, class offset_T>
inline GridEngine_t *ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::Grid()
{ return grads_->Grid();	}

template<class GridEngine_t, int BPops, class offset_T>
inline const GridEngine_t *ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::Grid()const
{ return grads_->Grid();	}


template<class GridEngine_t, int BPops, class offset_T>
inline bool ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::apply()
{	return apply_;	}


template<class GridEngine_t, int BPops, class offset_T>
void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::operator=
	(const ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > &rhs)
{
	if(this==&rhs) return;
	grads_=rhs.grads_;
	ListBlochParams::operator=(rhs);
}

template<class GridEngine_t, int BPops, class offset_T>
void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::operator=
	(const ListBlochParams<NullGrid, BPops, offset_T> &rhs)
{
	grads_=0;
	ListBlochParams::operator=(rhs);
}

template<class GridEngine_t, int BPops, class offset_T>
void ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T >::operator=
	(const ListBlochParams<GridEngine_t,BPops, offset_T> &rhs)
{
	grads_=0;
	ListBlochParams::operator=(rhs);
}

template<class GridEngine_t, int BPops, class offset_T>
std::ostream& operator<<(std::ostream &oo, ListBlochParams<GradientGrid<GridEngine_t> > &out)
{
	out.print(oo);
	return oo;
}

/************ ITERATORS ********************/

/*
template<int BPops, class offset_T>
ListBlochParamsIterator<NullGrid, BPops, offset_T>::ListBlochParamsIterator() :
	CurBlochP_(NULL), mylist_(NULL), notended_(false), curpos_(0), R_(0,0)
{}

template<int BPops, class offset_T>
ListBlochParamsIterator<NullGrid, BPops, offset_T>::ListBlochParamsIterator( ListBlochParams<NullGrid, BPops, offset_T> &in) :
	mylist_(&in),notended_(true), curpos_(0), R_(0,in.size(), 1)

{
	CurBlochP_=&(*mylist_)(0);
}

template<int BPops, class offset_T>
ListBlochParamsIterator<NullGrid, BPops, offset_T>::ListBlochParamsIterator( ListBlochParams<NullGrid, BPops, offset_T> &in, Range R) :
	mylist_(&in),notended_(true), curpos_(R.first(0)), R_(R)
{
	CurBlochP_=&(*mylist_)(0);
}



template<int BPops, class offset_T>
ListBlochParamsIterator<NullGrid, BPops, offset_T>::ListBlochParamsIterator(const ListBlochParamsIterator<NullGrid, BPops, offset_T> &in) :
	CurBlochP_(in.CurBlochP_), mylist_(in.mylist_), notended_(in.notended_), curpos_(in.curpos_), R_(in.R_)
{}

template<int BPops, class offset_T>
ListBlochParamsIterator<NullGrid, BPops, offset_T>::ListBlochParamsIterator(const ListBlochParamsIterator<NullGrid, BPops, offset_T> &in, Range R) :
	CurBlochP_(in.CurBlochP_), mylist_(in.mylist_), notended_(in.notended_), curpos_(in.curpos_), R_(R)
{}
*/

template<int BPops, class offset_T>
const void  ListBlochParamsIterator<NullGrid, BPops, offset_T>::IterErr() 
{
	BLEXCEPTION(" iterated too far...")
}



template<int BPops, class offset_T>
void ListBlochParamsIterator<NullGrid, BPops, offset_T>::calcMo()
{	CurBlochP_->calcMo();}


template<int BPops, class offset_T>
void ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator--()
{
	if(!curpos_) IterErr();
	curpos_--;
	CurBlochP_=&(*mylist_)(curpos_);
}

template<int BPops, class offset_T>
void ListBlochParamsIterator<NullGrid, BPops, offset_T>::reset()
{
	notended_=true;
	curpos_=0;
	CurBlochP_=&(*mylist_)(0);
}


template<int BPops, class offset_T>
BlochParams<BPops, offset_T>  &ListBlochParamsIterator<NullGrid, BPops, offset_T>::CurParam()
{	return *CurBlochP_;	}

/** inline methods for NullGrid itterator */

template<int BPops, class offset_T>
inline SolveUnits_ ListBlochParamsIterator<NullGrid, BPops, offset_T>::SolveUnits() const
{	return mylist_->MySolveUnits_;	}


/****** Temperature *****/
template<int BPops, class offset_T>
inline double &ListBlochParamsIterator<NullGrid, BPops, offset_T>::temperature()
{	return CurBlochP_->temperature_;	}

template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::temperature() const
{	return CurBlochP_->temperature_;	}

template<int BPops, class offset_T>
inline double &ListBlochParamsIterator<NullGrid, BPops, offset_T>::temperature(int i)
{	return mylist_->temperature(i);	}

template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::temperature(int i) const
{	return mylist_->temperature(i);	}

template<int BPops, class offset_T>
inline void ListBlochParamsIterator<NullGrid, BPops, offset_T>::temperature(double ne)
{	mylist_->temperature(ne);	}

/**** symbol ****/

template<int BPops, class offset_T>
inline std::string ListBlochParamsIterator<NullGrid, BPops, offset_T>::symbol(int i) const
{	return mylist_->symbol(i);	}

template<int BPops, class offset_T>
inline std::string ListBlochParamsIterator<NullGrid, BPops, offset_T>::symbol() const
{	return CurBlochP_->symbol();	}

/** moles ***/

template<int BPops, class offset_T>
inline double &ListBlochParamsIterator<NullGrid, BPops, offset_T>::moles()
{	return CurBlochP_->BasicBlochParams<offset_T>::moles();	}

template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::moles() const
{	return CurBlochP_->BasicBlochParams<offset_T>::moles();	}

template<int BPops, class offset_T>
inline double &ListBlochParamsIterator<NullGrid, BPops, offset_T>::moles(int i)
{	return mylist_->moles(i);	}

template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::moles(int i) const
{	return mylist_->moles(i);	}

template<int BPops, class offset_T>
inline void ListBlochParamsIterator<NullGrid, BPops, offset_T>::moles(double ne)
{	CurBlochP_->moles(ne);	}

/**** gamma ***/
template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::gamma() const
{	return CurBlochP_->gamma();	}

template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::gammaGauss() const
{	return CurBlochP_->gammaGauss();	}

template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::gamma(int i) const
{	return mylist_->gamma(i);	}

template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::gammaGauss(int i) const
{	return mylist_->gammaGauss(i);	}

/*** omega *****/
template<int BPops, class offset_T>
inline offset_T &ListBlochParamsIterator<NullGrid, BPops, offset_T>::omega()
{	return CurBlochP_->omega();	}


template<int BPops, class offset_T>
inline offset_T &ListBlochParamsIterator<NullGrid, BPops, offset_T>::omega(int i)
{	return mylist_->omega(i);	}

template<int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<NullGrid, BPops, offset_T>::omega() const
{	return CurBlochP_->omega();	}


template<int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<NullGrid, BPops, offset_T>::omega(int i) const
{	return mylist_->omega(i);	}



/********** Mo *********/
template<int BPops, class offset_T>
inline offset_T &ListBlochParamsIterator<NullGrid, BPops, offset_T>::Mo(int i)
{	return mylist_->Mo(i);	}

template<int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<NullGrid, BPops, offset_T>::Mo(int i) const
{	return mylist_->Mo(i);	}

template<int BPops, class offset_T>
inline void ListBlochParamsIterator<NullGrid, BPops, offset_T>::Mo(offset_T ne)
{	CurBlochP_->Mo(ne);	}

template<int BPops, class offset_T>
inline offset_T &ListBlochParamsIterator<NullGrid, BPops, offset_T>::Mo()
{	return CurBlochP_->BasicBlochParams<offset_T>::Mo();	}

template<int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<NullGrid, BPops, offset_T>::Mo() const
{	return CurBlochP_->BasicBlochParams<offset_T>::Mo();	}


template<int BPops, class offset_T>
inline double &ListBlochParamsIterator<NullGrid, BPops, offset_T>::TotMo()
{	return mylist_->TotMo();	}

template<int BPops, class offset_T>
inline double &ListBlochParamsIterator<NullGrid, BPops, offset_T>::totalMo()
{	return mylist_->TotMo();	}

template<int BPops, class offset_T>
inline double ListBlochParamsIterator<NullGrid, BPops, offset_T>::meanMo() const
{	return mylist_->TotMo()/double(mylist_->size());	}

/********** Bo *********/

template<int BPops, class offset_T>
inline offset_T &ListBlochParamsIterator<NullGrid, BPops, offset_T>::Bo()
{	return CurBlochP_->Bo_;	}

template<int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<NullGrid, BPops, offset_T>::Bo() const
{	return CurBlochP_->Bo_;	}

template<int BPops, class offset_T>
inline offset_T &ListBlochParamsIterator<NullGrid, BPops, offset_T>::Bo(int i)
{	return mylist_->Bo(i);	}

template<int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<NullGrid, BPops, offset_T>::Bo(int i) const
{	return mylist_->Bo(i);	}

template<int BPops, class offset_T>
inline void ListBlochParamsIterator<NullGrid, BPops, offset_T>::Bo(offset_T ne)
{	return CurBlochP_->Bo(ne);	}

template<int BPops, class offset_T>
  void ListBlochParamsIterator<NullGrid, BPops, offset_T>::
setInitialCondition(int IC, int Halfflag)
{
	CurBlochP_->setInitialCondition(IC, Halfflag);
}


template<int BPops, class offset_T>
inline BlochParams<BPops, offset_T>  &ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator()()
{
	return *CurBlochP_;
}

template<int BPops, class offset_T>
inline void ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator++()
{
//	std::cout<<"PRE:curpose: "<<curpos_<<" Range: "<<R_<<" R_.last(mylist_->size()):"<<R_.last(mylist_->size())<<std::endl;
	if(!notended_) IterErr();
	curpos_+=R_.stride();
//	std::cout<<"POST:curpose: "<<curpos_<<" Range: "<<R_<<" R_.last(mylist_->size()):"<<R_.last(mylist_->size())<<std::endl;
	if(curpos_>=R_.last(mylist_->size())){ notended_=false; return ;	}
	CurBlochP_=&(*mylist_)(curpos_);

}

template<int BPops, class offset_T>
inline void ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator++(int){ operator++();	}

template<int BPops, class offset_T>
inline ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator bool(){ return notended_;	}

/******************************* GRIDS ITERTATORS *******************/

//The Advance iterator..simply passes it off to the List Blach <NULL> and the Gradient Grid itterators
//as they better have the same number of elements

//this is the "++myIter" (faster)
template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParamsIterator< GridEngine_t, BPops, offset_T >::operator++()
{
	ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator++() ;
	Griditerator::operator++();
}

//this is the "myIter++" (a bit slower, but much easier to type and think about)
template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParamsIterator< GridEngine_t, BPops, offset_T >::operator++(int){ 	operator++();	}


//an operator the function as the "have we hit the end yet?"
// an examples....>  while(myIter){...}
template<class GridEngine_t, int BPops, class offset_T>
inline ListBlochParamsIterator< GridEngine_t, BPops, offset_T >::operator bool()
{
	return (ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator bool()) || ( Griditerator::operator bool() );
}

//reset the iterator after you have used it once (typically not used,
// but is it on occasion and thus saves a bit of initalization time)
template<class GridEngine_t, int BPops, class offset_T>
void ListBlochParamsIterator< GridEngine_t, BPops, offset_T >::reset()
{
	Griditerator::reset();
	ListBlochParamsIterator<NullGrid, BPops, offset_T>::reset();
}

/*template<class GridEngine_t, int BPops, class offset_T>
inline coord<> &ListBlochParamsIterator< GridEngine_t, BPops, offset_T >::Point()
{
	return Point();
}*/


/************* ITTERATOR FOR 'Bloch Params WITH Gradient GRIDS!! ************/
/*
//null Ctor
template<class GridEngine_t, int BPops, class offset_T>
ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::ListBlochParamsIterator() :
	mylist_(NULL)
{}


//the most use ctor..make an iterator from the List of params at hand
template<class GridEngine_t, int BPops, class offset_T>
ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::ListBlochParamsIterator(ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > &in) :
	ListBlochParamsIterator<NullGrid, BPops, offset_T>(in),
	GradientGridIter<GridEngine_t>(in.Grads()),
	mylist_(&in)
{
	CurBlochP_=&(*mylist_)(0);
}


template<class GridEngine_t, int BPops, class offset_T>
ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::ListBlochParamsIterator(const ListBlochParamsIterator &in) :
	ListBlochParamsIterator<NullGrid, BPops, offset_T>(in),
	GradientGridIter<GridEngine_t>(in),
	mylist_(in.mylist_)
{}

template<class GridEngine_t, int BPops, class offset_T>
ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::ListBlochParamsIterator( ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T > &in, Range R) :
	ListBlochParamsIterator<NullGrid, BPops, offset_T>(in,R),
	GradientGridIter<GridEngine_t>(in, R)
{}


template<class GridEngine_t, int BPops, class offset_T>
ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::ListBlochParamsIterator( ListBlochParams<GradientGrid<GridEngine_t>, BPops, offset_T > &in, Range R) :
	ListBlochParamsIterator<NullGrid, BPops, offset_T>(in,R),
	GradientGridIter<GridEngine_t>(in, R)
{}
*/
/*
//the offset operators are the main change from the ListBlochParamsIter<NullGrid>
// class  beucase it is assumed we have a gradient, the offset changes given
// any particular point in space
template<class GridEngine_t, int BPops, class offset_T>
inline double ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::offset()
{
	if(mylist_->Grads() && mylist_->apply()){
		return ListBlochParamsIterator<NullGrid, BPops, offset_T>::offset()+(GradientGridIter<GridEngine_t>::offset()*gammaGauss());
	}else{
		return ListBlochParamsIterator<NullGrid, BPops, offset_T>::offset();
	}
}


template<class GridEngine_t, int BPops, class offset_T>
inline double ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::offset(int i)
{
	if(mylist_->Grads() && mylist_->apply()){
		return ListBlochParamsIterator<NullGrid, BPops, offset_T>::offset(i)+(GradientGridIter<GridEngine_t>::offset(i)*gammaGauss());
	}else{
		return ListBlochParamsIterator<NullGrid, BPops, offset_T>::offset(i);
	}
}

template<class GridEngine_t, int BPops, class offset_T>
inline rmatrix ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::Jacobian()
{
	static rmatrix out(3,3,0);
	out(0,0)=(T2())?-1.0/T2():0.0;		out(0,1)=offset();		out(0,2)=0.0;
	out(1,0)=-offset();		out(1,1)=(T2())?-1.0/T2():0.0;		out(1,2)=0.0;
	out(2,0)=0.0;			out(2,1)=0.0;		out(2,2)=(T1())?1.0/T1():0.0;
	return out;
}

template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::Jacobian(rmatrix &out)
{
	out(0,0)-=(T2())?1.0/T2():0.0;		out(0,1)+=offset();
	out(1,0)-=offset();		out(1,1)-=(T2())?1.0/T2():0.0;
													out(2,2)-=(T1())?1.0/T1():0.0;
}

template<class GridEngine_t, int BPops, class offset_T>
inline rmatrix ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::Jacobian(int i)
{
	static rmatrix out(3,3,0);
	out(0,0)=(T2(i))?-1.0/T2(i):0.0;		out(0,1)=offset(i);		out(0,2)=0.0;
	out(1,0)=-offset(i);		out(1,1)=(T2(i))?-1.0/T2(i):0.0;		out(1,2)=0.0;
	out(2,0)=0.0;			out(2,1)=0.0;		out(2,2)=(T1(i))?1.0/T1(i):0.0;
	return out;
}

template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::Jacobian(rmatrix &out, int i)
{
	out(0,0)-=(T2(i))?1.0/T2(i):0.0;		out(0,1)+=offset(i);
	out(1,0)-=offset(i);		out(1,1)-=(T2(i))?1.0/T2(i):0.0;
													out(2,2)-=(T1(i))?1.0/T1(i):0.0;
}
*/

//The Advance iterator..simply passes it off to the List Blach <NULL> and the Gradient Grid itterators
//as they better have the same number of elements

//this is the "++myIter" (faster)
template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::operator++()
{
	ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator++();
	GradientGridIter<GridEngine_t>::operator++();
}

template<class GridEngine_t, int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<GradientGrid<GridEngine_t>, BPops, offset_T >::offset(int i)
{ return gammaGauss(i)*grads_->offset(i);	}

template<class GridEngine_t, int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<GradientGrid<GridEngine_t>, BPops, offset_T >::offset()
{ return gammaGauss()*grads_->offset();	}

//TIME DEPENDANT
template<class GridEngine_t, int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<GradientGrid<GridEngine_t>, BPops, offset_T >::offset(int i, double t)
{ return gammaGauss(i)*grads_->offset(i,t);	}

template<class GridEngine_t, int BPops, class offset_T>
inline offset_T ListBlochParamsIterator<GradientGrid<GridEngine_t>, BPops, offset_T >::offset(double t)
{ return gammaGauss()*grads_->offset(t);	}


//this is the "myIter++" (a bit slower, but much easier to type and think about)
template<class GridEngine_t, int BPops, class offset_T>
inline void ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::operator++(int){ 	operator++();	}


//an operator the function as the "have we hit the end yet?"
// an examples....>  while(myIter){...}
template<class GridEngine_t, int BPops, class offset_T>
inline ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::operator bool()
{
	return (ListBlochParamsIterator<NullGrid, BPops, offset_T>::operator bool()) || ( GradientGridIter<GridEngine_t>::operator bool() );
}

//reset the iterator after you have used it once (typically not used,
// but is it on occasion and thus saves a bit of initalization time)
template<class GridEngine_t, int BPops, class offset_T>
void ListBlochParamsIterator< GradientGrid<GridEngine_t>, BPops, offset_T >::reset()
{
	GradientGridIter<GridEngine_t>::reset();
	ListBlochParamsIterator<NullGrid, BPops, offset_T>::reset();
}

END_BL_NAMESPACE


#endif




