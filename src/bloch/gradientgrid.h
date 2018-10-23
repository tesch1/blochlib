

#ifndef _gradgrid_h_
#define _gradgrid_h_



#include "container/grids/grids.h"
#include "container/grids/coords.h"
#include "container/Vector/Vector.h"
#include <string>


/* This is a small class that holds the gradient strength
   the spin (i.e. "1H", "13C", etc) to which the gradient is applied
   in a list (i.e. there can be more then one spin)

    the class "GradientGrid" uses this as its little data stricture..

    SGradData--> holds a single gradient set and a single name
    GradData-->  holds the list of SGradData... (BUT IS NOT NEEDED FOR GRADIENTS!...a mild mistake, but there
    			is no sence is getting rid of good list like code....)
    GradientGrid --> uses the SGradData (i.e. same gradient everywhere)
*/

BEGIN_BL_NAMESPACE


class SGradData {
	private:
		coord<> G_;
		std::string name_;

	public:
		SGradData():
			G_(0), name_("1H") {}

		SGradData(coord<> &inG):
			G_(inG), name_("1H") {}

		SGradData(coord<> &inG, std::string inname):
			G_(inG), name_(inname) {}

		SGradData(double Gx, double Gy, double Gz):
			G_(Gx, Gy, Gz), name_("1H") {}

		SGradData(double Gx, double Gy, double Gz, std::string inname):
			G_(Gx, Gy, Gz), name_(inname) {}

		SGradData(const SGradData &cp):
			G_(cp.G_), name_(cp.name_){}

		SGradData &operator=(const SGradData &rhs);

		inline double Gx() const {	return G_.x();	}
		inline double Gy() const {	return G_.y();	}
		inline double Gz() const {	return G_.z();	}

		inline void Gx(double inGx)  {	return G_.x(inGx);	}
		inline void Gy(double inGy)  {	return G_.y(inGy);	}
		inline void Gz(double inGz)  {	return G_.z(inGz);	}

		inline std::string symbol() const			{ 	return name_;	}
		inline void symbol(std::string in)	{	name_=in;	}

		inline coord<> &G() {	return G_;	}
		inline void G(coord<> &inG_){	G_=inG_;	}
		inline void G(double Gx, double Gy, double Gz)	{	G_=coord<>(Gx, Gy, Gz);	}

		void print(std::ostream &oo);
		void print();
};


std::ostream &operator<<(std::ostream &oo, SGradData &out);


class GradData {
	private:
		Vector<SGradData > data_;


	public:
		GradData():
			data_() {}

		GradData(GradData &cp):
			data_(cp.data_){}

		GradData(coord<> &inG):
			data_(1, SGradData(inG)) {}

		GradData(coord<> &inG, std::string inname):
			data_(1, SGradData(inG, inname)) {}

		GradData(double Gx, double Gy, double Gz):
			data_(1, SGradData(Gx, Gy, Gz)) {}

		GradData(double Gx, double Gy, double Gz, std::string name):
			data_(1, SGradData(Gx, Gy, Gz, name)){ }

		inline int size()const {	return data_.size();	}

		inline SGradData &operator()(int i)
		{
			return data_(i);
		}

		inline SGradData &operator[](int i)
		{
			return data_(i);
		}

		GradData &operator=(const GradData &rhs);

		inline coord<> &G(std::string nm)
		{

			int i=Exists(nm);
			if(i!=-1) return data_(i).G();
			return ZeroType<coord<> >::zero();
		}

		inline void G(coord<> &in, std::string nm)
		{
			int i=Exists(nm);
			if(i!=-1) data_(i).G(in);
		}

		inline coord<> dBdr(std::string nm) 			{ return G(nm);		}
		inline void dBdr(coord<> &in, std::string nm) 	{ G(in, nm);			}



		inline double Gx(std::string nm)
		{
			int i=Exists(nm);
			if(i!=-1) return data_(i).Gx();
			return 0.;
		}

		inline void Gx(double inGx, std::string nm)
		{
			int i=Exists(nm);
			if(i!=-1)  data_(i).Gx(inGx);
		}

		inline double Gy(std::string nm)
		{
			int i=Exists(nm);
			if(i!=-1) return data_(i).Gy();
			return 0.;
		}

		inline void Gy(double inGy, std::string nm)
		{
			int i=Exists(nm);
			if(i!=-1)  data_(i).Gy(inGy);
		}

		inline double Gz(std::string nm)
		{
			int i=Exists(nm);
			if(i!=-1) return data_(i).Gz();
			return 0.;
		}

		inline void Gz(double inGz, std::string nm)
		{
			int i=Exists(nm);
			if(i!=-1)  data_(i).Gz(inGz);
		}

		inline int Exists(std::string nm)
		{
			for(int i=0;i<size();i++)
			{
				if(nm==data_(i).symbol()) return i;
			}
			return -1;
		}

		GradData &operator+(SGradData &rhs);
		GradData &operator+=(SGradData &rhs);

		void add(coord<> &inG, std::string inname);
		void add(double Gx, double Gy, double Gz, std::string inname);
		void add(SGradData &in);

		void print(std::ostream &oo);
		void print();

};


std::ostream &operator<<(std::ostream &oo, GradData &out);


/* Uses the Grid class as its base, simply adds on a gradient
	for each point in the grid in the grid...the grad simply acts
	like an offset term for each point in the grid

	It is assumed that the Units of the Gradient are in the standard G/cm (Gauss/cm)
	And the Grid sizes are in Centemeters
		(i.e. 5mm is ~ the width of a liquid NMR tube)

	*/

enum Graddirection{Gradz,Gradx,Grady};

//pre dec for typedef
template<class GridEngine_t>
class GradientGridIter;

template<class GridEngine_t>
class GradientGrid : public GridEngine_t, public SGradData {

	friend class GradientGridIter<GridEngine_t>;

	private:
		static const std::string Gunits_;
		static const std::string gridunits_;
		static double dott_;

	public:

		typedef typename GridEngine_t::iterator iterator;
		typedef typename GridEngine_t::Griditerator Griditerator;
		typedef double Offset_T;

		GradientGrid(double min, double max, int dims):
			GridEngine_t(min, max, dims), SGradData()
		{}

		GradientGrid(const GradientGrid &dup):
			GridEngine_t(dup),SGradData(dup)
		{}

		GradientGrid(GridEngine_t &dup):
			GridEngine_t(dup), SGradData()
		{}

		GradientGrid(double minr, double maxr, double minth, double maxth, int dims):
			GridEngine_t(minr, maxr, minth, maxth, dims), SGradData()
		{}

		GradientGrid(double minr, double maxr, double minth, double maxth,double minph, double maxph, int dims):
			GridEngine_t(minr, maxr, minth, maxth, minph, maxph, dims), SGradData()
		{}

		GradientGrid(double maxr, double minth, double maxth, int dims):
			GridEngine_t(maxr, minth, maxth, dims), SGradData()
		{}

		GradientGrid(double maxr, double minth, double maxth,double minph, double maxph, int dims):
			GridEngine_t(maxr, minth, maxth, minph, maxph, dims), SGradData()
		{}

		~GradientGrid(){}

//x(), y(), z() inherited from Engin_t
//G inherited from GradData
		inline double Xoffset(int i)		{ return x(i)*Gx();	}
		inline double Yoffset(int i)		{ return y(i)*Gy();	}
		inline double Zoffset(int i)		{ return z(i)*Gz();	}

		inline double offset(int i)			{ return dot(Point(i),G());}

//x(), y(), z() inherited from Engin_t
//G inherited from GradData
		inline double Xoffset(int i, double t)		{ return x(i,t)*Gx();	}
		inline double Yoffset(int i, double t)		{ return y(i,t)*Gy();	}
		inline double Zoffset(int i, double t)		{ return z(i,t)*Gz();	}

		inline double offset(int i, double t)			{ return dot(Point(i,t),G());}

		void operator=(const GradientGrid &rhs)
		{
			if(&rhs==this) return;
			GridEngine_t::operator=(rhs);
			SGradData::operator=(rhs);
		}

		void SetGrid(GridEngine_t &in)
		{
			GridEngine_t::operator=(in);
		}
		void setGrid(GridEngine_t &in)
		{	SetGrid(in);	}

		GridEngine_t *Grid(){	return this;	}
		const GridEngine_t *Grid() const {	return this;	}

		inline int size()  {	return GridEngine_t::size();	}


		void printData(std::ostream &oo, Graddirection dir=Gradz)
		{
			iterator myIt((*this));
			while(myIt){
				if(dir==Gradz){
					oo<<myIt.x()<<" "<<myIt.y()<<" "<<myIt.Zoffset()<<std::endl;
				}else if(dir==Grady){
					oo<<myIt.x()<<" "<<myIt.Yoffset()<<" "<<myIt.z()<<std::endl;
				}else if(dir==Gradx){
					oo<<myIt.Xoffset()<<" "<<myIt.y()<<" "<<myIt.z()<<std::endl;
				}
				++(myIt);
			}
		}


};


template<class Engine_t>
const std::string GradientGrid<Engine_t>::Gunits_="G/cm";

template<class Engine_t>
const std::string GradientGrid<Engine_t>::gridunits_="cm";

template<class Engine_t>
double GradientGrid<Engine_t>::dott_=0;


/**** Gradient grid itterator ****/
//the main iteration is taken care of by the GridEngine_t::iterator base class


template<class GridEngine_t>
class GradientGridIter: public GridEngine_t::iterator
{

	protected:

		GradientGrid<GridEngine_t> *TheList_;

	public:
		GradientGridIter():
			GridEngine_t::iterator()
		{}

		GradientGridIter(GradientGrid<GridEngine_t> &in):
			GridEngine_t::iterator(in), TheList_(&in)
		{}

		GradientGridIter(GradientGrid<GridEngine_t> *in):
			GridEngine_t::iterator(in), TheList_(in)
		{}


		~GradientGridIter()
		{
			TheList_=NULL;
		}

		//x(), y(), z() inherited from Engin_t
		//G inherited from GradData
		inline double Xoffset()		{ return x()*TheList_->Gx();	}
		inline double Yoffset()		{ return y()*TheList_->Gy();	}
		inline double Zoffset()		{ return z()*TheList_->Gz();	}

		inline double offset()		{ return dot(Point(),TheList_->G()); }

	//x(), y(), z() inherited from Engin_t
	//G inherited from GradData
		inline double Xoffset(double t)		{ return x( t)*TheList_->Gx();	}
		inline double Yoffset(double t)		{ return y( t)*TheList_->Gy();	}
		inline double Zoffset(double t)		{ return z( t)*TheList_->Gz();	}

		inline double offset(double t)		{ return dot(Point( t),TheList_->G()); }
};


END_BL_NAMESPACE




#endif

