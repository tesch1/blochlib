

#ifndef _spin_tensors_index__h_
#define _spin_tensors_index__h_ 1

//#include "fullbinarytree.h"
#include "blochlib.h"
#include <algorithm>

using namespace BlochLib;

/*

This is a set of classes that generate and maintain massive lists of
various combinations of spin tensors to certain 'orders'
i.e.
1st order--> Ie, Ix, Iy, Iz
2nd Order--> Ix^2, Iy^2, Iz^2, Ix*Iy.... (al powers of 2)
3rd order-->all powers of three
etc

there are 2 brands...

1) the Cartiesian tensors (Ix, Iy, Iz representation)
2) the Spherical Tensors (T00, T10, T11, T1m1..)





*/

//options for a single Spin Index
class SpinIdxOps
{
	public:
		enum{
				Commutators=0x0000001,
				Multiply=	0x0000002,
				Add=		0x0000004
		}MulType;

		enum{	T=			0x0000010,
				spherical=	0x0000010,
				I=			0x0000020,
				carteisian=	0x0000020,
				Hamiltonians=0x0000040

			}matrixType;

};

class Tensors{
	public:
		enum{	x=			1,
				y=			2,
				z=			3
			} axis;


		enum{	T1=			0x1000,
				T2=			0x2000
			} Ts;

		enum{	T1m1=		4 | T1,
				T10=		5 | T1,
				T11=		6 | T1
			} T11s;

		enum{	T2m2=		7 | T2,
				T2m1=		8 | T2,
				T20=		9 | T2,
				T21=		10 | T2,
				T22=		11 | T2
			} T22s;
		enum{	Hamiltonians	= SpinIdxOps::Hamiltonians	};

};

//the is the option class
class TensorGenOps{
	public:

	//what 'brand' of spin tensors to permute
		enum{	T=SpinIdxOps::T,	//spherical tensors (both first and second if possible)
				T2=Tensors::T2,	//spherical tensors second rank only
				T1=Tensors::T1,	//spherical tensors first rank only
				I=SpinIdxOps::I,	//simple cartesian spin ops
				carteisian=SpinIdxOps::I,
				spherical=SpinIdxOps::T,
				Hamiltonians=Tensors::Hamiltonians
			}Brand;

		enum{	//PermuteSpins=	0x00000001, 	//permute all the spins in the System
				TotalSpace=		0x0001000,		//or simply choose 'total' space
				ChooseSpins=	0x0002000,		//choose a vector of spins
				SingleSpin=		0x0004000		//a single spins
				//TwoSpins=		0x00000010			//choose a pair to permute
			}SpinPermutions;

		enum{	Commutators=	SpinIdxOps::Commutators,	//spinmatrices are 'commutated' [T, [T, T]...]]
				Multiply=		SpinIdxOps::Multiply,		//spinmats are simply multiplies
				Add=			SpinIdxOps::Add			//spinmats are added together
			}MulType;
};



//data holder for the Cartesian types
class SpinIdxCartData{
	private:

	//depending on the order of the tensor (1st,...Nth) contains the
	//order of which spin in a SpinSys is each tensor
		Vector<int> spins_;

	//depending on the order of the tensor (1st,...Nth) contains the
	//order of axis multiplication
		Vector<int> axis_;

	//the multiplication type
		int MulOps;

	public:
		static const int Type=SpinIdxOps::I; //this is a cartesian item


		inline SpinIdxCartData()
		{}

		inline SpinIdxCartData(const SpinIdxCartData &cp):
			spins_(cp.spins_),
			axis_(cp.axis_),
			MulOps(cp.MulOps)
		{}

		inline SpinIdxCartData(int spins,const Vector<int> &axis, int ops=SpinIdxOps::Multiply ):
			spins_(1, spins),
			axis_(axis),
			MulOps(ops)
		{}

		inline SpinIdxCartData(const Vector<int> &spins,const Vector<int> &axis, int ops=SpinIdxOps::Multiply ):
			spins_(spins), axis_(axis),
			MulOps(ops)
		{}


		void setSpin(int i);
		void setSpin(const Vector<int> &in);

		void setOptions(int ops);

		inline Vector<int> spins() const {	return spins_;	}

		bool operator==(const SpinIdxCartData &rhs);
		bool operator!=(const SpinIdxCartData &rhs);

		//Camparisons for the Sorting operations (sorts by the 'length' of axis;
		friend bool operator<(const SpinIdxCartData &lhs,const SpinIdxCartData &rhs);
		friend bool operator>(const SpinIdxCartData &lhs,const SpinIdxCartData &rhs);

		void operator=(const SpinIdxCartData &rhs);

		inline int  size() const {	return axis_.size();	}


		void setAxes(const Vector<int> &in);
		inline Vector<int> axes() const {	return axis_;	}

		std::string name() const ;

		matrix operator()(SpinSys &in);
		matrix tensor(SpinSys &in);

	//i/o pieces
		int binarySize();

	//format is...
	// "SCD<binarysize><numele><numsp1><spinarray1><numaxis><axis1>...<axis2>"
		void read(std::fstream &in); //binary file reader
		void write(std::fstream &out); //binary file writer
};



std::ostream &operator<<(std::ostream &oo, const SpinIdxCartData &out);

//binary out operator
std::fstream &operator<<(std::fstream &oo, SpinIdxCartData &out);

//binary read operator
std::fstream &operator>>(std::fstream &oo, SpinIdxCartData &out);


/***************************************************************/
/******** SPHERICAL DATA ELEMENT **************/
/***************************************************************/

//data holder for the Spherical types
class SpinIdxSphData{
	friend class TensorGen;
	private:

	//depending on the order of the tensor (1st,...Nth) contains the
	//order of which spin in a SpinSys is each tensor
	//
	//if spin[i]=-1. the 'full' system matrix is implied
	//
	// NOTE: both spins1_ and spins2_ must be the same length
	//
		Vector<int> spins1_; //used if ranks ones
		Vector<int> spins2_; //needed for the second rank tensors

	//depending on the order of the tensor (1st,...Nth) contains the
	//order of axis multiplication
	// Sphericals are treated like using hte SpinIdxOps class above
	//
		Vector<int> axis_;

		int MulOps; //this is the multiplication type for the tensors

	public:
		static const int Type=SpinIdxOps::T; //this is a cartesian item


		inline SpinIdxSphData()
		{}

		inline SpinIdxSphData(const SpinIdxSphData &cp):
			spins1_(cp.spins1_),
			spins2_(cp.spins2_),
			axis_(cp.axis_),
			MulOps(cp.MulOps)
		{}

		inline SpinIdxSphData(int spins,const Vector<int> &axis, int ops=SpinIdxOps::Multiply ):
			spins1_(1, spins),
			spins2_(1, spins),
			axis_(axis),
			MulOps(ops)
		{}

		inline SpinIdxSphData(int spins, int spins2, const Vector<int> &axis, int ops=SpinIdxOps::Multiply ):
			spins1_(1, spins),
			spins2_(1, spins2),
			axis_(axis),
			MulOps(ops)
		{}

		inline SpinIdxSphData(const Vector<int> &spins, const Vector<int> &spin2, const Vector<int> &axis, int ops=SpinIdxOps::Multiply ):
			spins1_(spins),spins2_(spin2), axis_(axis),
			MulOps(ops)
		{
			RunTimeAssert(spins1_.size()==spins2_.size());
		}


		void setSpin(int i);
		void setSpin(int i, int j);
		void setSpin(const Vector<int> &in, const Vector<int> &sp2);

		void setOptions(int ops);

		inline Vector<int> spins1() const {	return spins1_;	}
		inline Vector<int> spins2() const {	return spins2_;	}

		bool operator==(const SpinIdxSphData &rhs);
		bool operator!=(const SpinIdxSphData &rhs);

		//Camparisons for the Sorting operations (sorts by the 'length' of axis;
		friend bool operator>(const SpinIdxSphData &lhs,const SpinIdxSphData &rhs);
		friend bool operator<(const SpinIdxSphData &lhs,const SpinIdxSphData &rhs);


		void operator=(const SpinIdxSphData &rhs);

		inline int  size() const {	return axis_.size();	}

		void setAxes(const Vector<int> &in);
		inline Vector<int> axes() const {	return axis_;	}

		std::string name() const ;

		matrix operator()(SpinSys &in);
		matrix tensor(SpinSys &in);

	//i/o
		int binarySize();

	//format is...
	// "SSD<binarysize><numele><numsp1><spinarray1><numsp2><spinarr2><numaxis><axis1>...<axis2>"
		void write(std::fstream &oo);
		void read(std::fstream &oo);

};

std::ostream &operator<<(std::ostream &oo, const SpinIdxSphData &out);

//binary out operator
std::fstream &operator<<(std::fstream &oo, SpinIdxSphData &out);

//binary read operator
std::fstream &operator>>(std::fstream &oo, SpinIdxSphData &out);

/***************************************************************/
/******** HamilTONian Solid Sys DATA ELEMENT **************/
/***************************************************************/

//this is a container for a single hamiltonian type
// where we need to generate something like
// [Hdip, Hcsa]
// i.e. specific hamiltonian rather then the
// tensorial componenets
class HamilTypes{
	public:
		enum{
			Csa=0x0001,
			Dip=0x0002,
			J=	0x0004,
			Q=	0x0008
		}Htypes;
};


class SolidSysData{
	private:

	// a single element has 2 nesesary lables.. 1) what type of hamiltonian
	//i.e. Dipole or CSA, etc
	//and the index of the sepcifi interaction (the label in the list of spins)
	// of a SolidSys
		Vector<int> Htypes_; //hamiltonian types storage as above
		Vector<int> indexs_; //indexs of the specifici interaction


		int MulOps; //this is the multiplication type for the tensors

		//pointer to the SoldSys

	public:
		static const int Type=TensorGenOps::Hamiltonians; //this is a cartesian item


		inline SolidSysData()
		{}

		inline SolidSysData(const SolidSysData &cp):
			Htypes_(cp.Htypes_),
			indexs_(cp.indexs_),
			MulOps(cp.MulOps)
		{}

		inline SolidSysData(int Hyt, int ind, int ops=SpinIdxOps::Multiply ):
			Htypes_(1, Hyt),
			indexs_(1, ind),
			MulOps(ops)
		{}

		inline SolidSysData(const Vector<int> &hty, const Vector<int> &inds, int ops=SpinIdxOps::Multiply ):
			Htypes_(hty),indexs_(inds),
			MulOps(ops)
		{
			RunTimeInsist(Htypes_.size()==indexs_.size(), "input arrays must be equal sizes");
		}


		void setHtypes(int i);
		void setHtypes(const Vector<int> &sp2);

		void setIndex(const Vector<int> &in);

		void setOptions(int ops);

		inline Vector<int> htypes() const {	return Htypes_;	}
		inline Vector<int> indexs() const {	return indexs_;	}

		bool operator==(const SolidSysData &rhs);
		bool operator!=(const SolidSysData &rhs);

		//Camparisons for the Sorting operations (sorts by the 'length' of axis;
		friend bool operator>(const SolidSysData &lhs,const SolidSysData &rhs);
		friend bool operator<(const SolidSysData &lhs,const SolidSysData &rhs);


		void operator=(const SolidSysData &rhs);

		inline int  size() const {	return Htypes_.size();	}

		std::string name(SolidSys &sys) const ;

	//this gets the hamiltonians...it assumes that the
	// angles for the powder or spinning have already been set
		matrix operator()(SolidSys &in);
		matrix tensor(SolidSys &in);

	//i/o
		int binarySize();

	//format is...
	// "SOS<binarysize><numele><numhtype><htypes><numindx><index>"
		void write(std::fstream &oo);
		void read(std::fstream &oo);

};

std::ostream &operator<<(std::ostream &oo, const SolidSysData &out);

//binary out operator
std::fstream &operator<<(std::fstream &oo, SolidSysData &out);

//binary read operator
std::fstream &operator>>(std::fstream &oo, SolidSysData &out);


/***************************************************************/
/******** Tensor Auto-Generator classes (the users main interface) **************/
/***************************************************************/


//The ,master generator class...
// the Default options are
// "Cartesian" "SingleSpin" "Multiply"

class TensorGen{
	private:

		Vector<SpinIdxSphData> sphere_;	//place to store spherical elements
		Vector<SpinIdxCartData> cart_; 	//store the cartiesion elements
		Vector<SolidSysData> solidsys_; 	//store the Solid Sys data elemetns

		int options_;	//stoarge for the input options

		Vector<int>  spinperms1_; //place to store the spins for 'ChoosSpins, Single Spin, and TwoSpins' options
		Vector<int>  spinperms2_; //place to store the spins for 'ChoosSpins, Single Spin, and TwoSpins' options

	//the order of the highest multiple
	// if order_=0..NO generation
	//    order_=1..thing ~Ix^1 or T20^1
	//    order_=2..things ~Ix^2 or T20^2
	//  ...
	// i.e. is determins the large multiplication of a series of tensors
		int order_;

	//these are the Base lists of generator arrays

		static const int cartarr_[3];
		static const int totalspharr_[8];
		static const int doublespharr_[5];
		static const int singlespharr_[3];

		int size_;

	//generates the initial permutation lists
	//returns the length of the list
		void genPermutationList(int length, int *, int *, int *);
		int permutationLength();

	//same functions except for the Hamiltonians' type
		int permutationLength(SolidSys &sys);
		void genPermutationList(SolidSys &, int length, int *, int *);


	//this generates the 'next k-subset' of a list of length 'n'
		void nextKSubset(int n, int k, int *subset, bool &more);


	//gets the 'multiplication option int' from the 'ops'
		int MulOption();

	//the private function generators
		void generateVectors();

	//the private function generators for Hamiltonian types
		void generateVectors(SolidSys &sys);

	public:
		TensorGen(): order_(0) {}

		inline TensorGen(int  order, int ops=TensorGenOps::I | TensorGenOps::SingleSpin | TensorGenOps::Multiply):
			options_(ops),
			spinperms1_(1, 0), //choose the '0th' spin as the default
			spinperms2_(1, 0), //choose the '0th' spin as the default
			order_(order)
		{
		//	generate();	//sets up the sphere_ and cart_ lists
		}

		inline TensorGen(int  order,  Vector<int> choos1, int ops=TensorGenOps::I | TensorGenOps::ChooseSpins  | TensorGenOps::Multiply):
			options_(ops),
			spinperms1_(choos1),
			spinperms2_(choos1),
			order_(order)
		{
		//	generate();	//sets up the sphere_ and cart_ lists
		}

		inline TensorGen(int  order,  Vector<int> choos1, Vector<int> choos2, int ops=TensorGenOps::I | TensorGenOps::ChooseSpins  | TensorGenOps::Multiply):
			options_(ops),
			spinperms1_(choos1),
			spinperms2_(choos2),
			order_(order)
		{
		//	generate();	//sets up the sphere_ and cart_ lists
		}


	//master generating function
		void generate();

	//master generating for a Solid Sys Data
		void generate(SolidSys &sys);

	//set the order
		void setOrder(int in){	order_=in;	}

	//set the parameters
		inline void setSpin(int i){	setSpins(i);	}
		void setSpins(int i);
		void setSpins(int i, int j);
		void setSpins(Vector<int> i);
		void setSpins(Vector<int> i, Vector<int> j);

	//set the options
		inline void setOptions(int in){	options_=in;	}

	//get the string name of the operator
		std::string name(int i);

	//get the string name of the operator
		std::string name(SolidSys &sys, int i);

	//get the operator from the list
		matrix tensor(SpinSys &sys, int i);
		matrix operator()(SpinSys &sys, int i);

	//get the operator from the list
		matrix tensor(SolidSys &sys, int i);
		matrix operator()(SolidSys &sys, int i);

	//size
		inline int size()const {	return size_;	}

	//size of the subelemenets
		int size(int i)const;

	//the approximate size guess
		inline int approxSize(int in){	return sum_factorial(in);	}

	//sorts the list
		void sort();

	//i/o
		int binarySize();

	//binary i/o...format
	//BEGIN TensorGen
	//<cart or sph><spin ops><numele><the indivdual data elements>
	//END TensorGen
		void write(std::fstream &oo);
		void read(std::fstream &oo);
};


//binary out operator
std::fstream &operator<<(std::fstream &oo, TensorGen &out);

//binary read operator
std::fstream &operator>>(std::fstream &oo, TensorGen &out);












#endif



