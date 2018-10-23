

/* blochParams.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-30-01
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
 	blochParamsBasic.h-->the Base class forthe basic storage blocks
 	maitains the things that are always true for a Bloch Parameter in
 	any circustance..
 */



#ifndef _bloch_params_Basic_h_
#define _bloch_params_Basic_h_ 1

#include "bloch/Isotope.h"
#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE

class BPoptions{
	public:
		enum Type{
			Particle=0x0001,
			Density=0x0002
		};

		enum Field{
			HighField=0x0010,
			ZeroField=0x0020
		};
};

//Some initial condition setting flags...
class InitialCondition{
	public:
		enum RnotR{
			Random=0x10000,
			NotRandom=0x20000
		};

		enum IC{
			AllUp=0x00001 | NotRandom,
			AllDown=0x00002 | NotRandom,
			HalfUpHalfDown=0x00004 | NotRandom,
			RandomDistribution=0x00008 | Random,
			RandomUpDown=0x00010 | Random,
			RandomUp=0x00020 | Random,
			RandomDown=0x00040 | Random
		};

		static int getFlag(std::string initc){
			if(initc=="RandomUpDown"){
				return RandomUpDown;
			}else if(initc=="Random" || initc=="RandomDistribution"){
				return RandomDistribution;
			}else if(initc=="RandomUp" ){
				return RandomUp;
			}else if(initc=="RandomDown" ){
				return RandomDown;
			}else if(initc=="AllDown"){
				return AllDown;
			}else{
				return AllUp;
			}
		}

		static std::string name(int IC){
			switch(IC){
				case AllDown: return "AllDown";
				case HalfUpHalfDown: return "HalfUpHalfDown";
				case RandomDistribution: return "RandomDistribution";
				case RandomUp: return "RandomUp";
				case RandomUpDown: return "RandomUpDown";
				case RandomDown: return "RandomDown";
				default: return "AllUp";
			}
		}


};


template<class Offset_T=double>
class BasicBlochParams;


//The BASIC BLoch Param....those things that do not change regardless
//of fields, or prticle/density parts

template<>
class BasicBlochParams<double > :
	public Isotope
{
	private:
		double omega_;		//Zeeman frequency
		double moles_;		//number of moles of sample
		double Mo_;			//initiall bulk magnitization amount of the sample

	public:
		double Bo_;			//Magnetic firld strength in TESLA
		double Ho_ ;		// = 1/permVac Bo
		double temperature_;	//temperature of the sample (Kelvin)

		typedef double Offset_T;

		BasicBlochParams():
			Isotope(),
			moles_(1.0),// T2_(0.0), T1_(0.0),
			Mo_(1.0),
			Bo_(4.7), Ho_(1./permVac *Bo_), temperature_(300)
		{

		}

		BasicBlochParams(std::string iso):
			Isotope(iso),
			moles_(1.0), //T2_(0.0), T1_(0.0),
			Mo_(1.0),
			Bo_(4.7), Ho_(1./permVac *Bo_), temperature_(300)
		{

		}

		BasicBlochParams(const BasicBlochParams &cp):
			Isotope(cp),
			moles_(cp.moles_),// T2_(cp.T2_), T1_(cp.T1_),
			Mo_(cp.Mo_),
			Bo_(cp.Bo_), Ho_(1./permVac *Bo_), temperature_(cp.temperature_)
		{

		}

		void setInitialCondition(int IC);

		inline double &Bo() { return Bo_;	}
		inline double Bx() { return 0.0;	}
		inline double By() { return 0.0;	}
		inline double &Bz() { return Bo_;	}

		inline double &Ho() { return Ho_;	}
		inline double Hx() { return 0.0;	}
		inline double Hy() { return 0.0;	}
		inline double &Hz() { return Ho_;	}

		inline double &omega()  {	return omega_;	}

		inline double &temperature() { return temperature_;	}
		inline double &Mo() { return Mo_;	}

		inline double &moles() { return moles_;	}

		inline double Bo()  const { return Bo_;	}
		inline double Bx() const { return 0.0;	}
		inline double By() const { return 0.0;	}
		inline double Bz() const { return Bo_;	}

		inline double Ho()  const { return Ho_;	}
		inline double Hx() const { return 0.0;	}
		inline double Hy() const { return 0.0;	}
		inline double Hz() const { return Ho_;	}

		inline double omega()   const {	return omega_;	}
		inline double temperature()  const { return temperature_;	}

		inline double Mo() const  { return Mo_;	}
		inline double moles() const  { return moles_;	}

		void operator=(const BasicBlochParams &rhs);
		BasicBlochParams &operator=(const std::string &rhs);

		bool operator==(const BasicBlochParams &rhs);
		bool operator!=(const BasicBlochParams &rhs);
};

template<>
class BasicBlochParams<coord<> > :
	public Isotope
{
	private:
		coord<> omega_;		//Zeeman frequency
		double moles_;		//number of moles of sample
	//	double T2_;			// T2 relaxation time
	//	double T1_;			// T1 relaxation time
		coord<> Mo_;			//initiall bulk magnitization amount of the sample

	public:
		coord<> Bo_;			//Magnetic firld strength in TESLA
		coord<> Ho_ ;		// = 1/permVac Bo
		double temperature_;	//temperature of the sample (Kelvin)

		typedef coord<> Offset_T;

		BasicBlochParams():
			Isotope(),
			moles_(1.0), //T2_(0.0), T1_(0.0),
			Mo_(0.0,0.0,1.0),
			Bo_(0.0,0.0,4.7), Ho_(1./permVac *Bo_), temperature_(300)
		{

		}

		BasicBlochParams(std::string iso):
			Isotope(iso),
			moles_(1.0), //T2_(0.0), T1_(0.0),
			Mo_(0.0,0.0,1.0),
			Bo_(0.0,0.0,4.7), Ho_(1./permVac *Bo_), temperature_(300)
		{

		}

		BasicBlochParams(const BasicBlochParams &cp):
			Isotope(cp),
			moles_(cp.moles_),// T2_(cp.T2_), T1_(cp.T1_),
			Mo_(cp.Mo_),
			Bo_(cp.Bo_), Ho_(1./permVac *Bo_), temperature_(cp.temperature_)
		{

		}

		inline coord<> &Bo() { return Bo_;	}
		inline double &Bx() { return Bo_.x();	}
		inline double &By() { return Bo_.y();	}
		inline double &Bz() { return Bo_.z();	}

		inline coord<> &Ho() { return Ho_;	}
		inline double &Hx() { return Ho_.x();	}
		inline double &Hy() { return Ho_.y();	}
		inline double &Hz() { return Ho_.z();	}

		inline coord<> &omega()  {	return omega_;	}

		inline double &temperature() { return temperature_;	}
		inline coord<> &Mo() { return Mo_;	}

		inline double &moles() { return moles_;	}
//		inline double &T2(){ return T2_;	}
//		inline double &T1(){ return T1_;	}

		inline coord<> Bo()  const { return Bo_;	}
		inline double Bx() const { return Bo_.x();	}
		inline double By() const { return Bo_.y();	}
		inline double Bz() const { return Bo_.z();	}

		inline coord<> Ho()  const { return Ho_;	}
		inline double Hx() const { return Ho_.x();	}
		inline double Hy() const { return Ho_.y();	}
		inline double Hz() const { return Ho_.z();	}

		inline coord<> omega()   const {	return omega_;	}
		inline double temperature()  const { return temperature_;	}

		inline coord<> Mo() const  { return Mo_;	}
		inline double moles() const  { return moles_;	}
//		inline double T2() const { return T2_;	}
//		inline double T1() const { return T1_;	}

//		void T2(double inT2);
//		void T1(double inT1);

		void operator=(const BasicBlochParams &rhs);
		BasicBlochParams &operator=(const std::string &rhs);

		bool operator==(const BasicBlochParams &rhs);
		bool operator!=(const BasicBlochParams &rhs);
};



END_BL_NAMESPACE


#endif

