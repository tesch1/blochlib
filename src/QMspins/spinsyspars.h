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
/* spinsyspars.h *-*/

/*this class holds the list of isotopes types...it serves as the base class for
  the spin operators that generates all the spinoperators for each spin */

#ifndef _SpinSyspars_h_
#define _SpinSyspars_h_ 1

#include "bloch/Isotope.h"
#include "container/Vector/Vector.h"
#include <string>

BEGIN_BL_NAMESPACE


class spin_sysPars {
	private:
		int nspins;			// Number of spins in the system
		Vector<Isotope> Isotopes;	// The spin isotopes in the system
		static std::string DefIso;		// Default isotope type
        Isotope retIso;

		const void SpinErr() const
		{
			std::cerr<<std::endl<<"Error: spin_sysPars"<<std::endl;
			std::cerr<<" Desired spin does not exsist..."<<std::endl;
			std::cerr<<" Create Before acsessing.."<<std::endl;
		}
public:



		spin_sysPars();
		spin_sysPars(int spins);
		spin_sysPars(const spin_sysPars& cp):
			nspins(cp.nspins), Isotopes(cp.Isotopes) {};

		~spin_sysPars();

		spin_sysPars &Params() ;
		spin_sysPars Params() const ;

		spin_sysPars& operator=(const spin_sysPars& sys);

		bool operator==(const spin_sysPars& sys) const;
		bool operator!=(const spin_sysPars &sys) const;

		//see is spin with int 'spin' exsits...
		inline bool checkSpin(int spin, bool toexit=true) const
		{
			if((spin>=0)&&(spin<nspins)) return true;
			SpinErr();
			if(toexit) exit(-1);
			return false;
		}

//number os spins in the system
		inline int spins() const{ return nspins;	}
		inline int size() const { return nspins;	}
		void resize(int i);

//Here are the carry over from the Isotopr class

		int HS() const;	//total spin space dimension
		inline int HS(int spin) const
		{
			if(checkSpin(spin))	return Isotopes[spin].HS();
			return 0;
		}

		inline double weight(int spin) const
		{
			if(checkSpin(spin))	return Isotopes[spin].weight();
			return 0.;
		}

		inline double weight() const;

		inline double mass(int spin) const
		{
			if(checkSpin(spin))	return Isotopes[spin].mass();
			return 0.;
		}

		inline double mass() const;

		inline std::string element(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].element();
			return "";
		}

		inline std::string name(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].name();
			return "";
		}

		//total system momentum
		std::string momentum() const;
		inline std::string momentum(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].momentum();
			return "";
		}

		//sets the isotope for spin 'spin'
		void isotope(int spin, const std::string&);
		void isotope(int spin, const Isotope&);

		const Isotope& isotope(int spin) const
		{
			if(checkSpin(spin)) return Isotopes.data()[spin];
			return retIso;
		}

		inline std::string symbol(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].symbol();
			return "";
		}

		//total qn number
		double qn() const;

		inline double qn(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].qn();
			return 0.;
		}

		inline double gamma(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].gamma();
			return 0.;
		}

		inline double receptivity(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].receptivity();
			return 0.;
		}

		inline double relativeFrequency(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].relativeFrequency();
			return 0.;
		}

		inline double gammaGauss(int spin) const
		{
			if(checkSpin(spin)) return Isotopes[spin].gammaGauss();
			return 0.;
		}

		inline double gamma(const std::string& iso) const
		{
			return (Isotope(iso)).gamma();
		}

		inline double gammaGauss(const std::string &iso) const
		{
			return (Isotope(iso)).gammaGauss();
		}


		bool homonuclear() const;			//is the system homonuclear?
		bool heteronuclear() const;			//is the system heteronuclear?
		bool spinhalf() const;				//contains only spin 1/2's?


		int isotopes() ; 				//number of isotopes (different) in the system

		bool isotopes(const std::string& I) const;	//determins if an isotope exsists with label "I"

		inline static std::string IsoDefault(){	return DefIso;	}				//returns the default isotope label
		void IsoDefault(const std::string& DI){	DefIso=DI;		}		//sets the default isotope label

		std::ostream& print(std::ostream& out) const;
		friend std::ostream& operator<<(std::ostream& out, const spin_sysPars& sys);

		Vector<std::string> printstrings() const;

		Isotope &operator()(int i);
		Isotope operator()(int i) const;

		Isotope &operator[](int i);
		Isotope operator[](int i) const;

};

END_BL_NAMESPACE


#endif
