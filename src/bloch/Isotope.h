 /* Isotope.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-4-01
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
	Isotope.h--> atomic isotopes data storage container....

	thanks to 'Gamma'...

I'd like to thank to the makers of Gamma

  S.A. Smith, T.O. Levante, B.H. Meier and R.R. Ernst
  Computer Simulations in Magnetic Resonance:
  An Object-Oriented Programming Approach
  J. Magn. Reson., 106a, 75-105, (1994S, Series A)

http://gamma.magnet.fsu.edu/

for the hard data

*/

#ifndef _Isotope_h_
#define _Isotope_h_ 1


#include "utils/constants.h"
#include "utils/utils.h"
#include <string>
#include <fstream>


BEGIN_BL_NAMESPACE


// Slew of the spin contant parameters
/*          An array of known isotope spin Hilbert spaces (2*I+1)           */

extern const int NIso;

extern const int IsoSpins[131];

/*          An array of known isotope atomic numbers (1=H, 6=C, ...)        */

extern const int IsoNumbers[131];

/*     An array of known isotope atomic masses  (1=1H, 2=2H, 13=13C, ...)   */

extern const int IsoMasses[131];

/*   An array of known isotope atomic weights (1H = 1.008665 g/mole, ...)   */

extern const double IsoWeights[131];

/*        An array of known isotope receptivities (based on 13C = 1)        */

extern const double IsoRecepts[131];

/*  An array of known isotope relative frequencies (from 1H @ 400.130 MHz)  */

extern const double IsoRelFreqs[131];

/*                       An array of known isotope names                     */

extern const std::string IsoNames[87];

/*                       Known Isotope Indices For Previous Names Array   */

extern const int IsoNamesIdx[131];


/* 					Gamma Factors in Hz(rads)/Telsa and Hz(rads)/Gauss  */
extern  const double IsoGammaPerTesla[131];
extern const double IsoGammaPerGauss[131];



extern const std::string IsoElements[85];
extern const int IsoIndexElements[131] ;

/*				Symbol names "1H", "2H", etc			*/
extern const std::string IsoSymbol[131];


class Isotope{
	private:
		static double  RelativeRF_1H ;		//realtive base frequency for proton
											//becuase it is static, it will be defined accross
											//all instanses of 'Isotope'
		int Index;	//the index that pts to the correct array elements in the data above

		//finds the correct index of the array's above give a std::string "1H", "13C", etc
		int getCorrectIndex(std::string sname);
		void setCorrectIndex(std::string sname);

	public:
		void SpinType(std::string sname);
		void SpinType(char *sname);

		Isotope();
		Isotope(const std::string &SpinName);
		Isotope(const Isotope &copy);

		Isotope &operator=(const Isotope &rhs);
		Isotope &operator=(const std::string &rhs);
		inline bool operator==(const Isotope &rhs)const{ return Index==rhs.Index;	}
		inline bool operator!=(const Isotope &rhs)const{ return Index!=rhs.Index;	}

		int HS() const;			//hilbert space size
		double spin() const;		//spin quantum number
		double qn() const;		//spin quantum number
		int number() const;		//atomic number
		double mass() const;		//amu mass
		double weight() const;	//g/mole weight
		double receptivity() const;	//relative sensitivity
		std::string name() const;		//"Hydrogen", "Carbon", etc
		std::string element() const;	//"H", "C", etc
		std::string symbol() const;	//"1H", "2H", etc
		double gamma() const;			//gamma factor
		double gammaGauss() const;
		double relativeFrequency() const;	//frequency relative to "1H" frequency

		void Spin(std::string in);

		std::string momentum() const;
		friend std::ostream &operator<<(std::ostream &oo, const Isotope &out);
		void print(std::ostream &oo);
		void print(){ print(std::cout);	}
};

END_BL_NAMESPACE



#endif


