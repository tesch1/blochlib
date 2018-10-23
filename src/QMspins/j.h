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
 j.h --> spin interactions scaler coupling
*/

#ifndef _j_h_
#define _j_h_ 1

#include "container/Vector/Vector.h"
#include "container/complex.h"
#include "QMspins/spinsys.h"
#include "container/matrix/matrix.h"
#include "utils/constants.h"
#include "utils/utils.h"
#include "QMspins/spin_ten.h"
#include "QMspins/space_ten.h"


BEGIN_BL_NAMESPACE


class J {
	private:

	//the coupling
		double iso_;

	//the spin number in a series if there is one
		int on1_;
		int on2_;
	//strong or weak coupling
		int strong_;

	public:
		J();
		J(const J &cp):
			iso_(cp.iso_), on1_(cp.on1_), on2_(cp.on2_), strong_(cp.strong_),
			T00(cp.T00)
		{};

		J(double iso, int on1, int on2, int strong=0):
			iso_(iso), on1_(on1), on2_(on2), strong_(strong)
		{};

		~J() {};

		//returns a std::string like 'J(spin1, spin2)'
		std::string name();

		void reset();
		hmatrix T00; //the spin matrix
		void setSpinMats(SpinSys &A);
		void setCrystalAs(){}

		inline bool is_zero() { return iso_==0 ;}

		J operator=(const J &another);
		bool operator!=(const J &another);
		bool operator==(const J &another);

		hmatrix H(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0)
		{	return get_H(A);	}

		hmatrix Hamiltonian(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0)
		{	return get_H(A);	}

		hmatrix H(SpinSys &A,Rotations &therotation)
		{	return get_H(A);	}

		hmatrix Hamiltonian(SpinSys &A,Rotations &therotation)
		{	return get_H(A);	}


		hmatrix get_H(SpinSys &);


		void setParam(const std::string, double);
		double getParam(const std::string);
		inline int spinON1() 	const{	return on1_;	}
		inline int spinON2() 	const{	return on2_;	}
		inline int on1()		const{	return on1_;	}
		inline int on2()		const{	return on2_;	}
		inline double couple()	const{	return iso_;	}
		inline double iso()		const{	return iso_;	}
		inline int isstrong()	const{	return strong_;	}
		inline int strong()		const{	return strong_;	}

		inline int &spinON1() 	{	return on1_;	}
		inline int &spinON2() 	{	return on2_;	}
		inline int &on1()		{	return on1_;	}
		inline int &on2()		{	return on2_;	}
		inline double &couple()	{	return iso_;	}
		inline double &iso()		{	return iso_;	}
		inline int &isstrong()	{	return strong_;	}
		inline int &strong()		{	return strong_;	}

		inline void spinON1(int in) 	{	on1_=in;	}
		inline void spinON2(int in) 	{	on2_=in;	}
		inline void on1(int in)		{	on1_=in;	}
		inline void on2(int in)		{	on2_=in;	}
		inline void couple(double in)	{	iso_=in;	}
		inline void iso(double in)		{	iso_=in;	}
		inline void strong(int in)		{	strong_=in;	}

		friend std::ostream& operator<<(std::ostream &otr,const J &);
		void write(std::ostream &otr);
		void display();

};

//this function is kept outside the class becuase it must be asscessed
//without knowledge of what 'Js' exsist
Vector<J> read_js(const char *filename);
Vector<J> read_js(std::string filename);
Vector<J> read_js(std::ifstream &filename);
Vector<J> read_js(const Vector<std::string> &ss);

END_BL_NAMESPACE


#endif

