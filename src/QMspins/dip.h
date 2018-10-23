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
 dip.h -->  spin interactions dip
*/

#ifndef _dip_h_
#define _dip_h_ 1

#include "container/Vector/Vector.h"
#include "container/complex.h"
#include "QMspins/spinsys.h"
#include "container/matrix/matrix.h"
#include "utils/constants.h"
#include "utils/utils.h"
#include "QMspins/spin_ten.h"
#include "QMspins/space_ten.h"


BEGIN_BL_NAMESPACE


class Dip {
	private:
		double del_;
		//orientation with respect to another tensor
		double si_;
		double chi_;
		double psi_;
		//the spin number in a series if there is one
		int on1_;
		int on2_;

		//the A2m rotation bits for moving from PAS to crystal
		Vector<complex> crystalAs;

		rmatrix PAS_;

	public:
		Dip();
		Dip(const Dip &cp):
			del_(cp.del_), si_(cp.si_), chi_(cp.chi_), psi_(cp.psi_),
			on1_(cp.on1_), on2_(cp.on2_), crystalAs(cp.crystalAs){} ;


		Dip(double del, int on1, int on2, double si=0, double chi=0, double psi=0):
			del_(del), si_(si), chi_(chi), psi_(psi),on1_(on1), on2_(on2),crystalAs(5,0.0)
		{
			setCrystalAs();
		}

		~Dip() {};

		//returns a std::string like 'D(spin1, spin2)'
		std::string name();

		void reset();
		//the might extrnal spin_system

		hmatrix T20; //the spin matrix
		void setSpinMats(SpinSys &A);


		Dip operator=(const Dip &another);
		bool operator==(const Dip &another);
		bool operator!=(const Dip &another);


		friend std::ostream& operator<<(std::ostream &otr,const Dip &);

		inline bool is_zero()	{return del_==0.0;}
		inline bool is_hetero(SpinSys &in){	return in(on1_)!=in(on2_);	}

		inline int &spinON1()		{return on1_;}
		inline int &spinON2()		{return on2_;}
		inline double &couple()	{return del_;}
		inline int &on1()		{	return on1_;		}
		inline int &on2()		{	return on2_;		}
		inline double &del()		{	return del_;	}
		inline double &si()		{	return si_;		}
		inline double &alpha()	{	return si_;		}
		inline double &chi()		{	return chi_;	}
		inline double &beta()		{	return chi_;	}
		inline double &psi()		{	return psi_;	}
		inline double &gamma()	{	return psi_;	}

		inline int spinON1()	const	{return on1_;}
		inline int spinON2()	const	{return on2_;}
		inline double couple()	const{return del_;}
		inline int on1()		const{	return on1_;		}
		inline int on2()		const{	return on2_;		}
		inline double del()		const{	return del_;	}
		inline double si()		const{	return si_;		}
		inline double alpha()	const{	return si_;		}
		inline double chi()		const{	return chi_;	}
		inline double beta()		const{	return chi_;	}
		inline double psi()		const{	return psi_;	}
		inline double gamma()	const{	return psi_;	}

		inline void spinON1(int in)		{	on1_=in;	setCrystalAs();	}
		inline void spinON2(int in)		{	on2_=in;	setCrystalAs();	}
		inline void couple(double in)	{	del_=in;	setCrystalAs();	}
		inline void on1(int in)			{	on1_=in;	setCrystalAs();			}
		inline void on2(int in)			{	on2_=in;	setCrystalAs();			}
		inline void del(double in)		{	del_=in;	setCrystalAs();		}
		inline void si(double in)		{	si_=in;		setCrystalAs();		}
		inline void alpha(double in)	{	si_=in;	setCrystalAs();			}
		inline void chi(double in)		{	chi_=in;	setCrystalAs();		}
		inline void beta(double in)		{	chi_=in;	setCrystalAs();		}
		inline void psi(double in)		{	psi_=in;	setCrystalAs();		}
		inline void gamma(double in)	{	psi_=in;	setCrystalAs();	}

		void setCrystalAs();

		hmatrix H(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots(phi, theta, gam, alpha, beta);
			return get_H(A, rots);
		}

		hmatrix Hamiltonian(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots(phi, theta, gam, alpha, beta);
			return get_H(A,rots);
		}

		hmatrix H(SpinSys &A,Rotations &therotation)
		{	return get_H(A, therotation);	}

		hmatrix Hamiltonian(SpinSys &A,Rotations &therotation)
		{	return get_H(A, therotation);	}


		hmatrix get_H(SpinSys &A,double alpha, double beta, double theta, double phi, double gam=0.0);
		hmatrix get_H(SpinSys &, Rotations &therotation);
		hmatrix get_HINT(SpinSys &A,double wr,double alpha1, double alpha2,
		                double beta,double theta, double phi);


		void setParam(const std::string, double);
		double getParam(const std::string);
		void display();
		void write(std::ostream &otr);
};

Vector<Dip> read_dip(const char *filename);
Vector<Dip> read_dip(std::string filename);
Vector<Dip> read_dip(std::ifstream &filename);
Vector<Dip> read_dip(const Vector<std::string> &ss);


END_BL_NAMESPACE


#endif

