/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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
 csa.h --> spin interaction CSA
*/
#ifndef _csa_h_
#define _csa_h_ 1

#include "container/Vector/Vector.h"
#include "container/complex.h"
#include "QMspins/spinsys.h"
#include "container/matrix/matrix.h"
#include "utils/constants.h"
#include "utils/utils.h"
#include "QMspins/spin_ten.h"
#include "QMspins/space_ten.h"

BEGIN_BL_NAMESPACE


class Csa {
	private:
		//these are the general measured values
		double iso_;
		double del_;
		double eta_;
		double sig1_;
		double sig2_;
		double sig3_;
		//orientation with respect to another tensor
		double si_;
		double chi_;
		double psi_;
		//the spin number in a series if there is one
		int chain_;

		//the A2m rotation bits for moving from PAS to crystal
		Vector<complex> crystalAs;
		rmatrix PAS_;
	public:
		Csa();
		Csa(const Csa &cp):
			iso_(cp.iso_),del_(cp.del_), eta_(cp.eta_), sig1_(cp.sig1_), sig2_(cp.sig2_), sig3_(cp.sig3_),
			si_(cp.si_), chi_(cp.chi_), psi_(cp.psi_), chain_(cp.chain_), crystalAs(cp.crystalAs),
			T20(cp.T20), T00(cp.T00)
		{} ;

		Csa(double iso, double del, double eta, int chain,
			double si=0, double chi=0, double psi=0):
			iso_(iso), del_(del), eta_(eta),
			si_(si), chi_(chi), psi_(psi),
			crystalAs(5,0.0)
		{
			setCrystalAs();
		}

		~Csa() {};

		//returns a std::string like 'C(spin1)'
		std::string name();

		void reset();
		rdmatrix T20;  //the spin matrix for the CSA
		rdmatrix T00;	//oth rank spin mat
		void setSpinMats(SpinSys &A);

		inline bool is_zero() 	{ 	return del_==0 && iso_==0;}

		inline int &spinON()		{ 	return chain_; }
		inline double &iso_cop() 	{ 	return iso_;	}
		inline double &del_cop() 	{ 	return del_;	}
		inline double &eta_cop() 	{ 	return eta_;	}
		inline double &iso()		{	return iso_;	}
		inline double &del()		{	return del_;	}
		inline double &eta()		{	return eta_;	}
		inline int &on()			{	return chain_;	}
		inline double &si()		{	return si_;		}
		inline double &alpha()	{	return si_;		}
		inline double &chi()		{	return chi_;	}
		inline double &beta()		{	return chi_;	}
		inline double &psi()		{	return psi_;	}
		inline double &gamma()	{	return psi_;	}

		inline int spinON()		const{ 	return chain_; }
		inline double iso_cop() 	const{ 	return iso_;	}
		inline double del_cop() 	const{ 	return del_;	}
		inline double eta_cop() 	const{ 	return eta_;	}
		inline double iso()		const{	return iso_;	}
		inline double del()		const{	return del_;	}
		inline double eta()		const{	return eta_;	}
		inline int on()			const{	return chain_;	}
		inline double si()		const{	return si_;		}
		inline double alpha()	const{	return si_;		}
		inline double chi()		const{	return chi_;	}
		inline double beta()	const	{	return chi_;	}
		inline double psi()		const{	return psi_;	}
		inline double gamma()	const{	return psi_;	}

		inline void iso(double in)		{	 iso_=in; setCrystalAs();	}
		inline void del(double in)		{	 del_=in; setCrystalAs();	}
		inline void eta(double in)		{	 eta_=in; setCrystalAs();	}
		inline void on(int in)			{	 chain_=in;	}
		inline void si(double in)		{	 si_=in; setCrystalAs();		}
		inline void alpha(double in)	{	 si_=in; setCrystalAs();		}
		inline void chi(double in)		{	 chi_=in; setCrystalAs();	}
		inline void beta(double in)		{	 chi_=in; setCrystalAs();	}
		inline void psi(double in)		{	 psi_=in; setCrystalAs();	}
		inline void gamma(double in)	{	 psi_=in; setCrystalAs();	}


		void setCrystalAs();

		dmatrix H(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots( phi, theta,gam, alpha, beta);
			return get_H(A, rots);
		}

		dmatrix Hamiltonian(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots( phi, theta, gam, alpha, beta);
			return get_H(A, alpha, beta, theta, phi,gam);
		}

		dmatrix H(SpinSys &A,Rotations &therotation)
		{	 return get_H(A, therotation);	}

		dmatrix Hamiltonian(SpinSys &A,Rotations &therotation)
		{	return get_H(A, therotation);	}


		dmatrix get_H(SpinSys &,double alpha, double beta, double theta, double phi, double gam=0.0);
		dmatrix get_H(SpinSys &, Rotations &therotation);
		dmatrix get_HINT(SpinSys &, double wr, double alpha1, double alpha2, double beta, double theta, double phi);


		void setParam(const std::string, double);
		void setSigma();

		double getParam(const std::string param);

		Csa operator=(const Csa &another);
		friend std::ostream& operator<<(std::ostream &otr, const Csa &csa);
		bool operator==(const Csa &rhs);
		bool operator!=(const Csa &rhs);

		void display();
		void write(std::ostream &otr); //writes a line corresponding to the input
};

Vector<Csa> read_csa(const char *filename);
Vector<Csa> read_csa(std::string filename);
Vector<Csa> read_csa(std::ifstream & infile);
Vector<Csa> read_csa(const Vector<std::string> & infile);

END_BL_NAMESPACE



#endif

