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
	qua.cc--> methods spin interaction quadrupole
*/
#ifndef _qua_h_
#define _qua_h_

#include "container/Vector/Vector.h"
#include "container/complex.h"
#include "QMspins/spinsys.h"
#include "container/matrix/matrix.h"
#include "utils/constants.h"
#include "utils/utils.h"
#include "QMspins/spin_ten.h"
#include "QMspins/space_ten.h"

BEGIN_BL_NAMESPACE


class Qua {
	private:
		double Q_;
		double eta_;
		//orientation with respect to another tensor
		double si_;
		double chi_;
		double psi_;
		//the spin number in a series if there is one
		int on_;

		//the A2m rotation bits for moving from PAS to crystal
		Vector<complex> crystalAs;
		double Bfield_;

		rmatrix PAS_;

	public:
		Qua();
		Qua(const Qua &cp):
			Q_(cp.Q_), eta_(cp.eta_),
			si_(cp.si_), chi_(cp.chi_), psi_(cp.psi_),
			on_(cp.on_), crystalAs(cp.crystalAs),Bfield_(cp.Bfield_),
			order(cp.order)
		{} ;

		Qua(double Q, double eta, int on, double si=0, double chi=0, double psi=0, double Bf=0.):
			Q_(Q), eta_(eta), si_(si), chi_(chi), psi_(psi), on_(on),crystalAs(5,0.0), Bfield_(Bf),order(1)

		{
			setCrystalAs();
		};

		~Qua() {};

	//returns a std::string like Q(spin)
		std::string name();

		int order; //use only the second order bits (order==2) or all bits (order==1)
		           //to use second order at all, Bfield > 0
		void reset();
		hmatrix T20;  //the spin matrix, T20
		matrix T21T2m1;  //the second order spin mats
		matrix T22T2m2;  //the sec ord spin mats

		rdmatrix c2; //the rank 2 total spin matrix
		rdmatrix c4; //the rank 4 total spin matrix

		void setSpinMats(SpinSys &A);

		Qua operator=(const Qua &another);
		friend std::ostream& operator<<(std::ostream &otr,const Qua &);

		inline bool is_zero()	{	return Q_==0;}
		inline int &spinON()		{	return on_;		}
		inline int &on()			{	return on_;		}
		inline double &couple()	{	return Q_;		}
		inline double &get_Q()	{	return Q_;		}
		inline double &eta_cop()	{ 	return eta_;	}
		inline double &Q()		{	return Q_;		}
		inline double &eta()		{	return eta_;	}
		inline double &si()		{	return si_;		}
		inline double &alpha()	{	return si_;		}
		inline double &chi()		{	return chi_;	}
		inline double &beta()		{	return chi_;	}
		inline double &psi()		{	return psi_;	}
		inline double &gamma()	{	return psi_;	}
		inline double &Bfield()	{	return Bfield_;	}

		inline int spinON()		const{	return on_;		}
		inline int on()			const{	return on_;		}
		inline double couple()	const{	return Q_;		}
		inline double get_Q()	const{	return Q_;		}
		inline double eta_cop()	const{ 	return eta_;	}
		inline double Q()	const	{	return Q_;		}
		inline double eta()	const	{	return eta_;	}
		inline double si()	const	{	return si_;		}
		inline double alpha()	const	{	return si_;		}
		inline double chi()	const	{	return chi_;	}
		inline double beta() const		{	return chi_;	}
		inline double psi()	const	{	return psi_;	}
		inline double gamma() const	{	return psi_;	}
		inline double Bfield() const	{	return Bfield_;	}

		inline void Q(double in)		{	 Q_=in; setCrystalAs();	}
		inline void eta(double in)		{	 eta_=in; setCrystalAs();	}
		inline void on(int in)			{	 on_=in;	}
		inline void si(double in)		{	 si_=in; setCrystalAs();		}
		inline void alpha(double in)	{	 si_=in; setCrystalAs();		}
		inline void chi(double in)		{	 chi_=in; setCrystalAs();	}
		inline void beta(double in)		{	 chi_=in; setCrystalAs();	}
		inline void psi(double in)		{	 psi_=in; setCrystalAs();	}
		inline void gamma(double in)	{	 psi_=in; setCrystalAs();	}
		inline void Bfield(double in) 	{	Bfield_=in;	}

		//double freq(double a, double b, double g, double p, double t, double s);
		void setCrystalAs();

	/*** Total Hamiltonians ***/
		dmatrix H(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots(phi, theta, gam, alpha, beta);
			return get_H(A, rots);
		}

		dmatrix Hamiltonian(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots(phi, theta, gam, alpha, beta);
			return get_H(A, rots);
		}

		dmatrix H(SpinSys &A,Rotations &therotation)
		{	return get_H(A, therotation);	}

		dmatrix Hamiltonian(SpinSys &A,Rotations &therotation)
		{	return get_H(A, therotation);	}

		dmatrix get_H(SpinSys &,double alpha, double beta, double theta, double phi, double gam=0.0);
		dmatrix get_H(SpinSys &, Rotations &therotation);
		//hmatrix get_H(SpinSys &,double alpha, double beta, double gam, double theta, double phi);
		dmatrix get_HINT(SpinSys &A,double wr,double alpha1, double alpha2,
		                double beta,double theta, double phi, double B_field);

	/*** Second Order hamiltonainas ***/
		dmatrix H2(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots(phi, theta, gam, alpha, beta);
			return get_H2(A, rots);
		}

		dmatrix Hamiltonian2(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots(phi, theta, gam, alpha, beta);
			return get_H2(A, rots);
		}

		dmatrix H2(SpinSys &A,Rotations &therotation)
		{	return get_H2(A, therotation);	}

		dmatrix Hamiltonian2(SpinSys &A,Rotations &therotation)
		{	return get_H2(A, therotation);	}

		dmatrix get_H2(SpinSys &,double alpha, double beta, double theta, double phi, double gam=0.0);
		dmatrix get_H2(SpinSys &, Rotations &therotation);

	/*** first Order Hamiltonain***/
		dmatrix H1(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots(phi, theta, gam, alpha, beta);
			return get_H1(A, rots);
		}

		dmatrix Hamiltonian1(SpinSys &A,double alpha=0, double beta=0, double theta=0, double phi=0, double gam=0.0)
		{
			Rotations rots(phi, theta, gam, alpha, beta);
			return get_H1(A, rots);
		}

		dmatrix H1(SpinSys &A,Rotations &therotation)
		{	return get_H1(A, therotation);	}

		dmatrix Hamiltonian1(SpinSys &A,Rotations &therotation)
		{	return get_H1(A, therotation);	}

		dmatrix get_H1(SpinSys &,double alpha, double beta, double theta, double phi, double gam=0.0);
		dmatrix get_H1(SpinSys &, Rotations &therotation);

/*** parameter setting ***/
		void setParam(const std::string, double);
		double getParam(const std::string);

		bool operator==(const Qua &rhs);
		bool operator!=(const Qua &rhs);

		void display();
		void write(std::ostream &);
};

Vector<Qua> read_qua(const char *filename);
Vector<Qua> read_qua(std::string filename);
Vector<Qua> read_qua(std::ifstream &filename);
Vector<Qua> read_qua(const Vector<std::string> &ss);

END_BL_NAMESPACE


#endif

