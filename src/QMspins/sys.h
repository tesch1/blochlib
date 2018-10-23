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
 sys.h --> manages the spin system, AND spin interactions (csa, j, dip, qua)
*/

#ifndef _sys_h_
#define _sys_h_ 1

#include "QMspins/csa.h"
#include "QMspins/dip.h"
#include "QMspins/qua.h"
#include "QMspins/j.h"

BEGIN_BL_NAMESPACE


/* the is the top of the various spin parameters
 * it simply contains a 'gamma' spin system and functions to set up the
 * system.  Why devote a whole class to such a little item? you my ask
 * well all the hamiltonians depend on the total spin space, the spin of each
 * and the SpinSys keeps track of those things, so i can generate any operator
 * given the spin...that's why csa.h, dip.h, j.h, qua.h are it's child
 */



class SolidSys : public SpinSys{
	private:

		const void SizeErr();


		hmatrix curhamil;

		double Bfield_;

	public:
		Vector<Csa> csa;
		Vector<Dip> dip;
		Vector<J> jcop;
		Vector<Qua> qua;

		std::string fname;

		bool to_gamma;  //true=do a third power angle
										//false=do NOT """
		bool isdiag; //is the hamil diagonal?
		matrix roeq;		//eq density matrix
		rimatrix III;		//the total sytem identity matrix
		rdmatrix ZZZ;   //the '0' matrix

		SolidSys() : SpinSys()
		{
			theRotations.setPowderAngles(0,0);
			theRotations.setRotorAngles(0,0);
			Bfield_=(0.0);
			isdiag=true;
		};

		SolidSys(const SolidSys &cp);

		SolidSys(int si): SpinSys(si)
		{
			theRotations.setPowderAngles(0,0);
			theRotations.setRotorAngles(0,0);
			Bfield_=(0.0);
			isdiag=true;
		};

		SolidSys(const SpinSys &Ai): SpinSys(Ai)
		{
			theRotations.setPowderAngles(0,0);
			theRotations.setRotorAngles(0,0);
			Bfield_=(400.0);
			isdiag=true;
		};

		SolidSys operator=(const SolidSys &rhs);

		~SolidSys() {};

		SolidSys(const std::string in);
		SolidSys(const char * in);
		SolidSys(const Vector<std::string> &in);

		void read(const std::string in);
		void read(const char *in);
		void read(const char * in, std::ifstream &inf);
		void read(const Vector<std::string> &ss);

		void addDip(Dip lhs);
		void addCsa(Csa lhs);
		void addJ(J lhs);
		void addQ(Qua lhs);

		//inline int size() const { return A.spins(); }
		inline int size() const { return spins(); }

		inline double &Bfield(){	return Bfield_;	}
		inline double &Bo(){	return Bfield_;	}
		inline void SetBo(double Bf){	setBfield(Bf);	}
		void setBfield(double Bf);


	//sets 'Equilibrium den matrix' for the total spin system
		void setRoEQ();
		matrix GetRoEQ(){	return roeq;	}

	//sets 'III' and 'ZZZ'
		void setMats(const SpinSys &A);

	//this little guy runs through all the
	//spin types (Dip , CSA, J, Q) and presets
	//the spin matrices WITHIN THEIR CLASSES
		void setSpinMats();

	//sets the crystall angles inside all of the spins
	//calls 'set_cyrstalAs()' within each spin vector
		void setCrystalAs();

		//the rotation element classes...so we can pass it to the other spin types
		Rotations theRotations;

		void setPowderAngles(double theta, double phi, double gam=0.0);

		void setRotorAngles(double wr, double be,double chi=0.0);

		void setAngles(double wr, double be, double theta, double phi);


	//checks the CSA and Qaudropole
	//spins and sees if 'eta' is present,
	////if it is then we must do an average over
	//a third powder angle..gamma..
		bool check_gamma();

	//given an input std::string will see if that
	//type of interaction and spins
	//are acctually available (i.e. weather or
	//not a D01 (dipole between spin 0 and spin 1)
	//exsits
		bool check_spin_avail(std::string p);

	//these guys set a spin parameter
	//given a char type (D, Q, J, C)
	//a bit (iso, eta, del)
	//a spin number (0...numspins)
	//and the num to set the coupling to
	////D and J require two spin numbers
		void setSpinParam(std::string in, double nm);
		void setSpinParam(char type, std::string bit,int on1, int on2, double num);
		void setSpinParam(char type, std::string bit,int on1, double num);

	//these guys set a spin parameter
	//given a char type (D, Q, J, C)
	//a bit (iso, eta, del)
	//a spin number (0...numspins)
	//and the num to GET the coupling to
	//D and J require two spin numbers
		double getSpinParam(std::string in);
		double getSpinParam(char type, std::string bit,int on1, int on2);
		double getSpinParam(char type, std::string bit,int on1);


	//if the system only contains CSAs and Quads, then
	// the hamiltonian is Diagonal
		inline bool isDiagonal(){	return isdiag;	}

	//Hamiltonains....
		hmatrix &H(double wr , double beta, double th, double phi, double t1=0, double t2=0)
		{	return Hamiltonian(wr, beta, th, phi);	}

		hmatrix &H(double t1, double t2, double wr=0.0)
		{	return Hamiltonian(t1,t2,wr);	}

		hmatrix &Hamiltonian(double t1, double t2, double wr=0.0);
		hmatrix &H();
		hmatrix &Hamiltonian();

		hmatrix &Hamiltonian(double wr , double beta, double th, double phi, double t1=0, double t2=0);

		inline SpinSys *spinSys(){	return this;	}

/*
		inline rdmatrix		F0()const {	return A.F0();	}
		inline rimatrix		Fe()const{	return A.Fe();	};
	//	inline rmatrix		Fmi()const{	return A.Fmi();	};
	//	inline rmatrix		Fp()const{	return A.Fp();	};
	//	inline rdmatrix		Fz()const{	return A.Fz();	};
	//	inline smatrix		Fx()const{	return A.Fx();	};
	//	inline hmatrix		Fy()const{	return A.Fy();	};

		inline rdmatrix		&F0() {	return A.F0();	}
		inline rimatrix		&Fe(){	return A.Fe();	};
		inline rmatrix		&Fmi(){	return A.Fmi();	};
		inline rmatrix		&Fp(){	return A.Fp();	};
		inline rdmatrix		&Fz(){	return A.Fz();	};
		inline smatrix		&Fx(){	return A.Fx();	};
		inline hmatrix		&Fy(){	return A.Fy();	};

		inline rimatrix		Ie(int sp) const{	return A.Ie(sp);	};
	//	inline rmatrix		Imi(int sp)const{	return A.Imi(sp);	} ;
	//	inline rmatrix		Ip(int sp)const {	return A.Ip(sp);	};
	//	inline rdmatrix		Iz(int sp)const {	return A.Iz(sp);	};
	//	inline hmatrix		Iy(int sp)const{	return A.Iy(sp);	};
	//	inline smatrix		Ix(int sp)const{	return A.Ix(sp);	};

		inline rimatrix		&Ie(int sp) {	return A.Ie(sp);	};
		inline rmatrix		&Imi(int sp){	return A.Imi(sp);	} ;
		inline rmatrix		&Ip(int sp) {	return A.Ip(sp);	};
		inline rdmatrix		&Iz(int sp) {	return A.Iz(sp);	};
		inline hmatrix		&Iy(int sp){	return A.Iy(sp);	};
		inline smatrix		&Ix(int sp){	return A.Ix(sp);	};

	//interface functions with SpinSys
		inline std::string symbol(int i){	return A.symbol(i);	}
		inline double gamma(int i){	return A.gamma(i);	}

*/

	//display functions
		void display();
		void print(std::ostream &oo);
		void write(std::string fname);
		void write(std::ostream &otr);

		friend std::ostream& operator<<(std::ostream &otr, SolidSys &out);

};

typedef SolidSys SYS;

END_BL_NAMESPACE


#endif


