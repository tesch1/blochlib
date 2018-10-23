/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-28-01
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
	space_ten.h--> spacial spherical tensors up to rank 2
*/


#ifndef _space_ten_h_
#define _space_ten_h_

#include "container/matrix/matrix.h"
#include "container/Vector/Vector.h"
#include "utils/constants.h"
#include "utils/utils.h"

BEGIN_BL_NAMESPACE


complex A1(double,double, int);
complex A2(double, double, int);
complex wigner1(double, double,double, int, int);
complex wigner2(double, double, double,int, int);

double dnm2(double beta, int n, int m);
double dnm4(double beta, int n, int m);
double dnml(double beta, int n, int m, int l);

double A_Q_pas(double eta, int m);
double A_CSA_pas(double iso, double del, double eta,int l, int m);
double A_Dip_pas(double dip, int m);


rmatrix rotationMatrix3D(double a, double b, double g);
rmatrix rotationMatrix3D(double a, double b);

//this sets the powder averaging wigners...if we are going to use the
//g-compute method, then gamma is taken care of in the calculation
//at least for the fid part, but for the rest of things, gamma is not
//taken care of and must be included directly

class Rotations{
	public:
		double phi;
		double theta;
		double gamma;
		double alpha;
		double beta;
		double chi;
		matrix powderAs;

		rmatrix cartPowder;
		rmatrix cartSpin;

		bool rotate;

		void setPowderAngles(double phi, double theta, double gamma=0.);
		void setPowderAs(double phi, double theta, double gamma=0.);

		enum{Spherical=0, Cartesian=1, All=2};

		int RotationType;

//this little vector is the final rotation bit taking a crystal into the
//lab frame...becuase we are in the high field approx, we only want the final
//A20 bit, so this saves us an entire matrix multiplication an reduces it
//to the A20=Sum[D0m]*B2m...(the 'spinnerAs=D0m') also the 'gamma' for this rotation is zero
//(or rather it dos not enter into the calculation) becuase we are
//in the high field world
		Vector<complex> spinnerAs;
		void setRotorAngles(double alpha, double beta, double cin=0.0);
		void setSpinnerAs(double alpha, double beta, double cin=0.0);

	//this is a special function that is only needed
	//for the second order quadrupoles and does not need to
	//be calculated unless it is needed..thus it is given its own function
	//it calculates the entire spinner wigner matrix...
		matrix spinnerAsMat;
		void set_spinnerAsMat(double ain, double bin, double cin=0.0);
		void setSpinnerAsMat(double ain, double bin, double cin=0.0);

		Rotations();
		Rotations(int type);
		Rotations(double phi, double theta, double gamma=0, double alpha=0, double beta=0);
		Rotations(matrix &pAs, Vector<complex> &sAs);
		Rotations(const Rotations &cp);

		complex getPow(int i, int j);
		complex getSpn(int i);

		Rotations &operator=(const Rotations &rhs);

};


END_BL_NAMESPACE


#endif


