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
/* spinsys.h *****************************************C++ **

This combines 'spin_sysPars" and "multispinop" to create an entire spin system

it is meant for small spin systems (i.e. 10 or fewer spins in the systems)
becuase i maintains a list of ALL the spinoperators for each spin
in the TOTAL space so that they do not need to be generated more then once
each spin has Ix, Iy, Iz, Ip, Ie, and Imi so for a 10 spin 1/2 system
each matrix contains 1024x1024 elements (this is not entirely true, as the matrix
class does take shortcuts for special matrix types like identity and diagonal)
but for 10 spins there are ~52 million stored complex numbers (~848 Mb of info) and thus LOTS of
memory...so i recommend this class be used for SMALL spin systems

*/


#ifndef _spinsys_h_
#define _spinsys_h_ 1


#include "QMspins/multispinop.h"
#include "QMspins/spinsyspars.h"
#include "container/Vector/Vector.h"
#include "container/matrix/matrix.h"

BEGIN_BL_NAMESPACE



class SpinSys : public spin_sysPars {
	private:
		//rimatrix			Ies_;		//no need for a vector here...it's always the same we use "Fe"
		Vector<rmatrix>		Imis_;
		Vector<rmatrix>		Ips_;
		Vector<rdmatrix>	Izs_;
		Vector<hmatrix>		Iys_;
		Vector<smatrix>		Ixs_;

		Vector<bool> Imig;
		Vector<bool> Ipg;
		Vector<bool> Izg;
		Vector<bool> Iyg;
		Vector<bool> Ixg;

		bool totgenerated_;

		//total space spinops
		rdmatrix			F0_;		//zero matrix
		rimatrix			Fe_;
		rmatrix				Fmi_;
		rmatrix				Fp_;
		rdmatrix			Fz_;
		hmatrix				Fy_;
		smatrix				Fx_;

		void gen_Tots() ;

		const void SpinOpErr()const{
			std::cerr<<std::endl<<"Error::SpinSys::SpinOp(int)"<<std::endl;
			std::cerr<<" Operator you desire does not exists..."<<std::endl;
			std::cerr<<" death..."<<std::endl;
			exit(-1);
		}


	public:
		SpinSys():
			spin_sysPars(),
			Imis_(), Ips_(), Izs_(), Iys_(), Ixs_(),
			Imig(), Ipg(), Izg(), Iyg(), Ixg(),
			totgenerated_(false),
			F0_(),Fe_(), Fmi_(), Fp_(), Fz_(), Fy_(), Fx_()
		{};

		SpinSys(int nsp):
			spin_sysPars(nsp),
			Imis_(nsp), Ips_(nsp), Izs_(nsp), Iys_(nsp), Ixs_(nsp),
			Imig(nsp,false), Ipg(nsp,false), Izg(nsp,false), Iyg(nsp,false), Ixg(nsp,false),
			totgenerated_(false)

		{};

		SpinSys(const spin_sysPars &in):
			spin_sysPars(in),
			Imis_(in.spins()), Ips_(in.spins()), Izs_(in.spins()), Iys_(in.spins()), Ixs_(in.spins()),
			Imig(in.spins(),false), Ipg(in.spins(),false), Izg(in.spins(),false), Iyg(in.spins(),false), Ixg(in.spins(),false),
			totgenerated_(false),
			F0_(in.HS(),in.HS(), ZeroType<double>::zero()),
			Fe_(in.HS(),in.HS()), Fmi_(in.HS(),in.HS()), Fp_(in.HS(),in.HS()), Fz_(in.HS(),in.HS()),
			Fy_(in.HS(),in.HS()), Fx_(in.HS(),in.HS())
		{};


		SpinSys(const SpinSys &cp):
			spin_sysPars(cp),
			Imis_(cp.Imis_), Ips_(cp.Ips_), Izs_(cp.Izs_), Iys_(cp.Iys_), Ixs_(cp.Ixs_),
			Imig(cp.Imig), Ipg(cp.Ipg), Izg(cp.Izg), Iyg(cp.Iyg), Ixg(cp.Ixg),
			totgenerated_(cp.totgenerated_),
			F0_(cp.F0_),
			Fe_(cp.Fe_), Fmi_(cp.Fmi_), Fp_(cp.Fp_), Fz_(cp.Fz_), Fy_(cp.Fy_), Fx_(cp.Fx_)
		{};

		SpinSys &operator=(const spin_sysPars &pars);

		SpinSys &operator=(const SpinSys &rhs);


		void GenSpinOps();

		void resize(int i);

		inline rdmatrix		F0()const		{	return F0_;	}
		inline rimatrix		Fe()const		{	return Fe_;	}
	//	inline rmatrix		Fmi()const		{	if(totgenerated_){return Fmi_ ;}else{gen_Tots(); return Fmi_;}	}
	//	inline rmatrix		Fp()const		{	if(totgenerated_){return Fp_ ;}else{gen_Tots(); return Fp_;}	}
	//	inline rdmatrix		Fz()const		{	if(totgenerated_){return Fz_ ;}else{gen_Tots(); return Fz_;}	}
	//	inline smatrix		Fx() const		{	if(totgenerated_){return Fx_ ;}else{gen_Tots(); return Fx_;}	}
	//	inline hmatrix		Fy()const		{	if(totgenerated_){return Fy_ ;}else{gen_Tots(); return Fy_;}	}

		inline rdmatrix		&F0()		{	if(totgenerated_){return F0_ ;}else{gen_Tots(); return F0_;	} }
		inline rimatrix		&Fe()		{	if(totgenerated_){return Fe_ ;}else{gen_Tots(); return Fe_;} }
		inline rmatrix		&Fmi()		{	if(totgenerated_){return Fmi_ ;}else{gen_Tots(); return Fmi_;}		}
		inline rmatrix		&Fp()		{	if(totgenerated_){return Fp_ ;}else{gen_Tots(); return Fp_;}	}
		inline rdmatrix		&Fz()		{	if(totgenerated_){return Fz_ ;}else{gen_Tots(); return Fz_;} }
		inline smatrix		&Fx() 		{	if(totgenerated_){return Fx_ ;}else{gen_Tots(); return Fx_;}	}
		inline hmatrix		&Fy()		{	if(totgenerated_){return Fy_ ;}else{gen_Tots(); return Fy_;}	}


	//	rimatrix	Ie(int sp) const;
	//	rmatrix		Imi(int sp) const;
	//	rmatrix		Ip(int sp) const;
	//	rdmatrix	Iz(int sp) const;
	//	hmatrix		Iy(int sp) const;
	//	smatrix		Ix(int sp) const;

		rimatrix	&Ie(int sp);
		rmatrix		&Imi(int sp);
		rmatrix		&Ip(int sp);
		rdmatrix	&Iz(int sp);
		hmatrix		&Iy(int sp);
		smatrix		&Ix(int sp);



		spin_sysPars &parameters() ;
		spin_sysPars parameters() const;
};


//typedef SpinSys SpinSys;


END_BL_NAMESPACE


#endif


