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
/* spinsys.cc *****************************************C++ **

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

#ifndef _spinsys_cc_
#define _spinsys_cc_ 1



#include "QMspins/spinsys.h"
#include "container/Vector/Vector.h"
#include "container/matrix/matrix.h"
#include "QMspins/multispinop.h"


BEGIN_BL_NAMESPACE



SpinSys &SpinSys::operator=(const spin_sysPars &pars)
{
	parameters()=(pars);
	totgenerated_=false;
	Ixs_.resize(spins());
	Izs_.resize(spins());
	Iys_.resize(spins());
	Ips_.resize(spins());
	Imis_.resize(spins());
	Ixg.resize(spins());
	Izg.resize(spins());
	Iyg.resize(spins());
	Ipg.resize(spins());
	Imig.resize(spins());
	gen_Tots();
	return *this;
}

SpinSys &SpinSys::operator=(const SpinSys &rhs)
{
	if(this==&rhs) return *this;
	spin_sysPars::operator=(rhs);
	totgenerated_=rhs.totgenerated_;
	Imis_=rhs.Imis_; 	Imig=rhs.Imig;
	Ips_=rhs.Ips_;		Ipg=rhs.Ipg;
	Izs_=rhs.Izs_;		Izg=rhs.Izg;
	Ixs_=rhs.Ixs_;		Ixg=rhs.Ixg;
	Iys_=rhs.Iys_;		Iyg=rhs.Iyg;
	F0_=rhs.F0_;
	Fe_=rhs.Fe_;
	Fmi_=rhs.Fmi_;
	Fp_=rhs.Fp_;
	Fx_=rhs.Fx_;
	Fy_=rhs.Fy_;
	Fz_=rhs.Fz_;
	//GenSpinOps();
	return *this;
}


spin_sysPars &SpinSys::parameters()
{
	return this->Params();
}

void SpinSys::resize(int i)
{
	if(i==spins()) return;
	spin_sysPars::resize(i);
	Ixs_.resize(i); Ixg.resize(i);
	Izs_.resize(i); Izg.resize(i);
	Iys_.resize(i); Iyg.resize(i);
	Ips_.resize(i); Ipg.resize(i);
	Imis_.resize(i); Imig.resize(i);
	if(totgenerated_) gen_Tots();
}

//rimatrix	SpinSys::Ie(int i) const
//{
//	if(!checkSpin(i, false)) SpinOpErr();
//	if(!totgenerated_) Fe_.resize(HS(), HS());
//	return Fe_;
//}


rimatrix	&SpinSys::Ie(int i)
{
	if(!checkSpin(i, false)) SpinOpErr();
	if(!totgenerated_) Fe_.resize(HS(), HS());
	return Fe_;
}

rmatrix		&SpinSys::Imi(int i)
{
	if(!checkSpin(i, false)) SpinOpErr();
	if(Imig[i]) return Imis_[i];
	Imis_[i]=mps_Imi(Params(), i);
	return Imis_[i];
}

rmatrix		&SpinSys::Ip(int i)
{
	if(!checkSpin(i, false)) SpinOpErr();
	if(Ipg[i]) return Ips_[i];
	Ips_[i]= mps_Ip(Params(), i);
	return Ips_[i];
}

rdmatrix	&SpinSys::Iz(int i)
{
	if(!checkSpin(i, false)) SpinOpErr();
	if(Izg[i]) return Izs_[i];
	Izs_[i]= mps_Iz(Params(), i);
	return Izs_[i];
}

hmatrix		&SpinSys::Iy(int i)
{
	if(!checkSpin(i, false)) SpinOpErr();
	if(Iyg[i]) return Iys_[i];
	Iys_[i]= mps_Iy(Params(), i);
	return Iys_[i];
}

smatrix		&SpinSys::Ix(int i)
{
	if(!checkSpin(i, false)) SpinOpErr();
	if(Ixg[i]) return Ixs_[i];
	Ixs_[i]= mps_Ix(Params(), i);
	return Ixs_[i];
}


void SpinSys::gen_Tots()
{
	//F0_.resize(HS(), HS());
	//F0_.fill(ZeroType<double>::zero());
	F0_=rdmatrix(HS(), HS(),ZeroType<double>::zero());
	Fe_.resize(HS(),HS());
	Fx_.resize(HS(), HS());
	Fx_.fill(ZeroType<double>::zero());
	Fz_.resize(HS(), HS());
	Fz_.fill(ZeroType<double>::zero());
	Fmi_.resize(HS(), HS());
	Fmi_.fill(ZeroType<double>::zero());
	Fp_.resize(HS(), HS());
	Fp_.fill(ZeroType<double>::zero());
	Fy_.resize(HS(), HS());
	Fy_.fill(ZeroType<double>::zero());
	for(int i=0;i<spins();i++)
	{
		Fx_+=mps_Ix(Params(), i);
		Fy_+=mps_Iy(Params(), i);
		Fz_+=mps_Iz(Params(), i);
		Fmi_+=mps_Imi(Params(), i);
		Fp_+=mps_Ip(Params(), i);
	}
	totgenerated_=true;
}

void SpinSys::GenSpinOps()
{
	F0_.resize(HS(), HS());
	F0_.fill(ZeroType<double>::zero());
	Fe_.resize(HS(),HS());
	Fx_.resize(HS(), HS());
	Fx_.fill(ZeroType<double>::zero());
	Fz_.resize(HS(), HS());
	Fz_.fill(ZeroType<double>::zero());
	Fmi_.resize(HS(), HS());
	Fmi_.fill(ZeroType<double>::zero());
	Fp_.resize(HS(), HS());
	Fp_.fill(ZeroType<double>::zero());
	Fy_.resize(HS(), HS());
	Fy_.fill(ZeroType<double>::zero());
	for(int i=0;i<spins();i++)
	{
		Ixg(i)=true;
		Ixs_[i]=mps_Ix(Params(), i);
		Fx_+=Ixs_[i];
		Iyg(i)=true;
		Iys_[i]=mps_Iy(Params(), i);
		Fy_+=Iys_[i];
		Izg(i)=true;
		Izs_[i]=mps_Iz(Params(), i);
		Fz_+=Izs_[i];
		Imig(i)=true;
		Imis_[i]=mps_Imi(Params(), i);
		Fmi_+=Imis_[i];
		Ipg(i)=true;
		Ips_[i]=mps_Ip(Params(), i);
		Fp_+=Ips_[i];
	}
	totgenerated_=true;
}

END_BL_NAMESPACE


#endif

