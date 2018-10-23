

/* bloch_int_demag.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-8-01
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
 	bloch_int_demagmodulated.cc-->Demagnitizing Field interactions for the Bloch equations
 */


#ifndef _bloch_interatction_Modulated_demag_field_cc_
#define _bloch_interatction_Modulated_demag_field_cc_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"
#include "bloch/bloch_int_demagmodulate.h"

BEGIN_BL_NAMESPACE

/* dipole-dipole **Demagnetizing Field** coupling in the classical sence
	requires only the total magenization, the 'gradient' angle
	and the 'time factor' (based typically on td=1/muo*gamma*Mo)

		dm/dt=gamma (M x Bd)
		where Bd_i = 1/td*Del*[ mz-<Mz> + 1/3(M-<M>) ]

		where Del=[3*(s.z)^2 -1]

		where s is the direction of the magnetization modulation

		(i.e. the direction the gradient has been applied)

		This interaction assumes that



*/




ModulatedDemagField::ModulatedDemagField():
	i_on(true), td_(0), s_(0,0,1), delFact(1.0)
{}

ModulatedDemagField::ModulatedDemagField(double td, coord<> s):
	 i_on(true),  td_(td), s_(s/norm(s))
{
	double cosTH=s_.z();
	delFact=(3.0*cosTH*cosTH-1.0)/2.0;
}

ModulatedDemagField::ModulatedDemagField(const ModulatedDemagField &cp):
	  i_on(cp.i_on),  td_(cp.td_), s_(cp.s_), delFact(cp.delFact)
{}

ModulatedDemagField &ModulatedDemagField::operator=(const ModulatedDemagField &cp)
{
	if(&cp==this) return *this;
	td_=cp.td_;
	s_=cp.s_;
	delFact=cp.delFact;
	return *this;
}

ModulatedDemagField::~ModulatedDemagField()
{}

void ModulatedDemagField::setTd(double in){	td_=in;	}
void ModulatedDemagField::setDirection(const coord<> &in){	s_=in/norm(in);	delFact=(3.0*s_.z()*s_.z()-1.0)/2.0;}

double ModulatedDemagField::td() const {	return td_;	}
coord<> ModulatedDemagField::direction() const {	return s_;}

void ModulatedDemagField::off()	{	i_on=false;	}
void ModulatedDemagField::on()	{	i_on=true;	}


std::ostream &operator<<(std::ostream &oo, ModulatedDemagField &out)
{
	oo<<"ModulatedDemagField Interactions"<<std::endl;
	return oo;
}


END_BL_NAMESPACE


#endif


