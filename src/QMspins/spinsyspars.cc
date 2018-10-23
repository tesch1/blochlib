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

spin_sysPars.cc --> the list of parameters specfifc to a spin system


**/

#ifndef _SpinSysPars_cc_
#define _SpinSys_cc_ 1

#include "QMspins/spinsyspars.h"
#include "utils/utils.h"
#include <ctype.h>

BEGIN_BL_NAMESPACE


std::string spin_sysPars::DefIso = std::string("1H");	// Set default isotope type to 1H

spin_sysPars::spin_sysPars():
	retIso("1H")
{
	nspins   = 0;
}

spin_sysPars::spin_sysPars(int ns):
	retIso("1H")
{
	nspins   = ns;
	if(ns>0){
		Isotope IsoTmp(retIso);
		for(int i=0; i<ns; i++)	Isotopes.push_back(IsoTmp);
	}
}

void spin_sysPars::resize(int i)
{
	if(nspins==i) return;
	int oldsize=nspins;
	if(oldsize==0) oldsize=1;
	nspins = i;
	Isotopes.resizeAndPreserve(i);
	if(i>0){
		Isotope IsoTmp(retIso);
		for(int j=oldsize; j<nspins; j++)	Isotopes[j]=IsoTmp;
	}
}

spin_sysPars::~spin_sysPars() { }


spin_sysPars& spin_sysPars::operator=(const spin_sysPars& sys)
{
	if(this == &sys) return *this;
	nspins    = sys.nspins;
	Isotopes  = sys.Isotopes;
	return *this;
}

spin_sysPars &spin_sysPars::Params() {	return *this;	}
spin_sysPars spin_sysPars::Params()const {	return *this;	}

bool spin_sysPars::operator==(const spin_sysPars& sys) const
{
	if(this == &sys) return true;
	if(nspins != sys.spins()) return false;
	for(int i=0; i<nspins; i++)	{
		if(Isotopes[i] != sys.Isotopes[i]) return false;
	}
	return true;
}

bool spin_sysPars::operator!=(const spin_sysPars &sys) const
{
	if(nspins != sys.spins()) return true;
	for(int i=0; i<nspins; i++){
		if (Isotopes[i] != sys.Isotopes[i])   return true;
	}

	return false;
}


int spin_sysPars::HS() const
{
	int hs=1;
	for(int i=0; i<nspins; i++){
		hs *= Isotopes[i].HS();
	}
	return hs;
}


double spin_sysPars::qn() const
{
	double q = 0;
	for(int spin=0; spin<nspins; spin++)	q += Isotopes[spin].qn();
	return q;
}

double spin_sysPars::weight() const
{
	double w = 0;
	for(int spin=0; spin<nspins; spin++)	w += Isotopes[spin].weight();
	return w;
}

double spin_sysPars::mass() const
{
	double w = 0;
	for(int spin=0; spin<nspins; spin++)	w += Isotopes[spin].mass();
	return w;
}

void spin_sysPars::isotope(int spin, const std::string& symbol)
{
	isotope(spin,Isotope(symbol));
}


void spin_sysPars::isotope(int spin, const Isotope& Iso)
{
  if(checkSpin(spin))  Isotopes[spin] = Iso;
}


std::string spin_sysPars::momentum() const
{
	double d=qn();
	if (int(d)==d) return itost(int(d));
	else           return itost(int(2*d))+std::string("/2");
}



bool spin_sysPars::homonuclear() const
{
	Isotope iso = isotope(0);
	for (int i=1; i<nspins; i++)	if (isotope(i) != iso){ return false; }
	return true;
}


bool spin_sysPars::heteronuclear() const
{
	Isotope iso = isotope(0);
	for(int i=1; i<nspins; i++)	if (isotope(i) != iso) return true;
	return false;
}


bool spin_sysPars::spinhalf() const
{
	for(int i=0; i<nspins; i++)	if(qn(i) != 0.5) return false;
	return true;
}




int spin_sysPars::isotopes()
{
	int niso=0;
	if(nspins){
		Isotope *Isos; Isos=new Isotope[nspins];
		niso++;
		Isos[0] = isotope(0);
		int i=0, j=0, found=0;
		for(i=1; i<nspins; i++){
			found=0;
			for(j=0; j<niso && !found; j++){
				if(isotope(i) == Isos[j])	found = 1;
			}
			if(!found){
				Isos[niso] = isotope(i);
				niso++;
			}
		}
	}
	return niso;
}

bool spin_sysPars::isotopes(const std::string& I) const
{
	for(int i=0; i<nspins; i++)   if(symbol(i) == I) return true;
	return false;
}

Isotope &spin_sysPars::operator()(int i)
{
	if(checkSpin(i)) return Isotopes.data()[i];
	return retIso;
}

Isotope spin_sysPars::operator()(int i) const
{
	if(checkSpin(i)) return Isotopes.data()[i];
	Isotope loo("1H");
	return loo;
}

Isotope &spin_sysPars::operator[](int i)
{
	if(checkSpin(i)) return Isotopes.data()[i];
	return retIso;

}

Isotope spin_sysPars::operator[](int i) const
{
	if(checkSpin(i)) return Isotopes.data()[i];
	Isotope loo("1H");
	return loo;
}


std::ostream& spin_sysPars::print(std::ostream& ostr) const
{
	std::string s;
	ostr << "Spin     :";
	int i;
	for(i=0; i<nspins; i++) ostr << itost(i,10);
	ostr << "\nIsotope  :";
	for(i=0; i<nspins; i++){
		s = symbol(i);
		ostr << std::string(10-s.length(), ' ') << s;
	}
	ostr << "\nMomentum :";
	for (i=0; i<nspins; i++){
		s = momentum(i);
		ostr << std::string(10-s.length(), ' ') << s;
	}
	ostr << "\n";
	return ostr;
}


std::ostream& operator<< (std::ostream& ostr, const spin_sysPars& sys)
{
	sys.print(ostr);
	return ostr;
}


Vector<std::string> spin_sysPars::printstrings() const{
	Vector<std::string> StrArray;
	int i;
	std::string Stemp, s;

	Stemp = "Spin     :";
	for(i=0; i<nspins; i++) Stemp += itost(i,10);
	StrArray.push_back(Stemp);

	Stemp = "Isotope  :";

	for(i=0; i<nspins; i++) {
		s = symbol(i);
		Stemp += std::string(10-s.length(), ' ') + s;
	}
	StrArray.push_back(Stemp);
	Stemp = "Momentum :";
	for(i=0; i<nspins; i++)	{
		s = momentum(i);
		Stemp += std::string(10-s.length(), ' ') + s;
	}
	StrArray.push_back(Stemp);
	return StrArray;
}

END_BL_NAMESPACE


#endif
