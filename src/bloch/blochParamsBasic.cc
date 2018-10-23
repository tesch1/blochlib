


/* blochParams.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-30-01
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
 	blochParamsBasic.cc-->the Base class forthe basic storage blocks
 	maitains the things that are always true for a Bloch Parameter in
 	any circustance..
 */



#ifndef _bloch_params_Basic_cc_
#define _bloch_params_Basic_cc_ 1


#include "bloch/Isotope.h"
#include "container/matrix/matrix.h"
#include "bloch/blochParamsBasic.h"

BEGIN_BL_NAMESPACE


//The BASIC BLoch Param....those things that do not change regardless
//of fields, or prticle/density parts

void BasicBlochParams<double >::operator=(const BasicBlochParams<double> &rhs)
{
	if(&rhs==this) return;
	Isotope::operator=(rhs);
	omega_=rhs.omega_;		//Zeeman frequency
	moles_=rhs.moles_;		//number of moles of sample
	Mo_=rhs.Mo_;
	Bo_=rhs.Bo_;			//Magnetic firld strength in TESLA
	Ho_=rhs.Ho_ ;		// = 1/permVac Bo
	temperature_=rhs.temperature_;	//temperature of the sample (Kelvin)
}


BasicBlochParams<double> &BasicBlochParams<double>::operator=(const std::string &rhs)
{
	Isotope::operator=(rhs);
	moles_=1.0;		//number of moles of sample
	Mo_=1.0;
	Bo_=4.7;			//Magnetic firld strength in TESLA
	Ho_=1./permVac *Bo_ ;		// = 1/permVac Bo
	omega_=gamma()*Bo_;		//Zeeman frequency
	temperature_=300.0;	//temperature of the sample (Kelvin)
	return *this;
}

bool BasicBlochParams<double>::operator==(const BasicBlochParams<double> &rhs)
{
	if(this==&rhs) return true;
	if(Mo_==rhs.Mo_) return false;
	if(moles_!=rhs.moles_) return false;
	if(Isotope::operator!=(rhs)) return false;
	if(Bo_!=rhs.Bo_) return false;
	if(Ho_!=rhs.Ho_) return false;
	if(omega_!=rhs.omega_) return false;
	if(temperature_!=rhs.temperature_) return false;
	return true;
}

bool BasicBlochParams<double>::operator!=(const BasicBlochParams<double> &rhs)
{
	if(this==&rhs) return false;
	if(Isotope::operator==(rhs)) return false;
	if(Bo_==rhs.Bo_) return false;
	if(Mo_==rhs.Mo_) return false;
	if(Ho_==rhs.Ho_) return false;
	if(omega_==rhs.omega_) return false;
	if(temperature_==rhs.temperature_) return false;
	if(moles_==rhs.moles_) return false;
	return true;
}



/*********** COORD<> ****************/
void BasicBlochParams<coord<> >::operator=(const BasicBlochParams<coord<> > &rhs)
{
	if(&rhs==this) return;
	Isotope::operator=(rhs);
	omega_=rhs.omega_;		//Zeeman frequency
	moles_=rhs.moles_;		//number of moles of sample
	Mo_=rhs.Mo_;
	Bo_=rhs.Bo_;			//Magnetic firld strength in TESLA
	Ho_=rhs.Ho_ ;		// = 1/permVac Bo
	temperature_=rhs.temperature_;	//temperature of the sample (Kelvin)
}


BasicBlochParams<coord<> > &BasicBlochParams<coord<> >::operator=(const std::string &rhs)
{
	Isotope::operator=(rhs);
	moles_=1.0;		//number of moles of sample
	Mo_(0.0,0.0,1.0);
	Bo_(0.0,0.0,4.7);			//Magnetic firld strength in TESLA
	Ho_=1./permVac *Bo_ ;		// = 1/permVac Bo
	omega_=gamma()*Bo_;		//Zeeman frequency
	temperature_=300.0;	//temperature of the sample (Kelvin)
	return *this;
}

bool BasicBlochParams<coord<> >::operator==(const BasicBlochParams<coord<> > &rhs)
{
	if(this==&rhs) return true;
	if(moles_!=rhs.moles_) return false;
	if(Isotope::operator!=(rhs)) return false;
	if(Bo_.x()!=rhs.Bo_.x()) return false;
	if(Bo_.y()!=rhs.Bo_.y()) return false;
	if(Bo_.z()!=rhs.Bo_.z()) return false;

	if(Mo_.x()!=rhs.Mo_.x()) return false;
	if(Mo_.y()!=rhs.Mo_.y()) return false;
	if(Mo_.z()!=rhs.Mo_.z()) return false;

	if(Ho_.x()!=rhs.Ho_.x()) return false;
	if(Ho_.y()!=rhs.Ho_.y()) return false;
	if(Ho_.z()!=rhs.Ho_.z()) return false;

	if(omega_.x()!=rhs.omega_.x()) return false;
	if(omega_.y()!=rhs.omega_.y()) return false;
	if(omega_.z()!=rhs.omega_.z()) return false;

	if(temperature_!=rhs.temperature_) return false;
	return true;
}

bool BasicBlochParams<coord<> >::operator!=(const BasicBlochParams<coord<> > &rhs)
{
	if(this==&rhs) return false;
	if(Isotope::operator==(rhs)) return false;
	if(moles_==rhs.moles_) return false;

	if(Bo_.x()==rhs.Bo_.x()) return false;
	if(Bo_.y()==rhs.Bo_.y()) return false;
	if(Bo_.z()==rhs.Bo_.z()) return false;

	if(Mo_.x()==rhs.Mo_.x()) return false;
	if(Mo_.y()==rhs.Mo_.y()) return false;
	if(Mo_.z()==rhs.Mo_.z()) return false;

	if(Ho_.x()==rhs.Ho_.x()) return false;
	if(Ho_.y()==rhs.Ho_.y()) return false;
	if(Ho_.z()==rhs.Ho_.z()) return false;

	if(omega_.x()==rhs.omega_.x()) return false;
	if(omega_.y()==rhs.omega_.y()) return false;
	if(omega_.z()==rhs.omega_.z()) return false;

	if(temperature_==rhs.temperature_) return false;
	return true;
}



END_BL_NAMESPACE


#endif

