


/* blochParamsParticle.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-09-01
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
 	blochParamsParticle.cc-->the basic storage block for the SINGLE spin..
 */

#ifndef _bloch_params_Particle_cc_
#define _bloch_params_Particle_cc_ 1



#include "bloch/blochParams.h"
#include "utils/constants.h"
#include "utils/utils.h"
#include "bloch/Isotope.h"
#include "utils/random.h"


BEGIN_BL_NAMESPACE

//A small class that simply contains the 'slew' of parameters
// nessesary to perform a simple bloch simulation

//double BlochParams::Bo_=(4.7);
//double BlochParams::Ho_=(1./permVac *Bo_);
//double BlochParams::temperature_=(300);


BlochParams<BPoptions::Particle | BPoptions::HighField, double>::BlochParams():
	BlochParams<BPoptions::Density | BPoptions::HighField, double>()
{
	BasicBlochParams<double>::Mo()=1.0;
}

BlochParams<BPoptions::Particle | BPoptions::HighField, double>::BlochParams(std::string spin):
	BlochParams<BPoptions::Density | BPoptions::HighField, double>(spin)
{
	BasicBlochParams<double>::Mo()=1.0;
}


BlochParams<BPoptions::Particle | BPoptions::HighField, double>::
  BlochParams(const BlochParams<BPoptions::Particle | BPoptions::HighField, double> &copy):
	BlochParams<BPoptions::Density | BPoptions::HighField, double>(copy)
{}


void BlochParams<BPoptions::Particle | BPoptions::HighField, double>::calcMo()
{	BasicBlochParams<double>::Mo()=1.0;  }

//this is for use when performing 'dimensionless' integration routines
//the intial magnitization simply becomes 1...or if there are several
//species present in various concentrations

void BlochParams<BPoptions::Particle | BPoptions::HighField, double>::calcMoNorm()
{	BasicBlochParams<double>::Mo()=1.0;	}

void BlochParams<BPoptions::Particle | BPoptions::HighField, double>::Mo(double inMo)
{ BasicBlochParams<double>::Mo()=inMo;	}

void BlochParams<BPoptions::Particle | BPoptions::HighField, double>::Bo(double inBo){
	BasicBlochParams<double>::Bo()=inBo;
	BasicBlochParams<double>::Ho()=1/permVac*BasicBlochParams<double>::Bo();
	calcMo();
	BasicBlochParams<double>::omega()=BasicBlochParams<double>::Bo()*BasicBlochParams<double>::gamma();
}


//Initial Condition Setter..This should be used AFTER !! the 'calcMo()' function call
//this is also the main utility for the 'Bloch' class to set the 'grand' initial
//condition from a list of BlochParams
// 'IC' is the 'Initalondion' class enum
//.. 'HalfFlag' comes from the 'BlochBasic' class where it determins the sign of the
// Mo
void
  BlochParams<BPoptions::Particle | BPoptions::HighField, double>::
setInitialCondition(int IC, int HalfFlag)
{
	Random<UniformRandom<> > angleR(-1.0,1.0);
	switch(IC){
		case InitialCondition::RandomDistribution:
			BasicBlochParams<double>::Mo()=angleR();
			break;
		case InitialCondition::RandomUp:
			angleR.set(0.0, 1.0);
			BasicBlochParams<double>::Mo()=angleR();
			break;
		case InitialCondition::RandomDown:
			angleR.set(-1.0, 0.0);
			BasicBlochParams<double>::Mo()=angleR();
			break;
		case InitialCondition::RandomUpDown:
			angleR.set(0.0, 1.0);
			if(angleR()>0.5){
				BasicBlochParams<double>::Mo()=-BasicBlochParams<double>::Mo();
			}
			break;
		case InitialCondition::AllDown:
			BasicBlochParams<double>::Mo()=-abs(BasicBlochParams<double>::Mo());
			break;
		case InitialCondition::HalfUpHalfDown:
			if(HalfFlag) BasicBlochParams<double>::Mo()=-BasicBlochParams<double>::Mo();
			break;
		case InitialCondition::AllUp: //all are up upon a CalcMo call
		default:
		break;
	}
}


void BlochParams<BPoptions::Particle | BPoptions::HighField, double>::temperature(double Tempin)
{	BasicBlochParams<double>::temperature()=Tempin;  calcMo();	}

void BlochParams<BPoptions::Particle | BPoptions::HighField, double>::moles(double inmoles)
{	BasicBlochParams<double>::moles()=inmoles; calcMo();	}


void
  BlochParams<BPoptions::Particle | BPoptions::HighField, double>::operator=
  (const BlochParams<BPoptions::Particle | BPoptions::HighField, double> &rhs)
{
	if(this==&rhs) return;
	BlochParams<BPoptions::Density | BPoptions::HighField, double>::operator=(rhs);
}

BlochParams<BPoptions::Particle | BPoptions::HighField, double>
  &BlochParams<BPoptions::Particle | BPoptions::HighField, double>::operator=(const std::string &rhs)
{
	BlochParams<BPoptions::Density | BPoptions::HighField, double>::operator=(rhs);
	return *this;
}


bool BlochParams<BPoptions::Particle | BPoptions::HighField, double>::operator==
 (const BlochParams<BPoptions::Particle | BPoptions::HighField, double> &rhs)
{
	if(this==&rhs) return true;
	return BlochParams<BPoptions::Density | BPoptions::HighField, double>::operator!=(rhs);
}

bool BlochParams<BPoptions::Particle | BPoptions::HighField, double>::operator!=
 (const BlochParams<BPoptions::Particle | BPoptions::HighField, double> &rhs)
{
	if(this==&rhs) return false;
	return BlochParams<BPoptions::Density | BPoptions::HighField, double>::operator==(rhs);
}

void BlochParams<BPoptions::Particle | BPoptions::HighField, double>::print(std::ostream &oo) const
{
	oo<<"Paritcle BlochParam: Mo="<<BasicBlochParams<double>::Mo()
	//<<", T1="<<BasicBlochParams<double>::T1()
	 // <<", T2="<<BasicBlochParams<double>::T2()
	 <<", Spin="
	  <<BasicBlochParams<double>::symbol()<<", momentum="<<BasicBlochParams<double>::momentum()<<", gamma="<<BasicBlochParams<double>::gamma();

}

/********** COORDS ****************/
/********** COORDS ****************/
/********** COORDS ****************/

BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::BlochParams():
	BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >()
{
	BasicBlochParams<coord<> >::Mo()=coord<>(0.0,0.0,1.0);
}

BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::BlochParams(std::string spin):
	BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >(spin)
{
	BasicBlochParams<coord<> >::Mo()=coord<>(0.0,0.0,1.0);
}


BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::
  BlochParams(const BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> > &copy):
	BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >(copy)
{}


void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::calcMo()
{
	BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >::calcMo();
	BasicBlochParams<coord<> >::Mo()/=norm(BasicBlochParams<coord<> >::Mo());
}

//this is for use when performing 'dimensionless' integration routines
//the intial magnitization simply becomes 1...or if there are several
//species present in various concentrations

void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::calcMoNorm()
{	calcMo();	}

void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::Mo(const coord<>  &inMo)
{ BasicBlochParams<coord<> >::Mo()=inMo;	}

void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::Bo(const coord<>  &inBo){
	BasicBlochParams<coord<> >::Bo()=inBo;
	BasicBlochParams<coord<> >::Ho()=1/permVac*BasicBlochParams<coord<> >::Bo();
	calcMo();
	BasicBlochParams<coord<> >::omega()=BasicBlochParams<coord<> >::Bo()*BasicBlochParams<coord<> >::gamma();
}

void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::Mo( double  inMo)
{ BasicBlochParams<coord<> >::Mo()=coord<>(0.0,0.0,inMo);	}

void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::Bo( double inBo){
	BasicBlochParams<coord<> >::Bo()=coord<>(0.0,0.0,inBo);
	BasicBlochParams<coord<> >::Ho()=1/permVac*BasicBlochParams<coord<> >::Bo();
	calcMo();
	BasicBlochParams<coord<> >::omega()=BasicBlochParams<coord<> >::Bo()*BasicBlochParams<coord<> >::gamma();
}


void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::temperature(double  Tempin)
{	BasicBlochParams<coord<> >::temperature()=Tempin;  calcMo();	}

void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::moles(double inmoles)
{	BasicBlochParams<coord<> >::moles()=inmoles; calcMo();	}


void
  BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::operator=
  (const BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> > &rhs)
{
	if(this==&rhs) return;
	BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >::operator=(rhs);
}

BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >
  &BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::operator=(const std::string &rhs)
{
	BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >::operator=(rhs);
	return *this;
}


bool BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::operator==
 (const BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> > &rhs)
{
	if(this==&rhs) return true;
	return BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >::operator!=(rhs);
}

bool BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::operator!=
 (const BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> > &rhs)
{
	if(this==&rhs) return false;
	return BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >::operator==(rhs);
}

void BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >::print(std::ostream &oo) const
{
	oo<<"Paritcle BlochParam: Mo="<<BasicBlochParams<coord<> >::Mo()
	//<<", T1="<<BasicBlochParams<double>::T1()
	 // <<", T2="<<BasicBlochParams<double>::T2()
	 <<", Spin="
	  <<BasicBlochParams<coord<> >::symbol()<<", momentum="<<BasicBlochParams<coord<> >::momentum()
	  <<", gamma="<<BasicBlochParams<coord<> >::gamma();

}


END_BL_NAMESPACE
#endif


