

/* bloch_interac_meth.h ********/


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
 	bloch_interac_meth.h-->axilliary interactions for the Bloch equations METHODS
 */


#ifndef _bloch_interatction_meth_h_
#define _bloch_interatction_meht_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"

BEGIN_BL_NAMESPACE

// CLASSES::
//	class Interactions
//	  --> there can be up to 6 different interactions in the template argument
//  class BulkSus
//    --> a single interation simply the bulk susceptibilty of the sample...w
//	      which introduces an extra offset term depending on the amount of Mz present
//	class RadDamp
//    --> Radation damping the effectve back reaction field from the total magnetic field
//        present in the coil
//	class DemagField
//	  --> Dipolar field calculation...this is a beast of a caluclation as each spins
//	      requires the entire list
//	class DipoleDipole
//	  --> Dipola-Dipole field calculation...this is a beast of a caluclation as each spins
//	      requires the entire list


/* here are the remaing interactions that are passed as a template parameter to the
	'Bloch' class....There is also a class named "Interactions" that simply act as a
	transfer class for each interactions type
*/


/************ 1 Interaction *********/


template<	class Int_1>
template<class ParamIter>
inline rmatrix
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	return int1_->jacobian(t,M,pars,totM);
}

template<	class Int_1>
template<class ParamIter>
inline void
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1->jacobian(out,t,M,pars,totM);
}


template<	class Int_1>
template<class PList>
inline void
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
}


template<	class Int_1>
template<class PList>
inline void
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1>
template<class PList>
inline void
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
}

template<	class Int_1>
template<class PList>
inline coord<>
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1>
template<class PList>
inline coord<>
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars);
}

template<	class Int_1>
template<class PList>
inline rmatrix
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM);
}

template<	class Int_1>
template<class PList>
inline rmatrix
  Interactions<Int_1, BasicInter ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM);
}

/**********8 2 Interactions *********/

template<	class Int_1, class Int_2>
template<class ParamIter>
inline rmatrix
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	 //cout<<int1_->jacobian(t,M,pars,totM)+int2_->jacobian(t,M,pars,totM)<<endl;;
	return int1_->jacobian(t,M,pars,totM)+int2_->jacobian(t,M,pars,totM);
}

template<	class Int_1, class Int_2>
template<class ParamIter>
inline void
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1_->jacobian(out,t,M,pars,totM);
	int2_->jacobian(out,t,M,pars,totM);
}

template<	class Int_1, class Int_2>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
	int2_->preCalculate(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
}

template<	class Int_1, class Int_2>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM)+
			int2_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars)+
			int2_->function(t, M, dMdt, pars);
}

template<	class Int_1, class Int_2>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM)+
	 		int2_->magneticField(t, M, pars, totM);
}

template<	class Int_1, class Int_2>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,BasicInter,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM)+
	 int2_->evolutionMatrix(t, M, pars, totM);
}


/********** 3 Interactions ***********/


template<	class Int_1, class Int_2, class Int_3>
template<class ParamIter>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	return int1_->jacobian(t,M,pars,totM)
			+int2_->jacobian(t,M,pars,totM)
			+int3_->jacobian(t,M,pars,totM);
}

template<	class Int_1, class Int_2, class Int_3>
template<class ParamIter>
inline void
   Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1_->jacobian(out,t,M,pars,totM);
	int2_->jacobian(out,t,M,pars,totM);
	int3_->jacobian(out,t,M,pars,totM);
}

template<	class Int_1, class Int_2, class Int_3>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
	int2_->preCalculate(t, M, dMdt, pars, totM);
	int3_->preCalculate(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
	int3_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
	int3_->function(t, M, dMdt, pars);
}


template<	class Int_1, class Int_2, class Int_3>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM)+
			int2_->function(t, M, dMdt, pars, totM)+
			int3_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars)+
			int2_->function(t, M, dMdt, pars)+
			int3_->function(t, M, dMdt, pars);
}

/*
template<	class Int_1, class Int_2, class Int_3>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
	int3_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
	int3_->function(t, M, dMdt, pars);
}
*/

template<	class Int_1, class Int_2,class Int_3>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM)+
	 		int2_->magneticField(t, M, pars, totM)+
	 		 int3_->magneticField(t, M, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               BasicInter,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM)+
	 int2_->evolutionMatrix(t, M, pars, totM)+
	  int3_->evolutionMatrix(t, M, pars, totM);
}


/******** 4 Interactions.... **********/

template<	class Int_1, class Int_2, class Int_3, class Int_4>
template<class ParamIter>
inline rmatrix
   Interactions<Int_1, Int_2 ,Int_3,
               Int_4,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	return int1_->jacobian(t,M,pars,totM)
			+int2_->jacobian(t,M,pars,totM)
			+int3_->jacobian(t,M,pars,totM)
			+int4_->jacobian(t,M,pars,totM);
}

template<	class Int_1, class Int_2, class Int_3, class Int_4>
template<class ParamIter>
inline void
 Interactions<Int_1, Int_2 ,Int_3,
               Int_4,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1_->jacobian(out,t,M,pars,totM);
	int2_->jacobian(out,t,M,pars,totM);
	int3_->jacobian(out,t,M,pars,totM);
	int4_->jacobian(out,t,M,pars,totM);
}

template<	class Int_1, class Int_2, class Int_3, class Int_4>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
	int2_->preCalculate(t, M, dMdt, pars, totM);
	int3_->preCalculate(t, M, dMdt, pars, totM);
	int4_->preCalculate(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3, class Int_4>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
	int3_->function(t, M, dMdt, pars, totM);
	int4_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3, class Int_4>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
	int3_->function(t, M, dMdt, pars);
	int4_->function(t, M, dMdt, pars);
}

template<	class Int_1, class Int_2, class Int_3, class Int_4>
template<class PList>
inline coord<>
 Interactions<Int_1, Int_2 ,Int_3,
                Int_4,BasicInter,BasicInter,
                BasicInter,BasicInter,BasicInter>::
  function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM)+
			int2_->function(t, M, dMdt, pars, totM)+
			int3_->function(t, M, dMdt, pars, totM)+
			int4_->function(t, M, dMdt, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3, class Int_4>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars)+
			int2_->function(t, M, dMdt, pars)+
			int3_->function(t, M, dMdt, pars)+
			int4_->function(t, M, dMdt, pars);
}


template<	class Int_1, class Int_2,class Int_3, class Int_4>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM)+
	 		int2_->magneticField(t, M, pars, totM)+
	 		 int3_->magneticField(t, M, pars, totM)+
	 		 int4_->magneticField(t, M, pars, totM);
}

template<	class Int_1, class Int_2, class Int_3, class Int_4>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,BasicInter,BasicInter,
               BasicInter,BasicInter,BasicInter>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM)+
	 int2_->evolutionMatrix(t, M, pars, totM)+
	  int3_->evolutionMatrix(t, M, pars, totM)+
	  int4_->evolutionMatrix(t, M, pars, totM);
}



/******** 5 Interactions.... **********/

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
template<class ParamIter>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
 jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	return int1_->jacobian(t,M,pars,totM)
			+int2_->jacobian(t,M,pars,totM)
			+int3_->jacobian(t,M,pars,totM)
			+int4_->jacobian(t,M,pars,totM)
			+int5_->jacobian(t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
template<class ParamIter>
inline void
   Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1_->jacobian(out,t,M,pars,totM);
	int2_->jacobian(out,t,M,pars,totM);
	int3_->jacobian(out,t,M,pars,totM);
	int4_->jacobian(out,t,M,pars,totM);
	int5_->jacobian(out,t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
	int2_->preCalculate(t, M, dMdt, pars, totM);
	int3_->preCalculate(t, M, dMdt, pars, totM);
	int4_->preCalculate(t, M, dMdt, pars, totM);
	int5_->preCalculate(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
	int3_->function(t, M, dMdt, pars, totM);
	int4_->function(t, M, dMdt, pars, totM);
	int5_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
	int3_->function(t, M, dMdt, pars);
	int4_->function(t, M, dMdt, pars);
	int5_->function(t, M, dMdt, pars);

}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
template<class PList>
inline coord<>
   Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM)+
			int2_->function(t, M, dMdt, pars, totM)+
			int3_->function(t, M, dMdt, pars, totM)+
			int4_->function(t, M, dMdt, pars, totM)+
			int5_->function(t, M, dMdt, pars, totM);

}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars)+
			int2_->function(t, M, dMdt, pars)+
			int3_->function(t, M, dMdt, pars)+
			int4_->function(t, M, dMdt, pars)+
			int5_->function(t, M, dMdt, pars);

}


template<class Int_1, class Int_2,class Int_3, class Int_4, class Int_5>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM)+
	 		int2_->magneticField(t, M, pars, totM)+
	 		 int3_->magneticField(t, M, pars, totM)+
	 		 int4_->magneticField(t, M, pars, totM)+
	 		 int5_->magneticField(t, M, pars, totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,BasicInter,
               BasicInter,BasicInter,BasicInter>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM)+
	 int2_->evolutionMatrix(t, M, pars, totM)+
	  int3_->evolutionMatrix(t, M, pars, totM)+
	  int4_->evolutionMatrix(t, M, pars, totM)+
	  int5_->evolutionMatrix(t, M, pars, totM);
}


/******** 6 Interactions.... **********/

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5, class Int_6>
template<class ParamIter>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
 jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	return int1_->jacobian(t,M,pars,totM)
			+int2_->jacobian(t,M,pars,totM)
			+int3_->jacobian(t,M,pars,totM)
			+int4_->jacobian(t,M,pars,totM)
			+int5_->jacobian(t,M,pars,totM)
			+int6_->jacobian(t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5, class Int_6>
template<class ParamIter>
inline void
   Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1_->jacobian(out,t,M,pars,totM);
	int2_->jacobian(out,t,M,pars,totM);
	int3_->jacobian(out,t,M,pars,totM);
	int4_->jacobian(out,t,M,pars,totM);
	int5_->jacobian(out,t,M,pars,totM);
	int6_->jacobian(out,t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5, class Int_6>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
	int2_->preCalculate(t, M, dMdt, pars, totM);
	int3_->preCalculate(t, M, dMdt, pars, totM);
	int4_->preCalculate(t, M, dMdt, pars, totM);
	int5_->preCalculate(t, M, dMdt, pars, totM);
	int6_->preCalculate(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5, class Int_6>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
	int3_->function(t, M, dMdt, pars, totM);
	int4_->function(t, M, dMdt, pars, totM);
	int5_->function(t, M, dMdt, pars, totM);
	int6_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5, class Int_6>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
	int3_->function(t, M, dMdt, pars);
	int4_->function(t, M, dMdt, pars);
	int5_->function(t, M, dMdt, pars);
	int6_->function(t, M, dMdt, pars);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5, class Int_6>
template<class PList>
inline coord<>
   Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM)+
			int2_->function(t, M, dMdt, pars, totM)+
			int3_->function(t, M, dMdt, pars, totM)+
			int4_->function(t, M, dMdt, pars, totM)+
			int5_->function(t, M, dMdt, pars, totM)+
			int6_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5, class Int_6>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars);
			int2_->function(t, M, dMdt, pars)+
			int3_->function(t, M, dMdt, pars)+
			int4_->function(t, M, dMdt, pars)+
			int5_->function(t, M, dMdt, pars)+
			int6_->function(t, M, dMdt, pars);
}


template<class Int_1, class Int_2,class Int_3, class Int_4, class Int_5, class Int_6>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM)+
	 		int2_->magneticField(t, M, pars, totM)+
	 		 int3_->magneticField(t, M, pars, totM)+
	 		 int4_->magneticField(t, M, pars, totM)+
	 		 int5_->magneticField(t, M, pars, totM)+
	 		 int6_->magneticField(t, M, pars, totM);
}

template<class Int_1, class Int_2, class Int_3, class Int_4, class Int_5, class Int_6>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               BasicInter,BasicInter,BasicInter>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM)+
	 int2_->evolutionMatrix(t, M, pars, totM)+
	  int3_->evolutionMatrix(t, M, pars, totM)+
	  int4_->evolutionMatrix(t, M, pars, totM)+
	  int5_->evolutionMatrix(t, M, pars, totM)+
	  int6_->evolutionMatrix(t, M, pars, totM);

}


/******** 7 Interactions.... **********/

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class ParamIter>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
 jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	return int1_->jacobian(t,M,pars,totM)
			+int2_->jacobian(t,M,pars,totM)
			+int3_->jacobian(t,M,pars,totM)
			+int4_->jacobian(t,M,pars,totM)
			+int5_->jacobian(t,M,pars,totM)
			+int6_->jacobian(t,M,pars,totM)
			+int7_->jacobian(t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class ParamIter>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1_->jacobian(out,t,M,pars,totM);
	int2_->jacobian(out,t,M,pars,totM);
	int3_->jacobian(out,t,M,pars,totM);
	int4_->jacobian(out,t,M,pars,totM);
	int5_->jacobian(out,t,M,pars,totM);
	int6_->jacobian(out,t,M,pars,totM);
	int7_->jacobian(out,t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
	int2_->preCalculate(t, M, dMdt, pars, totM);
	int3_->preCalculate(t, M, dMdt, pars, totM);
	int4_->preCalculate(t, M, dMdt, pars, totM);
	int5_->preCalculate(t, M, dMdt, pars, totM);
	int6_->preCalculate(t, M, dMdt, pars, totM);
	int7_->preCalculate(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
	int3_->function(t, M, dMdt, pars, totM);
	int4_->function(t, M, dMdt, pars, totM);
	int5_->function(t, M, dMdt, pars, totM);
	int6_->function(t, M, dMdt, pars, totM);
	int7_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
	int3_->function(t, M, dMdt, pars);
	int4_->function(t, M, dMdt, pars);
	int5_->function(t, M, dMdt, pars);
	int6_->function(t, M, dMdt, pars);
	int7_->function(t, M, dMdt, pars);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM)+
			int2_->function(t, M, dMdt, pars, totM)+
			int3_->function(t, M, dMdt, pars, totM)+
			int4_->function(t, M, dMdt, pars, totM)+
			int5_->function(t, M, dMdt, pars, totM)+
			int6_->function(t, M, dMdt, pars, totM)+
			int7_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars);
			int2_->function(t, M, dMdt, pars)+
			int3_->function(t, M, dMdt, pars)+
			int4_->function(t, M, dMdt, pars)+
			int5_->function(t, M, dMdt, pars)+
			int6_->function(t, M, dMdt, pars)+
			int7_->function(t, M, dMdt, pars);
}


template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM)+
	 		int2_->magneticField(t, M, pars, totM)+
	 		 int3_->magneticField(t, M, pars, totM)+
	 		 int4_->magneticField(t, M, pars, totM)+
	 		 int5_->magneticField(t, M, pars, totM)+
	 		 int6_->magneticField(t, M, pars, totM)+
	 		 int7_->magneticField(t, M, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,BasicInter,BasicInter>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM)+
	 int2_->evolutionMatrix(t, M, pars, totM)+
	  int3_->evolutionMatrix(t, M, pars, totM)+
	  int4_->evolutionMatrix(t, M, pars, totM)+
	  int5_->evolutionMatrix(t, M, pars, totM)+
	  int6_->evolutionMatrix(t, M, pars, totM)+
	  int7_->evolutionMatrix(t, M, pars, totM);

}


/******** 8 Interactions.... **********/

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class ParamIter>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
 jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	return int1_->jacobian(t,M,pars,totM)
			+int2_->jacobian(t,M,pars,totM)
			+int3_->jacobian(t,M,pars,totM)
			+int4_->jacobian(t,M,pars,totM)
			+int5_->jacobian(t,M,pars,totM)
			+int6_->jacobian(t,M,pars,totM)
			+int7_->jacobian(t,M,pars,totM)
			+int8_->jacobian(t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class ParamIter>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1_->jacobian(out,t,M,pars,totM);
	int2_->jacobian(out,t,M,pars,totM);
	int3_->jacobian(out,t,M,pars,totM);
	int4_->jacobian(out,t,M,pars,totM);
	int5_->jacobian(out,t,M,pars,totM);
	int6_->jacobian(out,t,M,pars,totM);
	int7_->jacobian(out,t,M,pars,totM);
	int8_->jacobian(out,t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
	int2_->preCalculate(t, M, dMdt, pars, totM);
	int3_->preCalculate(t, M, dMdt, pars, totM);
	int4_->preCalculate(t, M, dMdt, pars, totM);
	int5_->preCalculate(t, M, dMdt, pars, totM);
	int6_->preCalculate(t, M, dMdt, pars, totM);
	int7_->preCalculate(t, M, dMdt, pars, totM);
	int8_->preCalculate(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
	int3_->function(t, M, dMdt, pars, totM);
	int4_->function(t, M, dMdt, pars, totM);
	int5_->function(t, M, dMdt, pars, totM);
	int6_->function(t, M, dMdt, pars, totM);
	int7_->function(t, M, dMdt, pars, totM);
	int8_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
	int3_->function(t, M, dMdt, pars);
	int4_->function(t, M, dMdt, pars);
	int5_->function(t, M, dMdt, pars);
	int6_->function(t, M, dMdt, pars);
	int7_->function(t, M, dMdt, pars);
	int8_->function(t, M, dMdt, pars);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM)+
			int2_->function(t, M, dMdt, pars, totM)+
			int3_->function(t, M, dMdt, pars, totM)+
			int4_->function(t, M, dMdt, pars, totM)+
			int5_->function(t, M, dMdt, pars, totM)+
			int6_->function(t, M, dMdt, pars, totM)+
			int7_->function(t, M, dMdt, pars, totM)+
			int8_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars);
			int2_->function(t, M, dMdt, pars)+
			int3_->function(t, M, dMdt, pars)+
			int4_->function(t, M, dMdt, pars)+
			int5_->function(t, M, dMdt, pars)+
			int6_->function(t, M, dMdt, pars)+
			int7_->function(t, M, dMdt, pars)+
			int8_->function(t, M, dMdt, pars);
}


template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM)+
	 		int2_->magneticField(t, M, pars, totM)+
	 		 int3_->magneticField(t, M, pars, totM)+
	 		 int4_->magneticField(t, M, pars, totM)+
	 		 int5_->magneticField(t, M, pars, totM)+
	 		 int6_->magneticField(t, M, pars, totM)+
	 		 int7_->magneticField(t, M, pars, totM)+
	 		 int8_->magneticField(t, M, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,BasicInter>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM)+
	 int2_->evolutionMatrix(t, M, pars, totM)+
	  int3_->evolutionMatrix(t, M, pars, totM)+
	  int4_->evolutionMatrix(t, M, pars, totM)+
	  int5_->evolutionMatrix(t, M, pars, totM)+
	  int6_->evolutionMatrix(t, M, pars, totM)+
	  int7_->evolutionMatrix(t, M, pars, totM)+
	  int8_->evolutionMatrix(t, M, pars, totM);

}

/******** 9 Interactions.... **********/

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class ParamIter>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
 jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	return int1_->jacobian(t,M,pars,totM)
			+int2_->jacobian(t,M,pars,totM)
			+int3_->jacobian(t,M,pars,totM)
			+int4_->jacobian(t,M,pars,totM)
			+int5_->jacobian(t,M,pars,totM)
			+int6_->jacobian(t,M,pars,totM)
			+int7_->jacobian(t,M,pars,totM)
			+int8_->jacobian(t,M,pars,totM)
			+int9_->jacobian(t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class ParamIter>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	int1_->jacobian(out,t,M,pars,totM);
	int2_->jacobian(out,t,M,pars,totM);
	int3_->jacobian(out,t,M,pars,totM);
	int4_->jacobian(out,t,M,pars,totM);
	int5_->jacobian(out,t,M,pars,totM);
	int6_->jacobian(out,t,M,pars,totM);
	int7_->jacobian(out,t,M,pars,totM);
	int8_->jacobian(out,t,M,pars,totM);
	int9_->jacobian(out,t,M,pars,totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->preCalculate(t, M, dMdt, pars, totM);
	int2_->preCalculate(t, M, dMdt, pars, totM);
	int3_->preCalculate(t, M, dMdt, pars, totM);
	int4_->preCalculate(t, M, dMdt, pars, totM);
	int5_->preCalculate(t, M, dMdt, pars, totM);
	int6_->preCalculate(t, M, dMdt, pars, totM);
	int7_->preCalculate(t, M, dMdt, pars, totM);
	int8_->preCalculate(t, M, dMdt, pars, totM);
	int9_->preCalculate(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	int1_->function(t, M, dMdt, pars, totM);
	int2_->function(t, M, dMdt, pars, totM);
	int3_->function(t, M, dMdt, pars, totM);
	int4_->function(t, M, dMdt, pars, totM);
	int5_->function(t, M, dMdt, pars, totM);
	int6_->function(t, M, dMdt, pars, totM);
	int7_->function(t, M, dMdt, pars, totM);
	int8_->function(t, M, dMdt, pars, totM);
	int9_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class PList>
inline void
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars)
{
	int1_->function(t, M, dMdt, pars);
	int2_->function(t, M, dMdt, pars);
	int3_->function(t, M, dMdt, pars);
	int4_->function(t, M, dMdt, pars);
	int5_->function(t, M, dMdt, pars);
	int6_->function(t, M, dMdt, pars);
	int7_->function(t, M, dMdt, pars);
	int8_->function(t, M, dMdt, pars);
	int9_->function(t, M, dMdt, pars);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM)
{
	return int1_->function(t, M, dMdt, pars, totM)+
			int2_->function(t, M, dMdt, pars, totM)+
			int3_->function(t, M, dMdt, pars, totM)+
			int4_->function(t, M, dMdt, pars, totM)+
			int5_->function(t, M, dMdt, pars, totM)+
			int6_->function(t, M, dMdt, pars, totM)+
			int7_->function(t, M, dMdt, pars, totM)+
			int8_->function(t, M, dMdt, pars, totM)+
			int9_->function(t, M, dMdt, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class PList>
inline coord<>
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
function(double t, coord<>  &M, coord<>  &dMdt, PList *pars)
{
	return int1_->function(t, M, dMdt, pars);
			int2_->function(t, M, dMdt, pars)+
			int3_->function(t, M, dMdt, pars)+
			int4_->function(t, M, dMdt, pars)+
			int5_->function(t, M, dMdt, pars)+
			int6_->function(t, M, dMdt, pars)+
			int7_->function(t, M, dMdt, pars)+
			int8_->function(t, M, dMdt, pars)+
			int9_->function(t, M, dMdt, pars);
}


template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->magneticField(t, M, pars, totM)+
	 		int2_->magneticField(t, M, pars, totM)+
	 		 int3_->magneticField(t, M, pars, totM)+
	 		 int4_->magneticField(t, M, pars, totM)+
	 		 int5_->magneticField(t, M, pars, totM)+
	 		 int6_->magneticField(t, M, pars, totM)+
	 		 int7_->magneticField(t, M, pars, totM)+
	 		 int8_->magneticField(t, M, pars, totM)+
	 		 int9_->magneticField(t, M, pars, totM);
}

template<class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
template<class PList>
inline rmatrix
  Interactions<Int_1, Int_2 ,Int_3,
               Int_4,Int_5,Int_6,
               Int_7,Int_8,Int_9>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return int1_->evolutionMatrix(t, M, pars, totM)+
	 int2_->evolutionMatrix(t, M, pars, totM)+
	  int3_->evolutionMatrix(t, M, pars, totM)+
	  int4_->evolutionMatrix(t, M, pars, totM)+
	  int5_->evolutionMatrix(t, M, pars, totM)+
	  int6_->evolutionMatrix(t, M, pars, totM)+
	  int7_->evolutionMatrix(t, M, pars, totM)+
	  int8_->evolutionMatrix(t, M, pars, totM)+
	  int9_->evolutionMatrix(t, M, pars, totM);

}

END_BL_NAMESPACE

#endif

