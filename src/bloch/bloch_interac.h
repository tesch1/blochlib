


/* bloch_interac.h ********/


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
 	bloch_interac.h-->axilliary interactions for the Bloch equations
 */


#ifndef _bloch_interatction_h_
#define _bloch_interatction_h_ 1

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
//  class DipoleDipole
//    --> Dipole-Dipole interaction (not averaged over the sample like DemagField)


/* here are the remaing interactions that are passed as a template parameter to the
	'Bloch' class....There is also a class named "Interactions" that simply act as a
	transfer class for each interactions type
*/


// There are several enums that tell the laster classes WHAT sort of input the equations need
// many of the interactions need the total magnitization (not just the one of its spin)


//this is the basic interactions (handeled by the 'bloch' itself so it is just a place holder)
class BasicInter{};

//this is the 'null' container lets 'Bloch' handle the interactions
//itself (normal ones like offset, T2, T1)

template<	class Int_1=BasicInter, class Int_2=BasicInter, class Int_3=BasicInter,
			class Int_4=BasicInter, class Int_5=BasicInter, class Int_6=BasicInter,
			class Int_7=BasicInter, class Int_8=BasicInter, class Int_9=BasicInter>
class Interactions;

typedef Interactions<BasicInter,BasicInter,BasicInter,BasicInter,BasicInter,BasicInter> NoInteractions;

//a single 'other' interaction
template<	class Int_1>
class Interactions<Int_1, BasicInter,BasicInter,
                    BasicInter,BasicInter,BasicInter,
                    BasicInter,BasicInter,BasicInter>
{
	private:
		Int_1 *int1_;
	public:

		static const int TotalMag=Int_1::TotalMag;
		static const int TotalVector=Int_1::TotalVector;
		static const int PreCalc=Int_1::PreCalc;

		Interactions():
			Int_1(0) {}
		Interactions(Int_1 &in)
		{
			int1_=&in;
		}

		~Interactions(){	int1_=NULL; 	}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);


		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};

//a double 'other' interaction
template<	class Int_1, class Int_2>
class Interactions<Int_1, Int_2 ,BasicInter,
                   BasicInter,BasicInter,BasicInter,
                   BasicInter,BasicInter,BasicInter>
{
	private:
		Int_1 *int1_;
		Int_2 *int2_;
	public:

		static const int TotalMag=Int_1::TotalMag + Int_2::TotalMag;
		static const int TotalVector=Int_1::TotalVector + Int_2::TotalVector;
		static const int PreCalc=Int_1::PreCalc + Int_2::PreCalc;


		Interactions():
			int1_(0), int2_(0) {}
		Interactions(Int_1 &in1, Int_2 &in2)
		{
			int1_=&in1;
			int2_=&in2;
		}

		~Interactions(){	int1_=NULL; int2_=NULL; 	}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};

//a trippple 'other' interaction
template<	class Int_1, class Int_2, class Int_3>
class Interactions<Int_1, Int_2 ,Int_3,
                   BasicInter,BasicInter,BasicInter,
                   BasicInter,BasicInter,BasicInter>
{
	private:
		Int_1 *int1_;
		Int_2 *int2_;
		Int_3 *int3_;
	public:
		static const int TotalMag=Int_1::TotalMag + Int_2::TotalMag + Int_3::TotalMag;
		static const int TotalVector=Int_1::TotalVector + Int_2::TotalVector+ Int_3::TotalVector;
		static const int PreCalc=Int_1::PreCalc + Int_2::PreCalc+ Int_3::PreCalc;


		Interactions():
			int1_(0), int2_(0), int3_(0) {}
		Interactions(Int_1 &in1, Int_2 &in2, Int_3 &in3)
		{
			int1_=&in1;
			int2_=&in2;
			int3_=&in3;
		}

		~Interactions(){	int1_=NULL; int2_=NULL; int3_=NULL;	}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);


		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);
/*
		template<class PList>
		inline void function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);
*/
		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};


//a quad 'other' interaction
template<	class Int_1, class Int_2, class Int_3, class Int_4>
class Interactions<Int_1, Int_2 ,Int_3,
                   Int_4,BasicInter,BasicInter,
                    BasicInter,BasicInter,BasicInter>
 {
	private:
		Int_1 *int1_;
		Int_2 *int2_;
		Int_3 *int3_;
		Int_4 *int4_;
	public:
		static const int TotalMag=Int_1::TotalMag + Int_2::TotalMag + Int_3::TotalMag + Int_4::TotalMag;
		static const int TotalVector=Int_1::TotalVector + Int_2::TotalVector+ Int_3::TotalVector+ Int_4::TotalVector;
		static const int PreCalc=Int_1::PreCalc + Int_2::PreCalc+ Int_3::PreCalc+ Int_4::PreCalc;


		Interactions():
			int1_(0), int2_(0), int3_(0), int4_(0)
		{}
		Interactions(Int_1 &in1, Int_2 &in2, Int_3 &in3, Int_4 &in4)
		{
			int1_=&in1;
			int2_=&in2;
			int3_=&in3;
			int4_=&in4;
		}

		~Interactions(){	int1_=NULL; int2_=NULL; int3_=NULL; int4_=NULL;	}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};


//a quint 'other' interaction
template<	class Int_1, class Int_2, class Int_3, class Int_4, class Int_5>
class Interactions<Int_1, Int_2 ,Int_3,
                   Int_4,Int_5,BasicInter,
                    BasicInter,BasicInter,BasicInter>
{
	private:
		Int_1 *int1_;
		Int_2 *int2_;
		Int_3 *int3_;
		Int_4 *int4_;
		Int_5 *int5_;

	public:
		static const int TotalMag=Int_1::TotalMag +
		   	                      Int_2::TotalMag +
		   	                      Int_3::TotalMag +
		   	                      Int_4::TotalMag +
		   	                      Int_5::TotalMag;
		static const int TotalVector=Int_1::TotalVector +
		                             Int_2::TotalVector+
		                             Int_3::TotalVector+
		                             Int_4::TotalVector+
		                             Int_5::TotalVector;
		static const int PreCalc=Int_1::PreCalc +
		                         Int_2::PreCalc+
		                         Int_3::PreCalc+
		                         Int_4::PreCalc+
		                         Int_5::PreCalc;


		Interactions():
			int1_(0), int2_(0), int3_(0), int4_(0), int5_(0)
		{}
		Interactions(Int_1 &in1, Int_2 &in2, Int_3 &in3, Int_4 &in4, Int_5 &in5)
		{
			int1_=&in1;
			int2_=&in2;
			int3_=&in3;
			int4_=&in4;
			int5_=&in5;
		}

		~Interactions(){	int1_=NULL; int2_=NULL; int3_=NULL; int4_=NULL; int5_=NULL;	}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};

//a sextet 'other' interaction
template<	class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6>
class Interactions<Int_1, Int_2 ,Int_3,
                   Int_4,Int_5,Int_6,
                    BasicInter,BasicInter,BasicInter>
{
	private:
		Int_1 *int1_;
		Int_2 *int2_;
		Int_3 *int3_;
		Int_4 *int4_;
		Int_5 *int5_;
		Int_6 *int6_;

	public:
		static const int TotalMag=Int_1::TotalMag +
		   	                      Int_2::TotalMag +
		   	                      Int_3::TotalMag +
		   	                      Int_4::TotalMag +
		   	                      Int_5::TotalMag +
		   	                      Int_6::TotalMag;
		static const int TotalVector=Int_1::TotalVector +
		                             Int_2::TotalVector+
		                             Int_3::TotalVector+
		                             Int_4::TotalVector+
		                             Int_5::TotalVector+
		                             Int_6::TotalVector;
		static const int PreCalc=Int_1::PreCalc +
		                         Int_2::PreCalc+
		                         Int_3::PreCalc+
		                         Int_4::PreCalc+
		                         Int_5::PreCalc+
		                         Int_6::PreCalc;

		Interactions():
			int1_(0), int2_(0), int3_(0),
			int4_(0), int5_(0), int6_(0)
		{}

		Interactions(Int_1 &in1, Int_2 &in2, Int_3 &in3,
		             Int_4 &in4, Int_5 &in5, Int_6 &in6)
		{
			int1_=&in1;
			int2_=&in2;
			int3_=&in3;
			int4_=&in4;
			int5_=&in5;
			int6_=&in6;
		}

		~Interactions()
		{	int1_=NULL; int2_=NULL; int3_=NULL; int4_=NULL; int5_=NULL;	int6_=NULL;	}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};

//a spetet 'other' interaction
template<	class Int_1, class Int_2, class Int_3,
class Int_4, class Int_5, class Int_6,
class Int_7>
class Interactions<Int_1, Int_2 ,Int_3,
                   Int_4,Int_5,Int_6,
                    Int_7,BasicInter,BasicInter>
{
	private:
		Int_1 *int1_;
		Int_2 *int2_;
		Int_3 *int3_;
		Int_4 *int4_;
		Int_5 *int5_;
		Int_6 *int6_;
		Int_7 *int7_;

	public:
		static const int TotalMag=Int_1::TotalMag +
		   	                      Int_2::TotalMag +
		   	                      Int_3::TotalMag +
		   	                      Int_4::TotalMag +
		   	                      Int_5::TotalMag +
		   	                      Int_6::TotalMag +
		   	                      Int_7::TotalMag;
		static const int TotalVector=Int_1::TotalVector +
		                             Int_2::TotalVector+
		                             Int_3::TotalVector+
		                             Int_4::TotalVector+
		                             Int_5::TotalVector+
		                             Int_6::TotalVector+
		                             Int_7::TotalVector;
		static const int PreCalc=Int_1::PreCalc +
		                         Int_2::PreCalc+
		                         Int_3::PreCalc+
		                         Int_4::PreCalc+
		                         Int_5::PreCalc+
		                         Int_6::PreCalc+
		                         Int_7::PreCalc;

		Interactions():
			int1_(0), int2_(0), int3_(0), int4_(0), int5_(0), int6_(0), int7_(0)
		{}

		Interactions(Int_1 &in1, Int_2 &in2, Int_3 &in3,
					Int_4 &in4, Int_5 &in5, Int_6 &in6,
					Int_7 &in7)
		{
			int1_=&in1;
			int2_=&in2;
			int3_=&in3;
			int4_=&in4;
			int5_=&in5;
			int6_=&in6;
			int7_=&in7;
		}

		~Interactions()
		{	int1_=NULL; int2_=NULL; int3_=NULL;
		    int4_=NULL; int5_=NULL;	int6_=NULL;
		    int7_=NULL;
		}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};

//a octet 'other' interaction
template<	class Int_1, class Int_2, class Int_3,
class Int_4, class Int_5, class Int_6,
class Int_7, class Int_8>
class Interactions<Int_1, Int_2 ,Int_3,
                   Int_4,Int_5,Int_6,
                    Int_7,Int_8,BasicInter>
{
	private:
		Int_1 *int1_;
		Int_2 *int2_;
		Int_3 *int3_;
		Int_4 *int4_;
		Int_5 *int5_;
		Int_6 *int6_;
		Int_7 *int7_;
		Int_8 *int8_;

	public:
		static const int TotalMag=Int_1::TotalMag +
		   	                      Int_2::TotalMag +
		   	                      Int_3::TotalMag +
		   	                      Int_4::TotalMag +
		   	                      Int_5::TotalMag +
		   	                      Int_6::TotalMag +
		   	                      Int_7::TotalMag +
		   	                      Int_8::TotalMag;
		static const int TotalVector=Int_1::TotalVector +
		                             Int_2::TotalVector+
		                             Int_3::TotalVector+
		                             Int_4::TotalVector+
		                             Int_5::TotalVector+
		                             Int_6::TotalVector+
		                             Int_7::TotalVector+
		                             Int_8::TotalVector;
		static const int PreCalc=Int_1::PreCalc +
		                         Int_2::PreCalc+
		                         Int_3::PreCalc+
		                         Int_4::PreCalc+
		                         Int_5::PreCalc+
		                         Int_6::PreCalc+
		                         Int_7::PreCalc+
		                         Int_8::PreCalc;

		Interactions():
			int1_(0), int2_(0), int3_(0),
			int4_(0), int5_(0), int6_(0),
			int7_(0), int8_(0)
		{}

		Interactions(Int_1 &in1, Int_2 &in2, Int_3 &in3,
					Int_4 &in4, Int_5 &in5, Int_6 &in6,
					Int_7 &in7, Int_8 &in8)
		{
			int1_=&in1;
			int2_=&in2;
			int3_=&in3;
			int4_=&in4;
			int5_=&in5;
			int6_=&in6;
			int7_=&in7;
			int8_=&in8;
		}

		~Interactions()
		{	int1_=NULL; int2_=NULL; int3_=NULL;
		    int4_=NULL; int5_=NULL;	int6_=NULL;
		    int7_=NULL; int8_=NULL;
		}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};

//a octet 'other' interaction
template<	class Int_1, class Int_2, class Int_3,
	class Int_4, class Int_5, class Int_6,
	class Int_7, class Int_8, class Int_9>
class Interactions
{
	private:
		Int_1 *int1_;
		Int_2 *int2_;
		Int_3 *int3_;
		Int_4 *int4_;
		Int_5 *int5_;
		Int_6 *int6_;
		Int_7 *int7_;
		Int_8 *int8_;
		Int_9 *int9_;

	public:
		static const int TotalMag=Int_1::TotalMag +
		   	                      Int_2::TotalMag +
		   	                      Int_3::TotalMag +
		   	                      Int_4::TotalMag +
		   	                      Int_5::TotalMag +
		   	                      Int_6::TotalMag +
		   	                      Int_7::TotalMag +
		   	                      Int_8::TotalMag +
		   	                      Int_9::TotalMag;
		static const int TotalVector=Int_1::TotalVector +
		                             Int_2::TotalVector+
		                             Int_3::TotalVector+
		                             Int_4::TotalVector+
		                             Int_5::TotalVector+
		                             Int_6::TotalVector+
		                             Int_7::TotalVector+
		                             Int_8::TotalVector+
		                             Int_9::TotalVector;
		static const int PreCalc=Int_1::PreCalc +
		                         Int_2::PreCalc+
		                         Int_3::PreCalc+
		                         Int_4::PreCalc+
		                         Int_5::PreCalc+
		                         Int_6::PreCalc+
		                         Int_7::PreCalc+
		                         Int_8::PreCalc+
		                         Int_9::PreCalc;

		Interactions():
			int1_(0), int2_(0), int3_(0),
			int4_(0), int5_(0), int6_(0),
			int7_(0), int8_(0), int9_(0)
		{}

		Interactions(Int_1 &in1, Int_2 &in2, Int_3 &in3,
					Int_4 &in4, Int_5 &in5, Int_6 &in6,
					Int_7 &in7, Int_8 &in8, Int_9 &in9)
		{
			int1_=&in1;
			int2_=&in2;
			int3_=&in3;
			int4_=&in4;
			int5_=&in5;
			int6_=&in6;
			int7_=&in7;
			int8_=&in8;
			int9_=&in9;
		}

		~Interactions()
		{	int1_=NULL; int2_=NULL; int3_=NULL;
		    int4_=NULL; int5_=NULL;	int6_=NULL;
		    int7_=NULL; int8_=NULL; int9_=NULL;
		}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, PList *pars);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};

END_BL_NAMESPACE

#include "bloch_interac_meth.h"

#endif


