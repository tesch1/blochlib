/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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
	rankType.h-->
here we determin the 'ranking' of various types
this is quite usefull when we multiply two template
types together...i.e. float*int==float!!  so we need
to choose the correct 'outtype'

thanks to the 'blitz++' folks for the nice and easy way
(better then writing them all out by hand...ug)

*/

#ifndef BL_rankType_h_
#define BL_rankType_h_ 1

BEGIN_BL_NAMESPACE

#define HaveComplex 1

//let it know about the complex numbers
template<class Ctype_T>
class Complex;

//#include<MasterHeader.h>  //contains our 'defines' for random bits

//***USE THIS FUNCTION TO PICK THE 'UP TYPE'*****************
#define OutType(A,B) typename testType<A,B>::Promotion

//Gives float if not double or complex, and complex or double otherwise
#define FloatType(A) typename AutoUpUp<A>::floatType
#define SumType(A) typename AutoUpUp<A>::sumType

//returns 1 if A type > B type 0 otherwise
//this should go in an IF statement as it is a booool
#define BiggerType(A,B) int(KnonwnTrait<A>::ranking)>=int(KnonwnTrait<B>::ranking)

//tells me if two number class types are equal
#define EqualType(A,B) int(KnonwnTrait<A>::ranking)==int(KnonwnTrait<B>::ranking)


/*
 * This compiler supports partial specialization, so type promotion
 * can be done the elegant way.  This implementation is after ideas
 * by Jean-Louis Leroy.
 */

//this little struct holds simply the type 'T' and the order
//we hold it in 'ranking' the bigger 'ranking' the more we bias it
//we set 'known' to 1 if we know it, 0 if not
template<class T>
struct KnonwnTrait {
    enum { ranking = 0,
           known = 0 };
};


//a macr to make life easier to set up 'known' types
#define EasyKnownTraitSetUp(T,rank)          \
    template<>                                \
    struct KnonwnTrait< T > {             \
        enum { ranking = rank,          \
           		known = 1 };           \
    };

EasyKnownTraitSetUp(int,10)
EasyKnownTraitSetUp(unsigned int,20)
EasyKnownTraitSetUp(long,30)
EasyKnownTraitSetUp(unsigned long,40)
EasyKnownTraitSetUp(float,50)
EasyKnownTraitSetUp(double,60)
//EasyKnownTraitSetUp(long double,70)
EasyKnownTraitSetUp(Complex<float>,80)
EasyKnownTraitSetUp(Complex<double>,90)

//for some types we ALWAYS want to 'move up'
//those would be the 'shorts, chars, and bool' types
//so here is yet another
template<class T>
struct AutoUp {
    typedef T typeOut;		//this is the automatic 'up' type
};

//a thing to make life a bit easier
#define EasyAutoKnown(T1,T2)     \
    template<>                            \
    struct AutoUp<T1> {        \
      typedef T2 typeOut;               \
   	};

EasyAutoKnown(bool, int)
EasyAutoKnown(char, int)
EasyAutoKnown(unsigned char, int)
EasyAutoKnown(short int, int)
EasyAutoKnown(short unsigned int, unsigned int)

//storage for promotion to T1
template<class T1, class T2, int up_to_T1>	//T1 wins
struct Up2 {
    typedef T1 Promotion;
};

//storage for promotion to T2 (note same struct as above, but set the 'up_to_T1==0'
template<class T1, class T2>
struct Up2<T1,T2,0> {
    typedef T2 Promotion;	//T2 wins
};


//the master comparer
template<class T1_orig, class T2_orig>
struct testType {
    // Handle promotion of small integers to int/unsigned int
    typedef typename AutoUp<T1_orig>::typeOut T1;
    typedef typename AutoUp<T2_orig>::typeOut T2;


    enum {
		// True if T1 is higher ranked
      T1wins =
        int(KnonwnTrait<T1>::ranking) > int(KnonwnTrait<T2>::ranking),

    // True if we know ranks for both T1 and T2
      knowBoth =
       int(KnonwnTrait<T1>::ranking) && int(KnonwnTrait<T2>::ranking),

    // True if we know T1 but not T2
      knowT1butNotT2 =  int(KnonwnTrait<T1>::ranking) && !(int(KnonwnTrait<T2>::ranking)),

    // True if we know T2 but not T1
      knowT2butNotT1 =  int(KnonwnTrait<T2>::ranking) && !(int(KnonwnTrait<T1>::ranking)),

    // True if T1 is bigger than T2
      T1big = sizeof(T1) >= sizeof(T2),

    // We know T1 but not T2: true
    // We know T2 but not T1: false
    // Otherwise, if T1 is bigger than T2: true
    //our default up type
      defUp = knowT1butNotT2 ? false :  (knowT2butNotT1 ? true : T1big)
    };

    // If we have both ranks, then use them.
    // If we have only one rank, then use the unknown type.
    // If we have neither rank, then promote to the larger type.

    enum {
      promoteToT1 = (int(knowBoth) ? int(T1wins): int(defUp)) ? 1 : 0
    };

    typedef typename Up2<T1,T2,promoteToT1>::Promotion Promotion;
};

//************FOR ALWAYS UP**************

//for some types we ALWAYS want to 'move up'
//those would be the 'shorts, chars, and bool' types
//so here is yet another
template<class T1>
struct AutoUpUp {
    typedef T1 typeOut;		//this is the automatic 'up' type
    typedef T1 floatType;
    typedef T1 sumType;
};

//a thing to make life a bit easier
#define EasyAutoKnownUp(T1,T2,T3)     \
    template<>                            \
    struct AutoUpUp< T1 > {        \
      typedef T1 typeOut;				\
      typedef T2 floatType;               \
      typedef T3 sumType;				\
   	};

EasyAutoKnownUp(bool, float, int)
EasyAutoKnownUp(char, float, int)
EasyAutoKnownUp(unsigned char, float, int)
EasyAutoKnownUp(short int, float,int)
EasyAutoKnownUp(short unsigned int, float, int)
EasyAutoKnownUp(int,float, int)
EasyAutoKnownUp(unsigned int,float, int)
EasyAutoKnownUp(long ,float, float)
EasyAutoKnownUp(unsigned long,float, float)
EasyAutoKnownUp(float,float, float)
EasyAutoKnownUp(double,double, double)
//EasyAutoKnownUp(long double,long double, long double)
EasyAutoKnownUp(Complex<double>,Complex<double>, Complex<double>)
EasyAutoKnownUp(Complex<float>,Complex<float>, Complex<float>)


//gives us the ability to specifiy a 'zero' for
//special types like Diagonal, and Identity
//it gives us the ability to return a 'const' and an '&' (address) value
//in the 'get' routines

template<class T>
class ZeroType {
	private:
		static T zz;
	public:
		typedef T numType;
		ZeroType(){}
		static T &zero(){	return zz;	}
};

template<class T>
T ZeroType<T>::zz=T(0);


template<class T>
class OneType {
	private:
		static T oo;
	public:
		typedef T numType;
		OneType(){}
		static T &one(){	return oo;	}
};

template<class T>
T OneType<T>::oo=T(1);


// 'other' trait types
// By default, the 'number' compnent
template<class NumType_t>
struct ObjectTrait {
    typedef NumType_t numtype;
};





// This macro is provided so that users can register their own
// multicomponent types.

#define MAKE_OBJECT_TRAIT(Obj,NumType_t)          \
  template<>                                                 \
  struct ObjectTrait<Obj > {                   \
    typedef NumType_t numtype;                            \
  };                                                         \


//do the common ones
MAKE_OBJECT_TRAIT(Complex<double>, Complex<double>)
MAKE_OBJECT_TRAIT(Complex<float>, Complex<float>)
MAKE_OBJECT_TRAIT(double, double)
MAKE_OBJECT_TRAIT(float, float)
MAKE_OBJECT_TRAIT(char, char)
MAKE_OBJECT_TRAIT(int, int)
MAKE_OBJECT_TRAIT(bool, bool)



END_BL_NAMESPACE


#endif
