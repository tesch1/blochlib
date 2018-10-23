
/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-26-02
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
	multipowder.h -->
	This has several classes associated with it

	its main purpose is to allow the easy changeing between
	powder averages in the middle of performing a simulation

	however, to properly perform the powder integration using different
	powder sets, at different times, a double (or triple, etc) integral is
	needed in order to avoid powder angle correlation.  THus making the
	simulation quite long, and not only that we have an
	unknown number of loops to perform (it depends on howmany 'powder()' functions
	are called in the main sequence) so there are a few helper classes to handle
	this odd problem...they basically allow one iterator for all the loops


*/

#ifndef __multipowder_h__
#define __multipowder_h__ 1

#include "blochlib.h"
#include "pulsedata.h"
using namespace BlochLib;
using namespace std;

/*************** Data container Element for Powder Lists *****/
//this is the data element for a 'std::map' like object
//escept here we can have duplicate keys and values
template<class key_t, class val>
class nonUniqueMapEle
{
	public:
		key_t key;
		val value;

		nonUniqueMapEle(){}
		nonUniqueMapEle(const key_t &in, const val &dat):
			key(in), value(dat)
		{}
};


/**********8 The 'multi' powder iterator *******/
//this class allows for the looping of multiple powder loops
// inside a 'single' loop  the iterator returns a 'std::map<sting, int>'
// that correspond to the current posisition of each powder loop

//It is ASSUMED that the input powder std::map has unquie elements
// and the value_t has a 'size()' function

template<class key_t, class value_t>
class multiMapIter
{
	private:

		bool notended;


	public:
		typedef Vector<nonUniqueMapEle<key_t,value_t> > AuxMap;

		AuxMap *auxDat; //our data ptr

		Vector<nonUniqueMapEle<key_t,int> > Indexes;
		Vector<nonUniqueMapEle<key_t,int> > Lengths;


		multiMapIter(AuxMap &in):
			auxDat(&in)
		{
			reset();
		}

		void reset(){
			if(auxDat!=NULL){
				Indexes.resize(auxDat->size());
				Lengths.resize(auxDat->size());
				for(int i=0;i<auxDat->size();++i){
					Indexes[i].key=(*auxDat)[i].key;
					Indexes[i].value=0;

					Lengths[i].key=(*auxDat)[i].key;
					Lengths[i].value=(*auxDat)[i].value.size();
				}
			}
			notended=true;
		}

		operator bool()
		{	return notended;	}

		void operator++()
		{
			if(Indexes.size()==0){
				notended=false; return;
			}
			if(Indexes[0].value==Lengths[0].value-1){
				notended=false; return;
			}
			if(Indexes.size()==1){
				Indexes[0].value++;	return;
			}

			for(int i=1;i<Indexes.size();++i){
				if(Indexes[i].value==Lengths[i].value-1)
				{	Indexes[i-1].value++; Indexes[i].value=0;	}
				else Indexes[i].value++;

			}
		}
};


/***** THe multi Powder Class ****/

class MultiPowder
{
	private:
		typedef std::map<std::string,powder> PowMap;
		typedef std::map<std::string,powder>::iterator PowMapIter;

		PowMap Powders; //contains the powder object by name


	public:

		MultiPowder(){}

		MultiPowder(Parameters &pset);

	//this will grab all the
	//subsections named 'powder1...powderN' and add them
	// to the std::map
		void parse(Parameters &pset);

	//add a powder to the list from a parameter set list
		void addPowder(std::string name,Parameters &PowderPset);

	//add a powder to the list from a powder
		void addPowder(std::string name,const powder &PowderPset);

	//returns the pointer to the powder in the list
		powder *getPowder(std::string name);

	//returns the 'iterator' vector that are used given
	// an input PulseData
		Vector<nonUniqueMapEle<std::string,powder> > findUsedPowders(const Vector<PulseData> &pdata);

		friend std::ostream &operator<<(std::ostream &oo,const MultiPowder &out);

};


#endif


