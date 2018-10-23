/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-11-01
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
	gentrain.h--> header for the base class of all the Time train engines
	contains the Base iterator class for the Time trains
*/

#ifndef _gentimetrain_h_
#define _gentimetrain_h_ 1

/* The base class for all subsiquent time trains... */

#include "container/Vector/Vector.h"

BEGIN_BL_NAMESPACE


class GenTimeIter;

class GenTimeEngine{
	friend class GenTimeIter;
	protected:
		int nele_;

	public:

		typedef GenTimeIter iterator;

		const void Additerr();
		const void Subiterr();
		const void Adderr();
		const void AddStepErr();
		const void SetErr();
		const void VecAssignErr();
		const void NstepErr();

		char RecordSep; //for output purposes either '\n' or ' '

		inline int size(){ return nele_;		}

		GenTimeEngine():
			 nele_(0), RecordSep('\n'){}

		GenTimeEngine(int numele):
			nele_(numele), RecordSep('\n'){}

		GenTimeEngine(const GenTimeEngine &copy):
			nele_(copy.nele_),RecordSep(copy.RecordSep){}

		GenTimeEngine &operator =(const GenTimeEngine &rhs);


		~GenTimeEngine(){}
};


class GenTimeIter{

	protected:
		GenTimeEngine *mye_;
		int curpos_; //current iter posistion
		int nextpos_;	//nest iter posistion
		bool notended_;	//have we finshed looping?
		int loops_;	//number of times to loop the list
		int curloop_; //current loop counter

	public:
		GenTimeIter():
			mye_(NULL),
			curpos_(0), nextpos_(1), notended_(true), loops_(0),curloop_(0)
		{}


		GenTimeIter(GenTimeIter &copy):
			mye_(copy.mye_),
 			curpos_(copy.curpos_), nextpos_(copy.nextpos_), notended_(copy.notended_),
			 loops_(copy.loops_),curloop_(copy.curloop_)
		{}


		GenTimeIter(GenTimeEngine &in):
			mye_(&in),
			curpos_(0), nextpos_(1), notended_(true), loops_(0),curloop_(0)
		{}

		~GenTimeIter(){ mye_=NULL;	};

		//reset the iterator
		void reset();

		GenTimeIter &operator=(const GenTimeIter &in);

		//have we gone threw the list yet?
		inline operator bool() const { return notended_; };

		//set the number of times to 'loop' around the list
		inline void setLoops(int i){ loops_=i;	}

		//standard acsess functions
		inline int loops()const { return loops_;	}
		inline int currentLoop() const { return curloop_;	}
};



END_BL_NAMESPACE



#endif
