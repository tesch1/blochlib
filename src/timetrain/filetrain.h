

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
	filetrain.h--> reads the timetrain from am ASCII file

	the file format Should be....
	<Time>  <substep>

	comments can be made in file with
	a line starting with '#'
	or the matlab type
	a line startin with '%'
	...
*/



#ifndef _filetrain_h_
#define _filetrain_h_


#include "timetrain/gentrain.h"
#include "utils/utils.h"
#include <string>
#include <iostream>


BEGIN_BL_NAMESPACE


class FileTimeEngine: public GenTimeEngine{

	private:
		std::string name_;
		const void FileErr();
		const void FileFormatErr();
		const void TimeOrderErr();
		const void SizeErr();

	public:


		FileTimeEngine():
			GenTimeEngine(1), name_("")
		{}

		FileTimeEngine(const FileTimeEngine &cp):
			GenTimeEngine(cp), name_(cp.name_)
		{}

		FileTimeEngine(const char *cp):
			GenTimeEngine(0), name_(std::string(cp))
		{}

		FileTimeEngine(std::string cp):
			GenTimeEngine(0), name_(cp)
		{}

		~FileTimeEngine(){};


/************************************
* The Master Function...the rest of everything here
* is 'fluff' To make the TimeTrain work, all you really need
* is this function....takes 3 vectors
* tO_--> The times starting at beginT---endT
* dtO_--> the time steps between each step in tO
* dstepO_--> the 'sub steps' for each dt step
*****************************/
		int TimeFunc(Vector<double> &tO_, Vector<double> &dtO_, Vector<int> &dstepO_);


/***********Functions for 'after' initialization modification
* theses are called from the Master Class "TimeTrain"
* and are thus 'overridden' by those function...typically
* after the modification the time lists must be recalulated
*************************************/
		inline void read(){}
		inline void read(std::string in){	name_=in;	}
		inline void read(const char *in){	name_=std::string(in);	}
		FileTimeEngine &operator =(const FileTimeEngine &rhs);

};


END_BL_NAMESPACE




#endif
