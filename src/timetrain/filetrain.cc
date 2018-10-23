

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
	<beginStepTime> <substep>
	...
*/



#ifndef _filetrain_cc_
#define _filetrain_cc_


#include "timetrain/gentrain.h"
#include "utils/utils.h"
#include "timetrain/filetrain.h"
#include <string>
#include <iostream>


BEGIN_BL_NAMESPACE


const void FileTimeEngine::FileErr()
{
	std::cerr<<std::endl<<"Error:: FileTimeTrain"<<std::endl;
	std::cerr<<" Either the file cannot be read OR you have given"<<std::endl;
	std::cerr<<" me no file name to read...."<<std::endl;
}


const void FileTimeEngine::FileFormatErr()
{
	std::cerr<<std::endl<<"Error:: FileTimeTrain"<<std::endl;
	std::cerr<<" The file format should contain only 2 columns"<<std::endl;
	std::cerr<<" <Time> <subSteps>...."<<std::endl;
	std::cerr<<" Ignoring the line...."<<std::endl;
}

const void FileTimeEngine::TimeOrderErr()
{
	std::cerr<<std::endl<<"Error:: FileTimeTrain"<<std::endl;
	std::cerr<<" the t[i-1] is bigger then the t[i]"<<std::endl;
	std::cerr<<" Time goes foward silly...although i'm sure you wish it did not"<<std::endl;
}


const void FileTimeEngine::SizeErr()
{
	std::cerr<<std::endl<<"Error:: FileTimeTrain"<<std::endl;
	std::cerr<<" The file contained no readable time slots"<<std::endl;
	std::cerr<<" it is as if you paused the world for sec..."<<std::endl;
}

/************************************
* The Master Function...the rest of everything here
* is 'fluff' To make the TimeTrain work, all you really need
* is this function....takes 3 vectors
* tO_--> The times starting at beginT---endT
* dtO_--> the time steps between each step in tO
* dstepO_--> the 'sub steps' for each dt step
*****************************/
int FileTimeEngine::TimeFunc(Vector<double> &tO_, Vector<double> &dtO_, Vector<int> &dstepO_)
{

	if(name_=="")
	{
		FileErr();
		return 0;
	}
	std::ifstream infile(name_.c_str());

	if(infile.fail())
	{
		FileErr();
		return 0;
	}
	tO_.resize(0);
	dtO_.resize(0);
	dstepO_.resize(0);
	Vector<std::string> tmp;
	char liner[1000];
	int temp=0, gotone=0;
	while((temp=infile.peek())!=EOF)
	{
		infile.getline(liner,1000,'\n');
		tmp=parse_param(liner);
		if(!tmp.empty() && tmp[0][0]!='#' && tmp[0][0]!='%')
		{
			if(tmp.size()==2)
			{
				tO_.push_back(atof(tmp[0].c_str()));
				if(gotone>0)
				{
					if(tO_[gotone-1]>=tO_[gotone])
					{
						TimeOrderErr();
						return 0;
					}
					dtO_.push_back(tO_[gotone]-tO_[gotone-1]);
					dstepO_.push_back(atoi(tmp[1].c_str()));
				}
			}
			gotone++;
		}
	}
	if(tO_.size()==0)
	{
		SizeErr();
		return 0;
	}
	nele_=tO_.size();
	return 1;
}




/***********Functions for 'after' initialization modification
* theses are called from the Master Class "TimeTrain"
* and are thus 'overridden' by those function...typically
* after the modification the time lists must be recalulated
*************************************/
FileTimeEngine &FileTimeEngine::operator =(const FileTimeEngine &rhs)
{
	if(this==&rhs) return *this;
	name_=rhs.name_;
	return *this;
}




END_BL_NAMESPACE



#endif
