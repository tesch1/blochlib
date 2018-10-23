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
 /* singlespinop.cc ****************************************************-*-c++-*

Generates simple ONE spin operators (given the quantum number of the system)

this is just a few functions that serve as starting points for the spin system
i do not recomend using these functions in any program....
instead use the "spin_sys" class with one spin as it stores theses output matrices
so they do not have to calculated again..

I'd like to thank to the makers of Gamma

  S.A. Smith, T.O. Levante, B.H. Meier and R.R. Ernst
  Computer Simulations in Magnetic Resonance:
  An Object-Oriented Programming Approach
  J. Magn. Reson., 106a, 75-105, (1994S, Series A)

http://gamma.magnet.fsu.edu/

for this intial set type of set up...

*/


#ifndef _singlespinop_cc
#define _singlespinop_cc 1

#include "QMspins/singlespinop.h"
#include "container/complex.h"
#include "utils/constants.h"
#include "utils/utils.h"
#include "container/matrix/matrix.h"
#include "QMspins/spinsyspars.h"
#include <map>
#include <list>

BEGIN_BL_NAMESPACE


//these functions generate matrices of size hsXhs...

//these static lists are maintatined so we only calculate the matrix ONCE
std::map<int, rimatrix> 	IeList;		//real identity matrix Ie lists
std::map<int, smatrix> 	IxList;		//real symmetric matrix Ix lists
std::map<int, hmatrix>	IyList;		//complex hermieitan matrix Iy lists
std::map<int, rdmatrix>	IzList;		//real diagonal matrix Iz list
std::map<int, rmatrix>	IpList;		//complex full matrix Iplus list
std::map<int, rmatrix>	ImiList;	//complex full matrix Iminus list




//IDENTITY MATRIX
//a real identity matrix
rimatrix sps_Ie(int hs)
{
	if(hs <1) return rimatrix();
	std::map<int,rimatrix>::iterator i;
	i = IeList.find(hs);
	if(i == IeList.end()){
		rimatrix mx(hs,hs);
		IeList.insert(std::pair<int,rimatrix>(hs, mx));
		i = IeList.find(hs);
	}
	return i->second;
}

//Ix MATRIX
//a symmetric (real) matrix
smatrix sps_Ix(int hs)
{
	if(hs <1) return smatrix();
	std::map<int,smatrix>::iterator i;
	i = IxList.find(hs);
	if(i == IxList.end())
	{
		smatrix mx(hs,hs,0);
		double q  = (hs-1.0)/2;
		double q1 = q*(q+1);
		double m, tmp;
		int k;
		for(k=0,m=q-1; k<hs-1; k++,m-=1)
		{
			tmp = sqrt(q1-m*(m+1))/2;
			mx.put(k+1,k,tmp);
		}
		IxList.insert(std::pair<int,smatrix>(hs, mx));
		return mx;
		i = IxList.find(hs);
	}
	return i->second;
}

//Ix MATRIX
//a hermitian (complex) matrix
hmatrix sps_Iy(int hs)
{
	if(hs <1) return hmatrix();
	std::map<int,hmatrix>::iterator i;
	i = IyList.find(hs);
	if(i == IyList.end())
	{
		hmatrix mx(hs,hs,0);
		double q  = (hs-1.0)/2;
		double q1 = q*(q+1);
		double m, tmp;
		int k;
		for(k=0,m=q-1; k<hs-1; k++,m-=1)
		{
			tmp = sqrt(q1-m*(m+1))/2;
			mx.put(k+1,k,complex(0,tmp));
		}
		IyList.insert(std::pair<int,hmatrix>(hs, mx));
		return mx;
		i = IyList.find(hs);
	}
	return i->second;
}

//Iz MATRIX
//a diagonal (real) matrix
rdmatrix sps_Iz(int hs)
{
	if(hs <1) return rdmatrix();
	std::map<int,rdmatrix>::iterator i;
	i = IzList.find(hs);
	if(i == IzList.end())
	{
		rdmatrix mx(hs,hs);
		double m, q=(hs-1.0)/2;
		int k;
		for(k=0,m=q; k<hs; k++,m-=1) mx.put(k,k,m);
		IzList.insert(std::pair<int,rdmatrix>(hs, mx));
		return mx;
		i = IzList.find(hs);
	}
	return i->second;
}

//Iplus MATRIX
//a full (real) matrix
rmatrix sps_Ip(int hs)
{
	if(hs <1) return rmatrix();
	std::map<int,rmatrix>::iterator i;
	i = IpList.find(hs);
	if(i == IpList.end())
	{
		rmatrix mx(hs,hs,0);
		double q  = (hs-1.0)/2;
		double q1 = q*(q+1);
		double m;
		int k;
		for(k=0,m=q-1; k<hs-1; k++,m-=1){	mx.put(k,k+1,sqrt(q1-m*(m+1))); }
		IpList.insert(std::pair<int,rmatrix>(hs, mx));
		return mx;
		i = IpList.find(hs);
	}
	return i->second;
}

//Iminus MATRIX
//a full (real) matrix
rmatrix sps_Imi(int hs)
{
	if(hs <1) return rmatrix();
	std::map<int,rmatrix>::iterator i;
	i = ImiList.find(hs);
	if(i == ImiList.end())
	{
		rmatrix mx(hs,hs,0);
		double q  = (hs-1.0)/2;
		double q1 = q*(q+1);
		double m;
		int k;
		for(k=0,m=q-1; k<hs-1; k++,m-=1)	mx.put(k+1,k,sqrt(q1-m*(m+1)));
		ImiList.insert(std::pair<int,rmatrix>(hs, mx));
		return mx;
		i = ImiList.find(hs);
	}
	return i->second;
}


rmatrix Imi(int hs){ return sps_Imi(hs);	}
rmatrix Ip(int hs){ return sps_Ip(hs);	}
rdmatrix Iz(int hs){ return sps_Iz(hs);	}
hmatrix Iy(int hs){ return sps_Iy(hs);	}
smatrix Ix(int hs){ return sps_Ix(hs);	}
rimatrix Ie(int hs){ return sps_Ie(hs);	}

//Generates all the entries for that 'hs' dimension
void GenSingleOps(int hs)
{
	sps_Ie(hs);
	sps_Ix(hs);
	sps_Iy(hs);
	sps_Iz(hs);
	sps_Ip(hs);
	sps_Imi(hs);
}

END_BL_NAMESPACE



#endif

