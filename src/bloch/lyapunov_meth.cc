
/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08-06-01
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
	lyapunov_meth.cc-->a class to calculate Lyapunov exponents from
	 a chunk of evolved Variational Equation data...the COMPILABLE methods

*/

#ifndef _lyapunov_meth_cc_
#define _lyapunov_meth_cc_ 1

#include "bloch/lyapunov.h"

BEGIN_BL_NAMESPACE


/*
  the lyapunov exponent is a measure of the deviation from the current trajectory
  if we move by small amount 'delat x'  if this number is positive the tragectory
  deviates very quickly, if it is negative the path converges...The direction of converges or
  deviation is 'ficticious' the basis is typically NOT the same as the differentical equation basis
  the reason for this is simple: at each point in the tragectory, the Variational equations form a
  tensorial basis for the entire field of the diff eqs.  the eigen vectors of this matrix equation
  form it basis, and becuase the Varitational Eqs CHANGE as the trajectory changes, this basis
  it typically NOT constant.

  The problem with calculation of the lyapunov exponents is a question of numerical stability.
  if a path diverges quickly, the the eigenvalue of that direction explodes to infinity
  and more then likely there will be a negative lyapunov exponent, this eigenvalue quickly goes to zero.
  thus we have typically an underflow error, AND an overflow error..not to mention calcualation of
  the eigen values becomes quite expenive..

  so we choose a renormalization metod.  concider a single vector in our Variational basis
  if along this direction the growth is positive, at each step we can get a value for the
  growth by takeing the length of the vector (we are assuming that the length was intially 1)
  then renormalize the vector (as well as making it orthogonal to the rest of the basis, to
  span the variational space as well possible).  The next iteration of the equations shold then
  expand our new basis by the same amount...

  so each lyapunov at each step is calculated by orthogonalize each variation vector, using
  the norm as the measure, and replacing the vectors back in the solver as new intial conditions
*/

Lyapunov<coord<> > &Lyapunov<coord<> >::operator=(Lyapunov<coord<> > &rhs)
{
	if(this==&rhs) return *this;
	data_=rhs.data_;
	curct_=rhs.curct_;
	lyps=rhs.lyps;
	ts=rhs.ts;
	dts=rhs.dts;
	calcstep_=rhs.calcstep_;
	return *this;
}


void Lyapunov<coord<> >::calcLyapunov(double ontime,double dt)
{
	static rmatrix Jt(3,3,0);
	static Vector<double> norms(3,0);
	if(!data_ || !size_ || data_->size() != (size_*4))
	{
		std::cerr<<std::endl<<"Error:: Lyapunov.calcLyapunov "<<std::endl;
		std::cerr<<" Cannot calculate any Lyapunov's as there has been"<<std::endl;
		std::cerr<<" no data ptr set, or the size of the system is zero..."<<std::endl;
		return;
	}
	if(curct_%calcstep_==0)
	{
		ts.push_back(ontime);
		dts.push_back(dt);
		int j=0;
		for(int i=size_;i<size_*4;i+=3, ++j)
		{
			//cout<<j<<" "<<i<<" "<<size_<<endl;
			Jt.putCol( 0, data_->get(i));
			Jt.putCol( 1, data_->get(i+1));
			Jt.putCol( 2, data_->get(i+2));
			//cout<<Jt<<endl;
			Jt=GramSchmidt(Jt, norms);
			//cout<<norms<<" lyps: "<<lyps[j]<<" Logn:  "<<log(norms)<<endl;
			lyps[j]+=log(norms);
			data_->put(i,Jt.col(0));
			data_->put(i+1,Jt.col(1));
			data_->put(i+2,Jt.col(2));
		}
		++pos_;
	}

	++curct_;
}

void Lyapunov<coord<> >::reset()
{
	lyps.fill(0);
	ts.resize(0);
	dts.resize(0);
	curct_=0;
	pos_=0;
	//reset the initial condition bits to the identity matrix!
	if(data_ && size_ && data_->size() == (size_*4))
	{
		for(int i=size_;i<size_*4;i+=3)
		{
			data_->put(i,coord<>(1,0,0));
			data_->put(i+1,coord<>(0,1,0));
			data_->put(i+2,coord<>(0,0,1));
		}
	}
}


/*
writes the data to a file....
in the following format

	time lyp1 lyp2 lyp3 ... lypN


*/


void Lyapunov<coord<> >::write(const char *fname)
{
	std::fstream ouf(fname, std::ios::binary| std::ios::out);
	if(ouf.fail())
	{
		std::cerr<<std::endl<<"Error Lyapunov.write() "<<std::endl;
		std::cerr<<" cannot open file for writing...."<<std::endl;
		return;
	}
	write(ouf);
}


void Lyapunov<coord<> >::write(std::fstream &ouf)
{
	if(ouf.fail())
	{
		std::cerr<<std::endl<<"Error Lyapunov.write() "<<std::endl;
		std::cerr<<" cannot open file for writing...."<<std::endl;
		return;
	}
	if(continprint_)
	{
		if((curct_-1)%calcstep_==0 && (curct_-1)!=0) //the first calcd point is always botched
		{
			ouf.write((char *)&ts(pos_-1), sizeof(double));
			static coord<> tmm;
			for(int i=0;i<lyps.size();++i)
			{
				tmm=lyps[i]/dts(pos_-1)/(curct_-1);
				ouf.write((char *)&tmm.x(), sizeof(double));
				ouf.write((char *)&tmm.y(), sizeof(double));
				ouf.write((char *)&tmm.z(), sizeof(double));
			}
			ouf<<std::endl;
		}
	}else{}

}



void Lyapunov<coord<> >::print(const char *fname)
{
	std::ofstream ouf(fname);
	if(ouf.fail())
	{
		std::cerr<<std::endl<<"Error Lyapunov.print() "<<std::endl;
		std::cerr<<" cannot open file for writing...."<<std::endl;
		return;
	}
	print(ouf);
}


void Lyapunov<coord<> >::print(std::ostream &ouf)
{
	if(ouf.fail())
	{
		std::cerr<<std::endl<<"Error Lyapunov.write() "<<std::endl;
		std::cerr<<" cannot open file for writing...."<<std::endl;
		return;
	}
	if((curct_-1)%calcstep_==0 && (curct_-1)!=0) //the first calcd point is always botched
	{
		ouf<<ts(pos_-1)<<" ";
		for(int i=0;i<lyps.size();++i)
		{
			ouf<<lyps[i]/dts(pos_-1)/(pos_-1)<<" ";
			//cout<<lyps[i]<< " dt: "<<dts(pos_-1)<<" pos: "<<(pos_-1)<<endl;

		}
		ouf<<std::endl;
	}
}


std::ostream &operator<<(std::ostream &oo, Lyapunov<coord<> > &rhs)
{
	rhs.print(oo);
	return oo;
}

std::fstream &operator<<(std::fstream &oo, Lyapunov<coord<> > &rhs)
{
	rhs.write(oo);
	return oo;
}

END_BL_NAMESPACE



#endif


