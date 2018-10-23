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

/*
 integrate.h --> integrates a function...using an adaptive simpson method...
 algorithm adapted from W.Gander and W.Gautschi "Adaptive Quarature revisitied" , 1998.
 (atapted from Matlab 6.0 'quad.m')
*/

#ifndef _integrate_h_
#define _integrate_h_ 1


#include "container/complex.h"
#ifndef ON_WINDOWS
	#include "blochconfig.h"
#endif

BEGIN_BL_NAMESPACE


template<class Func_T, class T=double>
class Integrate
{
	private:
		Func_T *mf;

		const void NoMf()
		{
			std::cerr<<std::endl<<"Error: Integrate(...)"<<std::endl;
			std::cerr<<" NO Function defined...(i.e. a class with "<<std::endl;
			std::cerr<<" T=Func_T::operator(x) "<<std::endl;
		}

		const void InfFunc()
		{
			std::cerr<<std::endl<<"Error: Integrate(...)"<<std::endl;
			std::cerr<<" Function/Integral evaluation at points is infinite"<<std::endl;
			std::cerr<<" cannot calulate integral"<<std::endl;
		}

		const void SingErr()
		{
			std::cerr<<std::endl<<"Error: Integrate(...)"<<std::endl;
			std::cerr<<" Sigularity Possible...step size became TOO small"<<std::endl;
			std::cerr<<" cannot finish calculating integral"<<std::endl;
		}

		const void SingErr2()
		{
			std::cerr<<std::endl<<"Error: Integrate(...)"<<std::endl;
			std::cerr<<" Sigularity Possible...went over max number of iterations"<<std::endl;
			std::cerr<<" cannot finish calculating integral"<<std::endl;
		}

	public:
		T begin;
		T end;
		T eps;
		T step;
		T minst;
		int maxct;
		int ct;

		Integrate():
			begin(0), end(1), eps(1.e-6), step(0.1), maxct(10000), ct(0)
		{
			mf=new Func_T;
		}

		Integrate(Func_T &in, T be=0, T en=1, T epss=1.e-6):
			begin(be), end(en), eps(epss ), step((be-en)/10), maxct(10000),ct(0)
		{
			mf=new Func_T;
			mf=&in;
		}

		~Integrate(){	mf=NULL;	delete mf;	}

		int FunctionCount()const{	return ct;	}
		int count() const {	return ct;	}

		T StartQuad(T a, T b)
		{
			begin=a;
			end=b;
			step=(a+b)/2.0;
			minst=eps/1024.0*abs(b-a);
			double fa, fb, fc;
			fa=mf->operator()(a);
			fc=mf->operator()(step);
			fb=mf->operator()(b);
			ct+=3;

	#ifdef HAVE_ISNAN
			if(hasnan(fa)){	fa=mf->operator()(begin+eps*(end-begin)); ct++;	}
			if(hasnan(fb)){	fb=mf->operator()(end-eps*(end-begin));	ct++;	}
	#endif
			return Quad(a,b, fa, fc, fb);
		}

		T Quad(T a, T b, T fa, T fc, T fb)
		{
			T Q, h, c ;
			h=b-a;
			c=(a+b)/2;

			if(abs(h) < minst || c==begin || c==end)
			{
				Q=fc*h;
				SingErr();
				return Q;
			}

			T y1, y2, fe, fd;
			y1=mf->operator()((a+c)/2.0);
			y2=mf->operator()((c+b)/2.0);
			ct+=2;
			if(ct>maxct)
			{
				Q=h*fc;
				SingErr2();
				return Q;
			}

			fd=y1;
			fe=y2;

			T Q1, Q2;
			Q1=(h/6.0)*(fa+4.0*fc+fb);
			Q2=(h/12.0)*(fa+4.0*fd+2.0*fc+4.0*fe+fb);
			Q=Q2+(Q2-Q1)/15.0;

			//cout<<ct<<" "<<a<<" "<<h<<" "<<Q<<endl;

			if(hasnan(Q))
			{
				InfFunc();
				return Q;
			}

			if(abs(Q2-Q) <= eps)	return Q;
			T Qa, Qb;
			Qa=Quad(a,c,fa,fd,fc);
			Qb=Quad(c,b,fc,fe,fb);
			Q=Qa+Qb;
			return Q;
		}

		T operator()()
		{	ct=0; return StartQuad(begin, end);	}

		T operator()(T t1, T t2)
		{	ct=0;return StartQuad(t1, t2);	}

		T integrate(Func_T &in, T t1, T t2)
		{
			mf=&in;
			ct=0;
			return StartQuad(t1, t2);
		}

		T integrate(T t1, T t2)
		{
			ct=0;
			return StartQuad(t1, t2);
		}

		template<class T2>
		Vector<T> integrate(const Vector<T2> &time)
		{
			if(time.size()<=1){
				std::cerr<<std::endl<<"Error: Integrate::integrate(vector)"<<std::endl;
				std::cerr<<" input time vector is too small (size<=1)"<<std::endl;
				std::cerr<<" returning empty vector"<<std::endl;
				return Vector<T>();
			}
			Vector<T> out(time.size()-1, 0);

			out[0]=integrate(T(time[0]), T(time[1]));
			for(int i=2;i<time.size();++i){
				out[i-1]=integrate(T(time[i-1]), T(time[i]))+out[i-2];
			}
			return out; //-sum(out)/out.size();
		}

		template<class T2>
		Vector<T> integrate(const Spread<T2> &time)
		{	return integrate(Vector<T2>(time));	}

		Vector<T> integrate(const Range &time)
		{	return integrate(Vector<T>(time));	}

};

END_BL_NAMESPACE


#endif

