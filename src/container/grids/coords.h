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
	coord.h--> a coordinate class
*/

/*simple little x,y,z class that stores a point in space*/

#ifndef _coord_h_
#define _coord_h_ 1

#include "container/Vector/Vector.h"
#include "container/operations.h"
#include "container/matrix/matrix.h"
#include "utils/blassert.h"
#include "math.h"

BEGIN_BL_NAMESPACE


#ifndef PI2
#define PI2 6.28318530717958647692
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif


/*here N is the cordiant depth (default is '3')
for N=2, you have coords labled x,y
for N=3, you have coords labled x,y,z
for N=4, yous have coords labled x,y,z,t

NOTE::size checking here is left OUT accept at initiallization
The base class "Vector" can have Bounds checking turned on BUT
there is a significant speed loss....


*/
// coord trait
template<class NumType_t, int N>
struct ObjectTrait<coord<NumType_t,N> > {
    typedef NumType_t numtype;
};

//foward decs
template<class T, int N>
class coordIter;

template<class T, int N>
class coordIterConst;

enum coordtype{cartesian, spherical,  cylindrical};

template<class T=double, int N=3>
class coord {
	private:
		//coordtype TheType;
		T data_[N];
	public:
		typedef T numtype;
		typedef coordIter<T,N> iterator;
		typedef coordIterConst<T,N> const_iterator;

		coord(){} //: TheType(cartesian) {};

		template<class T1>
		inline coord(const coord<T1, N> &in);

		template<class T1>
		inline coord(const Vector<T1> &in);

		template<class expr>
		inline coord(const VExpr<expr> &in);

		inline coord(T *data);

		template<class T1>
		inline coord(T1 iso)//, coordtype intype=cartesian):TheType(intype)
		{
			for(int i=0;i<N;++i) data_[i]=T(iso);
		}


		template<class T1>
		inline coord(T1 x, T1 y)  //, coordtype intype=cartesian): TheType(intype)
		{
			if(N<2){
				BLEXCEPTION(std::string(" To many values for this initializer")+
						std::string("\n  The coord depth is set to ")+itost(N))
			}
			data_[0]=T(x);
			data_[1]=T(y);
		}

		template<class T1>
		inline coord(T1 x, T1 y, T1 z)  //, coordtype intype=cartesian): TheType(intype){
		{
			if(N<3){
				BLEXCEPTION(std::string(" To many values for this initializer")+
						std::string("\n  The coord depth is set to ")+itost(N))
			}
			data_[0]=T(x);
			data_[1]=T(y);
			data_[2]=T(z);
		}

		template<class T1>
		inline coord(T1 x, T1 y, T1 z, T1 t)  //, coordtype intype=cartesian): TheType(intype){
		{
			if(N<4){
				BLEXCEPTION(std::string(" To many values for this initializer")+
						std::string("\n  The coord depth is set to ")+itost(N))
			}
			data_[0]=T(x);
			data_[1]=T(y);
			data_[2]=T(z);
			data_[3]=T(t);

		}

		template<class Expr, class Op>
		inline void assign(Expr expr, Op up);

		template<class T1, int N1>
		inline coord<T1,N1> alterSize(int N1)
		{
			coord<T1, N1> tm;//(TheType);
			for(int i=0;i<N1;i++){
				tm(i)=T1((*this)(i));
			}
			return tm;
		}

		inline Vector<T> Vect(){ return Vector<T>((*this).data_, N);	}

		//inline coordtype type(){ return TheType;	}

		inline T operator()(int i)const { return data_[i];	}
		inline T & operator()(int i)  { return data_[i];	}

		inline T operator[](int i)const { return data_[i];	}
		inline T & operator[](int i) { return data_[i];	}

		template<class T2>
		void fill(const T2 &in){		*this=in;		}

		template<class T1>
		inline void put(int i, const T1 &in){ data_[i]=T(in);	}

		inline void put(Range r, const T &in)
		{
			for(int i=r.first(0);i<r.last(N);++i)
			{
				data_[i]=in;
			}
		}

		template<class T2,int N2>
		inline void put(Range r, const coord<T2, N2> &in)
		{
			for(int i=r.first(0);i<r.last(N2);++i)
			{
				data_[i]=T(in[i]);
			}
		}

		template<class T2, int N2>
		inline void put(const coord<T2, N2> &in)
		{
			CompTimeAssert(N>=N2);
			for(int i=0;i<N2;++i)
			{
				data_[i]=t(in[i]);
			}
		}

		inline T  &get(int i) { return data_[i];	}
		inline T  get(int i) const { return data_[i];	}

		inline T *data(){ 	return data_;	}
		inline const T * data() const {	return data_;	}

		inline const_iterator const_iter() const {	return const_iterator(*this);	}
		inline iterator iter(){	return iterator(*this);	}

		inline int begin() const 	{ return 0;	}
		inline int end()	const	{	return N;	}

		inline int being(int gb)	const	{	return 0;	}
		inline int end(int ge)	const	{	return N;	}

		inline int size() const {	return N;	}
		inline int length() const {	return N;	}

		T operator()(char wh){
			switch(wh){
				case 'x':
					return data_[0];
					break;
				case 'y':
					return data_[1];
					break;
				case 'z':
					return data_[2];
					break;
				case 't':
					return data_[3];
					break;
				default:
					return ZeroType<T>::zero();
					break;
			}
		}

		T & operator()(char wh) const {
			switch(wh){
				case 'x':
					return data_[0];
					break;
				case 'y':
					return data_[1];
					break;
				case 'z':
					return data_[2];
					break;
				case 't':
					return data_[3];
					break;
				default:
					return ZeroType<T>::zero();
					break;
			}
		}


		inline void operator()(T xx, T yy, T zz)  {
			data_[0]=xx;
			data_[1]=yy;
			data_[2]=zz;
		}

		inline void operator()(T xx, T yy)  {
			data_[0]=xx;
			data_[1]=yy;
		}

		inline void operator()(T xx, T yy, T zz, T tt)  {
			data_[0]=xx;
			data_[1]=yy;
			data_[2]=zz;
			data_[3]=tt;
		}

		inline T  x()const{ return data_[0];	}
		inline T & x(){ return data_[0];	}

		inline T y()const{ return data_[1]; }
		inline T &y()  { return data_[1];	}

		inline T z()const{ return  data_[2];	}
		inline T &z() { return data_[2];	}

		inline T  r()const{ return data_[0];	}
		inline T & r(){ return data_[0];	}

		inline T theta()const{ return data_[1];	}
		inline T &theta() { return data_[1];	}

		inline T phi()const { return data_[2];	}
		inline T &phi() { return data_[2];	}

		inline T t()const{ return data_[3];	}
		inline T &t() { return  data_[3] ;	}

		template<class T1> inline void x(const T1 &inx){ data_[0]=T1(inx);	}
		template<class T1> inline void y(const T1 &iny){ data_[1]=T1(iny);	}
		template<class T1> inline void z(const T1 &inz){ data_[2]=T1(inz);	}
		template<class T1> inline void r(const T1 &inx){ data_[0]=T1(inx);	}
		template<class T1> inline void theta(const T1 &iny){ data_[1]=T1(iny);	}
		template<class T1> inline void phi(const T1 &inz){ data_[2]=T1(inz);	}
		template<class T1> inline void t(const T1 &intt){ data_[3]=T1(intt);	}

		template<class T1> void set(T1 num, char wh){
			switch(wh){
				case 'x':
					data_[0]=T(num);
					break;
				case 'y':
					data_[1]=T(num);
					break;
				case 'z':
					data_[2]=T(num);
					break;
				case 't':
					data_[3]=T(num);
					break;
				default:
					break;
			}
		}

		template<class T2> inline coord<T,N>& operator= (const coord<T2,N> &rhs);
		template<class T2> inline coord<T,N>& operator+=(const coord<T2,N> &in);
		template<class T2> inline coord<T,N>& operator-=(const coord<T2,N> &in);
		template<class T2> inline coord<T,N>& operator*=(const coord<T2,N> &in);
		template<class T2> inline coord<T,N>& operator/=(const coord<T2,N> &in);

		template<class T1> inline coord<T,N>& operator= (const VExpr<T1> &rhs);
		template<class T2> inline coord<T,N>& operator+=(const VExpr<T2> &in);
		template<class T2> inline coord<T,N>& operator-=(const VExpr<T2> &in);
		template<class T2> inline coord<T,N>& operator*=(const VExpr<T2> &in);
		template<class T2> inline coord<T,N>& operator/=(const VExpr<T2> &in);

		template<class T2> inline coord<T,N>& operator= (const Vector<T2> &rhs);
		template<class T2> inline coord<T,N>& operator+=(const Vector<T2> &in);
		template<class T2> inline coord<T,N>& operator-=(const Vector<T2> &in);
		template<class T2> inline coord<T,N>& operator*=(const Vector<T2> &in);
		template<class T2> inline coord<T,N>& operator/=(const Vector<T2> &in);


		template<class T1> inline coord<T,N>& operator= (const T1 &rhs);
		template<class T1> inline coord<T,N>& operator+=(const T1 &in);
		template<class T1> inline coord<T,N>& operator-=(const T1 &in);
		template<class T1> inline coord<T,N>& operator*=(const T1 &in);
		template<class T1> inline coord<T,N>& operator/=(const T1 &in);



//aplies a 'function' across an entire vector
//The input class MUST have the operator(T, int) defined

		template<class Func>
		void apply( Func &f)
		{
			for(int i=0;i<N;++i)
			{
				(*this)(i)=f((*this)(i), i);
			}
		}

//aplies a 'function' across an a Range in the vector
//The input class MUST have the operator(T, int) defined

		template<class Func>
		void apply( Func &f, const Range &r)
		{
			for(int i=0;i<r.length(N);++i)
			{
				(*this)(r(i))=f(get(r(i)), r(i));
			}
		}



//Simple 3-D Euler angle roation of a 3D point!!! NOT 2D
		void Rotate3D(double theta, double phi=0., double psi=0.){
			T tx=x(), ty=y(), tz=z();
			double cth=cos(theta), sth=sin(theta), cps=cos(psi), sps=sin(psi), cph=cos(phi), sph=sin(phi);
			put(0,tx*(cps*cph-cth*sph*sps)+
				ty*(cps*sph+cth*cph*sps)+
				tz*(sps*sth));
			put(1,tx*(-sps*cph-cth*sph*cps)+
				ty*(-sps*sph+cth*cph*cps)+
				tz*(cps*sth));
			put(2,tx*(sth*sph)-
				ty*(sth*cph)+
				tz*cth);
		}

//Simple rotationg about the axis for a the angle theta
		void rotate(double theta,const coord<> &axi){
			T tx=x(), ty=y(), tz=z();
			double c=cos(theta), s=sin(theta), t=(1.0-c);

			put(0,tx*(c+axi.x()*axi.x()*t)+
				ty*(axi.x()*axi.y()*t+axi.z()*s)+
				tz*(axi.x()*axi.z()*t-axi.y()*s));

			put(1,tx*(axi.x()*axi.y()*t-axi.z()*s)+
				ty*(c+axi.y()*axi.y()*t)+
				tz*(axi.y()*axi.z()*t+axi.x()*s));

			put(2,tx*(axi.x()*axi.z()*t+axi.y()*s)+
				ty*(axi.y()*axi.z()*t-axi.x()*s)+
				tz*(c+axi.z()*axi.z()*t));

		}

//Simple 2-D Euler angle roation of a 2D point
		void Rotate2D(double theta){
			T tx=(*this)(0), ty=(*this)(1);
			(*this)(0)=tx*(cps*cph-cth*sph*sps)+
				ty*(cps*sph+cth*cph*cps)+
				tz*(sps*sth);
			(*this)(1)=tx*(-sps*cph-cth*sph*cps)+
				ty*(-sps*sph-cth*cph*cps)+
				tz*(cps*sth);
		}

		coord cyl2cart(){
			coord tm(*this);
			//if(TheType==cartesian) return tm;
			double cp=cos(phi());
			tm.x(r()*cos(theta()));
			tm.y(r()*sin(theta()));
			//tm.z(r()*sin(phi());
			//tm.TheType=cylindrical;
			return tm;
		}

		coord cart2cyl(){
			coord tm(*this);
			//if(TheType==cylindrical) return tm;
			double cp=cos(phi());
			tm.r(sqrt(x()*x()+y()*y()));
			tm.theta(fmod(atan(y()/x()), PI2));
			if(tm.theta()<0)	tm.theta(PI2+tm.theta());
			//tm.z(z());
			//tm.TheType=cylindrical;
			return tm;
		}

		coord cart2sph(){
			coord tm(*this);
			//if(TheType==spherical) return tm;
			tm.r(sqrt(x()*x()+y()*y()+z()*z()));

			tm.phi(fmod(asin( z()/tm.r() ), PI2) );
			if(tm.phi()<0)	tm.phi(PI2+tm.phi());

			tm.theta(acos(x()/(tm.r()*cos(tm.phi() ) ) ) );
			if(tm.theta()<0)	tm.theta(PI2+tm.theta());
			//TheType=spherical;
			return tm;
		}

		coord sph2cart(){
			coord tm(*this);
			//if(TheType==spherical) return tm;
			T sp=cos(phi());
			tm.x(r()*cos(theta())*sp);
			tm.y(r()*sin(theta()) * sp);
			tm.z(r()*sin(phi() ) );
			//tm.TheType=spherical;
			return tm;
		}

		void ToSphere(){
			//if(TheType==spherical) return;
			//else if(TheType==cartesian){
				T tx=x(), ty=y(), tz=z();
				r(sqrt(tx*tx+ty*ty+tz*tz));
				phi(asin(tz/r()));
				if(phi()<0)	phi(PI2+phi());

				theta(acos(tx/(r()*cos(phi() ) ) ) );
				if(theta()<0)	theta(PI2+theta());

			//	TheType=spherical;
			//	return;
			//}else if(TheType==cylindrical){
			//	ToCart();
			//	ToSphere();
			//}
		}

		void ToCart(){
			//if(TheType==cartesian) return;
			//else if(TheType==spherical){
				T tr=r(), tt=theta(), tp=phi(), cph=cos(tp);
				x(tr*cos(tt)*cph);
				y(tr * sin(tt) * cph);
				z(tr*sin(tp));
			//	TheType=cartesian;
			//}else if(TheType==cylindrical){
			//	T tr=r(), tt=theta();
			//	x(tr*cos(tt));
			//	y(tr*sin(tt));
			//	TheType=cartesian;
			//}
		}

		void ToCyl(){
			//if(TheType==cylindrical) return;
			//else if(TheType==cartesian){
				T tx=x(), ty=y();
				r(sqrt(tx*tx+ty*ty));
				theta(fmod(atan(ty/tx), PI2));
				if(theta()<0)	theta(PI2+theta());
				//tm.z(z());
			//	TheType=cylindrical;
			//}else if(TheType==spherical){
			//	ToCart();
			//	ToCyl();
			//}
		}
};

//-------ITERATORS-------------------------------------------

template<class T, int N>
class coordIter {
private:
	T *data_;
public:
    typedef T numtype;

    explicit coordIter( coord<T, N>& x)
    	: data_(x.data())
    { }

    coordIter(const coordIter<T, N>& iter)
    	: data_(iter.data_)
    { }

    inline T operator[](int i) const
    {
      //  CompTimeAssert(i < N);
        return data_[i];
    }

    inline T& operator[](int i)
    {
      //  CompTimeAssert(i < N);
        return data_[i];
    }

    T operator()(int i) const
    {
     //  CompTimeAssert(i < N);
        return data_[i];
    }

    T& operator()(int i)
    {
      //  CompTimeAssert(i < N);
        return data_[i];
    }

	inline int begin() const 	{ return 0;	}
	inline int end()	const	{	return N;	}

	inline int begin(int gb)	const	{	return 0;	}
	inline int end(int ge)	const	{	return N;	}

    inline int length(int ) const   { return N; }
    inline int guessLength() const    { return N; }
    inline int guessLength(int ) const    { return N; }
    inline int guessLen() const    { return N; }
    inline int guessLen(int ) const    { return N; }

};

template<class T, int N>
class coordIterConst {
private:
	const T *data_;
public:
    typedef T numtype;

    explicit coordIterConst(const coord<T, N>& x)
    	: data_(x.data())
    { }

    coordIterConst(const coordIterConst<T, N>& iter)
    	: data_(iter.data_)
    { }

    inline T operator[](int i) const
    {
      //  CompTimeAssert(i < N);
        return data_[i];
    }


    T operator()(int i) const
    {
        //CompTimeAssert(i < N);
        return data_[i];
    }
	inline int begin() const 	{ return 0;	}
	inline int end()	const	{	return N;	}

	inline int begin(int gb)	const	{	return 0;	}
	inline int end(int ge)	const	{	return N;	}


    inline int length(int ) const   { return N; }
    inline int guessLength() const    { return N; }
    inline int guessLength(int ) const    { return N; }
    inline int guessLen() const    { return N; }
    inline int guessLen(int ) const    { return N; }

};


//---------------------------------------------------------

//assignment meta programs...------------------------------
template<class T, int N>
template<class Expr, class Op>
inline void coord<T, N>::assign(Expr expr, Op up)
{
    vecAssign<N, 0>::assign(*this, expr, up);
}

// Copy Constructors.....----------------------------------

template<class T, int N>
template<class T1>
inline coord<T, N>::coord(const coord<T1,N>& x)
{
    //TheType=x.TheType;
    coord<T,N>::assign(VExpr<typename coord<T,N>::const_iterator>(x.begin()), ApAssign<T, T1>());
}


// constructor from VExpr
template<class T, int N>
template<class expr>
inline coord<T, N>::coord(const VExpr<expr> &in)
{
  coord<T, N>::assign(in, ApAssign<T, typename VExpr<expr>::numtype>());
}

//constructor from Vector
template<class T, int N>
template<class T1>
inline coord<T, N>::coord(const Vector<T1> &expr)
{
  coord<T, N>::assign(VExpr<Vector<T1> >(expr), ApAssign<T, typename VExpr<Vector<T1> >::numtype>());
}

//constructor from ptr data
template<class T, int N>
inline coord<T, N>::coord(T *data)
{
	for(int i=0;i<N;++i) data_[i]=data[i];
}
//---------------------------------------------------------



//assignments
template<class T, int N>
template<class T1>
inline coord<T,N> &coord<T,N>::operator=(const coord<T1,N> &rhs)
{
	//TheType=rhs.TheType;
    (*this) = VExpr<typename coord<T1, N>::const_iterator>(x.const_iter());
    return *this;
}

template<class T, int N>
template<class T1>
inline coord<T,N> &coord<T,N>::operator=(const VExpr<T1> &rhs)
{
	coord<T, N>::assign(rhs, ApAssign<T, typename VExpr<T1>::numtype>());
	return *this;
}

template<class T, int N>
template<class T1>
inline coord<T,N> &coord<T,N>::operator=(const Vector<T1> &rhs)
{
   	(*this) = VExpr<Vector<T1> >(rhs);
   	return *this;
}

template<class T, int N>
template<class T1>
inline coord<T,N> &coord<T,N>::operator=(const T1 &rhs)
{
   	(*this) = VExpr<VExprConst<T1> >(rhs);
   	return *this;
}

/*pass the 'binary' self operations -->coord OP coord<-- */
#define CoordBinSelfOp(OP) \
		template<class T, int N> template<class T2> \
		inline coord<T,N>& coord<T,N>::operator OP (const coord<T2, N> &x){	\
			(*this) OP VExpr<typename coord<T2, N>::const_iterator>(x.const_iter());	\
			 return *this;	\
		}	\

/*pass the 'binary' self operations -->coord OP Number<-- */
#define CoordBinSelfOpNum(OP) \
		template<class T, int N> template<class T1> \
		inline coord<T,N>& coord<T,N>::operator OP (const T1 &x){	\
			(*this) OP VExpr< VExprConst<T1> >( VExpr< VExprConst<T1> >(x) );	\
			 return *this;	\
		}	\

//Vector bin ops -->coord OP Vector<--
#define CoordBinSelfOpVec(OP) \
		template<class T, int N> template<class T1> \
		inline coord<T,N>& coord<T,N>::operator OP (const Vector<T1> &in)	\
		{	\
			(*this) OP VExpr<typename coord<T1, N>::const_iterator>(x.const_iter());	\
			return *this;	\
		}	\

//VExpr bin ops.... -->coord OP VExpr<--
#define CoordBinSelfOpVexpr(OP, OpName) \
		template<class T, int N> template<class T2> \
		inline coord<T,N>& coord<T,N>::operator OP (const VExpr<T2> &in)	\
		{	\
			coord<T,N>::assign(in, OpName<T, typename T2::numtype>());	\
			return *this;	\
		}	\

/*
//for the builtin class types...(double, float, etc)
// avoids the 'ambigous over load' error
#define CoordBinSelfOp2BuiltIn(WHAT, OP) \
		template<class T, int N> \
		inline coord<T,N>& coord<T,N>::operator OP (const WHAT &in){	\
			Vector<T>::operator OP (in); return *this;	\
		}	\
*/


CoordBinSelfOpVexpr( +=, ApAdd2)
CoordBinSelfOpVexpr( -=, ApSubtract2)
CoordBinSelfOpVexpr( *=, ApMultiply2)
CoordBinSelfOpVexpr( /=, ApDivide2)

CoordBinSelfOp(+=)
CoordBinSelfOp(-=)
CoordBinSelfOp(*=)
CoordBinSelfOp(/=)

CoordBinSelfOpVec(+=)
CoordBinSelfOpVec(-=)
CoordBinSelfOpVec(*=)
CoordBinSelfOpVec(/=)

CoordBinSelfOpNum(+=)
CoordBinSelfOpNum(-=)
CoordBinSelfOpNum(*=)
CoordBinSelfOpNum(/=)

/*
CoordBinSelfOp2BuiltIn(complex, +=)
CoordBinSelfOp2BuiltIn(double, +=)
CoordBinSelfOp2BuiltIn(float, +=)
CoordBinSelfOp2BuiltIn(short, +=)
CoordBinSelfOp2BuiltIn(int, +=)
CoordBinSelfOp2BuiltIn(long, +=)
CoordBinSelfOp2BuiltIn(char, +=)
CoordBinSelfOp2BuiltIn(bool, +=)

CoordBinSelfOp2BuiltIn(complex, -=)
CoordBinSelfOp2BuiltIn(double, -=)
CoordBinSelfOp2BuiltIn(float, -=)
CoordBinSelfOp2BuiltIn(short, -=)
CoordBinSelfOp2BuiltIn(int, -=)
CoordBinSelfOp2BuiltIn(long, -=)
CoordBinSelfOp2BuiltIn(char, -=)
CoordBinSelfOp2BuiltIn(bool, -=)

CoordBinSelfOp2BuiltIn(complex, /=)
CoordBinSelfOp2BuiltIn(double, /=)
CoordBinSelfOp2BuiltIn(float, /=)
CoordBinSelfOp2BuiltIn(short, /=)
CoordBinSelfOp2BuiltIn(int, /=)
CoordBinSelfOp2BuiltIn(long, /=)
CoordBinSelfOp2BuiltIn(char, /=)
CoordBinSelfOp2BuiltIn(bool, /=)

CoordBinSelfOp2BuiltIn(complex, *=)
CoordBinSelfOp2BuiltIn(double, *=)
CoordBinSelfOp2BuiltIn(float, *=)
CoordBinSelfOp2BuiltIn(short, *=)
CoordBinSelfOp2BuiltIn(int, *=)
CoordBinSelfOp2BuiltIn(long, *=)
CoordBinSelfOp2BuiltIn(char, *=)
CoordBinSelfOp2BuiltIn(bool, *=)
*/
/* Here we pass the 'binary' operators of a coordinate operation to the VEctor class
   this clears up the ambiguity of the orlaoded function..
   it behaves like the Vector versions...*/


#define CoordBinOpVec(LHS,N1, RHS, N2,NAME, OP)	\
	template<class T1, class T2, int N> \
	inline VExpr<VBinExprOp<coordIterConst<T1,N>,coordIterConst<T2,N>,NAME<T1, T2> > >   operator OP (const LHS , N1 & a, const RHS , N1 & b) \
	{		\
		typedef VBinExprOp<coordIterConst<T1,N>,coordIterConst<T2,N>, NAME<T1,T2> > ExprT;		\
		return VExpr<ExprT>(ExprT(a.const_iter(),b.const_iter()));			\
	}		\

#define CoordBinOpVec2(LHS,N1, RHS, NAME, OP)	\
	template<class T1, class T2, int N> \
	inline VExpr<VBinExprOp<coordIterConst<T1,N>,VIterConst<T2>,NAME<T1, T2> > >   operator OP (const LHS , N1 & a, const RHS  & b) \
	{		\
		typedef VBinExprOp<coordIterConst<T1,N>,VIterConst<T2>, NAME<T1,T2> > ExprT;		\
		return VExpr<ExprT>(ExprT(a.const_iter(),b.const_iter()));			\
	}	\


#define CoordBinOpVec3(LHS, RHS,N2, NAME, OP)	\
	template<class T1, class T2, int N> \
	inline VExpr<VBinExprOp<VIterConst<T1>,coordIterConst<T2,N>,NAME<T1, T2> > >   operator OP (const LHS  & a, const RHS , N2 & b) \
	{		\
		typedef VBinExprOp<VIterConst<T1>,coordIterConst<T2,N>, NAME<T1,T2> > ExprT;		\
		return VExpr<ExprT>(ExprT(a.const_iter(),b.const_iter()));			\
	}		\

CoordBinOpVec(coord<T1,N>, coord<T2,N>, ApMultiply, *)
CoordBinOpVec(coord<T1,N>, coord<T2,N>, ApAdd, +)
CoordBinOpVec(coord<T1,N>, coord<T2,N>, ApSubtract, -)
CoordBinOpVec(coord<T1,N>, coord<T2,N>, ApDivide, /)

CoordBinOpVec2(coord<T1,N>, Vector<T2>,  ApMultiply, *)
CoordBinOpVec2(coord<T1,N>, Vector<T2>,  ApAdd, +)
CoordBinOpVec2(coord<T1,N>, Vector<T2>,  ApSubtract, -)
CoordBinOpVec2(coord<T1,N>, Vector<T2>,  ApDivide, /)

CoordBinOpVec3(Vector<T1>, coord<T2,N>, ApMultiply, *)
CoordBinOpVec3(Vector<T1>, coord<T2,N>, ApAdd, +)
CoordBinOpVec3(Vector<T1>, coord<T2,N>, ApSubtract, -)
CoordBinOpVec3(Vector<T1>, coord<T2,N>, ApDivide, /)


#define CoordBinOpExprRhs(RHS, NAME, OP)	\
	template<class T1, class T2, int N> \
	inline VExpr<VBinExprOp<VExpr<T1>,coordIterConst<T2,N>,NAME<typename T1::numtype, T2> > >   operator OP (const VExpr<T1> & a, const RHS ,N> & b) \
	{		\
		typedef VBinExprOp<VExpr<T1>,coordIterConst<T2,N>, NAME<typename T1::numtype,T2> > ExprT;		\
		return VExpr<ExprT>(ExprT(a,b.const_iter()));			\
	}		\

CoordBinOpExprRhs(coord<T2, ApMultiply, *)
CoordBinOpExprRhs(coord<T2, ApAdd, +)
CoordBinOpExprRhs(coord<T2, ApSubtract, -)
CoordBinOpExprRhs(coord<T2, ApDivide, /)



#define CoordBinOpExprLhs(LHS, NAME, OP)	\
	template<class T1, class T2, int N> \
	inline VExpr<VBinExprOp<coordIterConst<T1,N>,VExpr<T2>,NAME<T1, typename T2::numtype> > >   operator OP (const LHS ,N> & a, const VExpr<T2> & b) \
	{		\
		typedef VBinExprOp<coordIterConst<T1,N>,VExpr<T2>, NAME<T1,typename T2::numtype> > ExprT;		\
		return VExpr<ExprT>(ExprT(a.const_iter(),b));			\
	}		\

CoordBinOpExprLhs(coord<T1, ApMultiply, *)
CoordBinOpExprLhs(coord<T1, ApAdd, +)
CoordBinOpExprLhs(coord<T1, ApSubtract, -)
CoordBinOpExprLhs(coord<T1, ApDivide, /)



//for builtin types

#define CoordBinOpConstLhsBuiltIn(RHS, LHS, NAME, OP)	\
	template< class T1,  int N> \
	inline VExpr<VBinExprOp<coordIterConst<T1,N>,VExprConst<RHS>,NAME<T1, RHS> > >   operator OP (const LHS ,N> & a, const RHS & b) \
	{		\
		typedef VBinExprOp<coordIterConst<T1,N>,VExprConst<RHS>, NAME<T1,RHS> > ExprT;		\
		return VExpr<ExprT>(ExprT(a.const_iter(),VExprConst<RHS>(b)));			\
	}		\

CoordBinOpConstLhsBuiltIn(complex,coord<T1, ApMultiply, *)
CoordBinOpConstLhsBuiltIn(complex,coord<T1, ApAdd, +)
CoordBinOpConstLhsBuiltIn(complex,coord<T1, ApSubtract, -)
CoordBinOpConstLhsBuiltIn(complex,coord<T1, ApDivide, /)

CoordBinOpConstLhsBuiltIn(double,coord<T1, ApMultiply, *)
CoordBinOpConstLhsBuiltIn(double,coord<T1, ApAdd, +)
CoordBinOpConstLhsBuiltIn(double,coord<T1, ApSubtract, -)
CoordBinOpConstLhsBuiltIn(double,coord<T1, ApDivide, /)

CoordBinOpConstLhsBuiltIn(float,coord<T1, ApMultiply, *)
CoordBinOpConstLhsBuiltIn(float,coord<T1, ApAdd, +)
CoordBinOpConstLhsBuiltIn(float,coord<T1, ApSubtract, -)
CoordBinOpConstLhsBuiltIn(float,coord<T1, ApDivide, /)

CoordBinOpConstLhsBuiltIn(int,coord<T1, ApMultiply, *)
CoordBinOpConstLhsBuiltIn(int,coord<T1, ApAdd, +)
CoordBinOpConstLhsBuiltIn(int,coord<T1, ApSubtract, -)
CoordBinOpConstLhsBuiltIn(int,coord<T1, ApDivide, /)

CoordBinOpConstLhsBuiltIn(short,coord<T1, ApMultiply, *)
CoordBinOpConstLhsBuiltIn(short,coord<T1, ApAdd, +)
CoordBinOpConstLhsBuiltIn(short,coord<T1, ApSubtract, -)
CoordBinOpConstLhsBuiltIn(short,coord<T1, ApDivide, /)

CoordBinOpConstLhsBuiltIn(long,coord<T1, ApMultiply, *)
CoordBinOpConstLhsBuiltIn(long,coord<T1, ApAdd, +)
CoordBinOpConstLhsBuiltIn(long,coord<T1, ApSubtract, -)
CoordBinOpConstLhsBuiltIn(long,coord<T1, ApDivide, /)

CoordBinOpConstLhsBuiltIn(char,coord<T1, ApMultiply, *)
CoordBinOpConstLhsBuiltIn(char,coord<T1, ApAdd, +)
CoordBinOpConstLhsBuiltIn(char,coord<T1, ApSubtract, -)
CoordBinOpConstLhsBuiltIn(char,coord<T1, ApDivide, /)

CoordBinOpConstLhsBuiltIn(bool,coord<T1, ApMultiply, *)
CoordBinOpConstLhsBuiltIn(bool,coord<T1, ApAdd, +)
CoordBinOpConstLhsBuiltIn(bool,coord<T1, ApSubtract, -)
CoordBinOpConstLhsBuiltIn(bool,coord<T1, ApDivide, /)

#define CoordBinOpConstRhsBuiltIn(LHS, RHS, NAME, OP)	\
	template< class T2, int N> \
	inline VExpr<VBinExprOp<VExprConst<LHS>,coordIterConst<T2,N>,NAME<LHS, T2> > >   operator OP (const LHS & a, const RHS ,N> & b) \
	{		\
		typedef VBinExprOp<VExprConst<LHS>,coordIterConst<T2,N>, NAME<LHS,T2> > ExprT;		\
		return VExpr<ExprT>(ExprT(VExprConst<LHS>(a),b.const_iter()));			\
	}		\

CoordBinOpConstRhsBuiltIn(complex,coord<T2, ApMultiply, *)
CoordBinOpConstRhsBuiltIn(complex,coord<T2, ApAdd, +)
CoordBinOpConstRhsBuiltIn(complex,coord<T2, ApSubtract, -)
CoordBinOpConstRhsBuiltIn(complex,coord<T2, ApDivide, /)

CoordBinOpConstRhsBuiltIn(double,coord<T2, ApMultiply, *)
CoordBinOpConstRhsBuiltIn(double,coord<T2, ApAdd, +)
CoordBinOpConstRhsBuiltIn(double,coord<T2, ApSubtract, -)
CoordBinOpConstRhsBuiltIn(double,coord<T2, ApDivide, /)

CoordBinOpConstRhsBuiltIn(float,coord<T2, ApMultiply, *)
CoordBinOpConstRhsBuiltIn(float,coord<T2, ApAdd, +)
CoordBinOpConstRhsBuiltIn(float,coord<T2, ApSubtract, -)
CoordBinOpConstRhsBuiltIn(float,coord<T2, ApDivide, /)

CoordBinOpConstRhsBuiltIn(int,coord<T2, ApMultiply, *)
CoordBinOpConstRhsBuiltIn(int,coord<T2, ApAdd, +)
CoordBinOpConstRhsBuiltIn(int,coord<T2, ApSubtract, -)
CoordBinOpConstRhsBuiltIn(int,coord<T2, ApDivide, /)

CoordBinOpConstRhsBuiltIn(short,coord<T2, ApMultiply, *)
CoordBinOpConstRhsBuiltIn(short,coord<T2, ApAdd, +)
CoordBinOpConstRhsBuiltIn(short,coord<T2, ApSubtract, -)
CoordBinOpConstRhsBuiltIn(short,coord<T2, ApDivide, /)

CoordBinOpConstRhsBuiltIn(long,coord<T2, ApMultiply, *)
CoordBinOpConstRhsBuiltIn(long,coord<T2, ApAdd, +)
CoordBinOpConstRhsBuiltIn(long,coord<T2, ApSubtract, -)
CoordBinOpConstRhsBuiltIn(long,coord<T2, ApDivide, /)

CoordBinOpConstRhsBuiltIn(char,coord<T2, ApMultiply, *)
CoordBinOpConstRhsBuiltIn(char,coord<T2, ApAdd, +)
CoordBinOpConstRhsBuiltIn(char,coord<T2, ApSubtract, -)
CoordBinOpConstRhsBuiltIn(char,coord<T2, ApDivide, /)

CoordBinOpConstRhsBuiltIn(bool,coord<T2, ApMultiply, *)
CoordBinOpConstRhsBuiltIn(bool,coord<T2, ApAdd, +)
CoordBinOpConstRhsBuiltIn(bool,coord<T2, ApSubtract, -)
CoordBinOpConstRhsBuiltIn(bool,coord<T2, ApDivide, /)


/* Comparison Operations */


template<class T1, class T2, int N>
inline bool  operator == ( coord<T1,N> & a,  coord<T2 ,N> & b)
{
	for(int i=0;i<N;++i)
		if(a[i] != b[i]) return false;

	return true;
}

template<class T1, class T2, int N>
inline bool  operator != ( coord<T1,N> & a,  coord<T2 ,N> & b)
{
	for(int i=0;i<N;++i)
		if(a[i] == b[i]) return false;

	return true;
}

template<class T1, class T2, int N>
inline bool  operator >= ( coord<T1,N> & a,  coord<T2 ,N> & b)
{
	for(int i=0;i<N;++i)
		if(a[i] < b[i]) return false;

	return true;
}

template<class T1, class T2, int N>
inline bool  operator <= ( coord<T1,N> & a,  coord<T2 ,N> & b)
{
	for(int i=0;i<N;++i)
		if(a[i] > b[i]) return false;

	return true;
}

template<class T1, class T2, int N>
inline bool  operator < ( coord<T1,N> & a,  coord<T2 ,N> & b)
{
	for(int i=0;i<N;++i)
		if(a[i] >= b[i]) return false;

	return true;
}

template<class T1, class T2, int N>
inline bool  operator > ( coord<T1,N> & a,  coord<T2 ,N> & b)
{
	for(int i=0;i<N;++i)
		if(a[i] <= b[i]) return false;

	return true;
}

//comparisons for number types...
#define CoordCompNum(NUM) \
	template<class T1,  int N>	\
	inline bool  operator == (const coord<T1,N> & a, const NUM & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] != b) return false; 	\
		return true;	\
	} \
	template<class T1,  int N>	\
	inline bool  operator != (const coord<T1,N> & a, const NUM & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] == b) return false; \
		return true;	\
	} \
	template<class T1,  int N>	\
	inline bool  operator >= (const coord<T1,N> & a, const NUM & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] < b) return false; 	\
		return true;	\
	} \
	template<class T1,  int N>	\
	inline bool  operator <= (const coord<T1,N> & a, const NUM & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] > b) return false; \
		return true;	\
	} \
	template<class T1,  int N>	\
	inline bool  operator < (const coord<T1,N> & a, const NUM & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] >= b) return false; 	\
		return true;	\
	} \
	template<class T1, int N>	\
	inline bool  operator > (const coord<T1,N> & a, const NUM & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] <= b) return false; 	\
		return true;	\
	}	\
	\
	template<class T1,  int N>	\
	inline bool  operator == (const NUM & a, const coord<T1,N> & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] != b) return false; 	\
		return true;	\
	} \
	template<class T1,  int N>	\
	inline bool  operator != (const NUM & a, const coord<T1,N> & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] == b) return false; \
		return true;	\
	} \
	template<class T1,  int N>	\
	inline bool  operator >= (const NUM & a, const coord<T1,N> & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] < b) return false; 	\
		return true;	\
	} \
	template<class T1,  int N>	\
	inline bool  operator <= (const NUM & a, const coord<T1,N> & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] > b) return false; \
		return true;	\
	} \
	template<class T1,  int N>	\
	inline bool  operator < (const NUM & a, const coord<T1,N> & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] >= b) return false; 	\
		return true;	\
	} \
	template<class T1, int N>	\
	inline bool  operator > (const NUM & a, const coord<T1,N> & b)	\
	{	\
		for(int i=0;i<N;++i)	\
			if(a[i] <= b) return false; 	\
		return true;	\
	}	\

CoordCompNum( double)
CoordCompNum( float)
CoordCompNum( int)
CoordCompNum( unsigned int)
CoordCompNum( complex)
CoordCompNum( short)
CoordCompNum( bool)
CoordCompNum( char)
CoordCompNum( unsigned char)

// MAtrix Coord functions

//
//// These Are needed for the optimized version of the matrix multiplication
//

 /************8 SUB MAtrix * coord Template functions */
template<class T, class INstructure, int N, class T2>
SumType(OutType(T,T2)) _mul_mat_vec_(const _matrixReference<T, INstructure> &a,const coordIterConst<T2,N> &v,  int row)
{
	SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	for(int k=0;k<a.cols(0);k++)	 tmp+=a(row,k) * v(k);
	return tmp;
}

template<class T,int N, class T2>
SumType(OutType(T,T2)) _mul_mat_vec_(const _matrixReference<T, DiagonalMatrix> &a,const coordIterConst<T2,N> &v,  int row)
{
	//SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	//tmp+=a(i,k) * v(k);
	return a(row, row)*v(row);
}

template<class T, int N,class T2>
SumType(OutType(T,T2)) _mul_mat_vec_(const _matrixReference<T, IdentityMatrix> &a,const coordIterConst<T2,N>  &v,  int row)
{
	//SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	//tmp+=a(i,k) * v(k);
	return v(row);
}

template<class T, int N,class T2>
SumType(OutType(T,T2)) _mul_mat_vec_(const _matrixReference<T, TriDiagonalMatrix> &a,const coordIterConst<T2,N>  &v,  int row)
{
	//SumType(OutType(T,T2)) tmp=ZeroType<SumType(OutType(T,T2))>::zero();
	//tmp+=a(i,k) * v(k);
	if(row==0 && a.rows(v.size())>1){	return a(0,0)*v(0)+a(0,1)*v(1);	}
	else if(row==0){	return a(0,0)*v(0);	}

	if(row==a.rows(v.size())-1 && a.rows(v.size())-1>=0){	return a(row,row)*v(row)+a(row,row-1)*v(row-1);	}
	else if(row==a.rows(v.size())-1){	return a(row,row)*v(row);	}

	return a(row, row-1)*v(row-1)+a(row,row)*v(row)+a(row, row+1)*v(row+1);
}


// Matrix*coord
template<class T, class INstructure, class T2, int N>
inline VExpr<VExprBinOpMatMul< _matrixReference<T, INstructure>, coordIterConst<T2, N> > >
operator*(const _matrix<T, INstructure>& m1,const coord<T2,N>& v1)
{
	typedef VExprBinOpMatMul< _matrixReference<T, INstructure>, coordIterConst<T2,N> > expr;
	return VExpr<expr>(expr(m1.getRef(),v1.const_iter()));
}

//coord*matrix
 template<class T, class INstructure, class T2, int N>
 inline VExpr<VExprBinOpMatMul2<coordIterConst<T2,N>, _matrixReference<T, INstructure> > >
 operator*(const coord<T2,N>& v1,const _matrix<T, INstructure>& m1)
 {
     typedef VExprBinOpMatMul2<coordIterConst<T2, N> , _matrixReference<T, INstructure> > expr;
     return VExpr<expr>(expr(v1.const_iter(),m1.getRef()));
 }

/* Here we pass the 'unary' operators of a coordinate operation to the VEctor class
   this clears up the ambiguity of the orlaoded function..
   it behaves like the Vector versions...*/

#define CoordUnaryOp(name, op)														\
	template<class T1, int N>																	\
	inline VExpr<VBinUnaryOp<coordIterConst<T1,N>, name<T1> > > 	op(const coord<T1,N>& a)				\
	{																					\
		typedef VBinUnaryOp<coordIterConst<T1,N>,name<T1> > ExprT;							\
		return VExpr<ExprT>(ExprT(a.const_iter()));											\
	}		\


CoordUnaryOp(ApAbs, abs)
CoordUnaryOp(ApCos, cos)
CoordUnaryOp(ApSin, sin)
CoordUnaryOp(ApTan, tan)
CoordUnaryOp(ApExp, exp)
CoordUnaryOp(ApAcos, acos)
CoordUnaryOp(ApAsin, asin)
CoordUnaryOp(ApAtan, atan)
CoordUnaryOp(ApLog, log)
CoordUnaryOp(ApCosh, cosh)
CoordUnaryOp(ApSinh, sinh)
CoordUnaryOp(ApTanh, tanh)
CoordUnaryOp(ApAcosh, acosh)
CoordUnaryOp(ApAsinh, asinh)
CoordUnaryOp(ApAtanh, atanh)
CoordUnaryOp(ApCiel, ceil)
CoordUnaryOp(ApFloor, floor)
#ifdef HAVE_FINITE
CoordUnaryOp(ApNan, isnan)
#elif HAVE_ISNAN
CoordUnaryOp(ApNan, isnan)
#endif
CoordUnaryOp(ApSqrt, sqrt)
CoordUnaryOp(ApSqr, sqr)

#ifdef HAVE_ISNAN
	template<class T1, int N>
	inline bool hasnan(const coord<T1,N> &b)
	{
		for(int i=0;i<N;++i)
			if(hasnan(b[i])) return true ;
		return false;
	}

#endif


//special operator for the 'negative'
template<class T1, int N>
inline VExpr<VBinUnaryOp<coordIterConst<T1,N>, ApNeg<T1> > >
operator-(const coord<T1,N>& a)
{
	typedef VBinUnaryOp<coordIterConst<T1,N>,ApNeg<T1> > ExprT;
	return VExpr<ExprT>(ExprT(a.const_iter()));
}

/*************** NORM **********************/
template<class T1>
inline FloatType(T1) ApNorm(const coord<T1, 1>& a)
{
	return a[0];
}

template<class T1>
inline FloatType(T1) ApNorm(const coord<T1, 2>& a)
{
	return pow(a[0]*a[0]+a[1]*a[1], 0.5);
}

template<class T1>
inline FloatType(T1) ApNorm(const coord<T1, 3>& a)
{
	return pow(a[0]*a[0]+a[1]*a[1]+a[2]*a[2], 0.5);
}

template<class T1>
inline FloatType(T1) ApNorm(const coord<T1, 4>& a)
{
	return pow(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]+a[3]*a[3], 0.5);
}

template<class T1, int N>
inline FloatType(T1) ApNorm(const coord<T1, N>& a)
{
	FloatType(T1) ss=a[0]*a[0];
	for(int i=1;i<N;++i) ss+=a[i]*a[i];
	return pow(ss, 0.5);
}


template<class T1,int N>
inline FloatType(T1)
norm(const coord<T1, N>& a)
{
	return ApNorm(a);
}

/*************** Square NORM **********************/
template<class T1>
inline FloatType(T1) ApSqNorm(const coord<T1, 1>& a)
{
	return a[0];
}

template<class T1>
inline FloatType(T1) ApSqNorm(const coord<T1, 2>& a)
{
	return a[0]*a[0]+a[1]*a[1];
}

template<class T1>
inline FloatType(T1) ApSqNorm(const coord<T1, 3>& a)
{
	return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
}

template<class T1>
inline FloatType(T1) ApSqNorm(const coord<T1, 4>& a)
{
	return a[0]*a[0]+a[1]*a[1]+a[2]*a[2]+a[3]*a[3];
}

template<class T1, int N>
inline FloatType(T1) ApSqNorm(const coord<T1, N>& a)
{
	FloatType(T1) ss=a[0]*a[0];
	for(int i=1;i<N;++i) ss+=a[i]*a[i];
	return ss;
}


template<class T1,int N>
inline FloatType(T1)
square_norm(const coord<T1, N>& a)
{
	return ApSqNorm(a);
}


//******************Min***********************

template<class T1, int N>
inline T1 ApMin(const coord<T1, N>& a){
	T1 ss=a[0];
	for(int i=1;i<N;i++){
		if(ss>=a[i]) ss=a[i];
	}
	return ss;
}

template<class T1, int N>
inline T1 min(const coord<T1,N>& a){
	return ApMin(a);
}

template<class T, int N>
inline T min(T in1, const coord<T, N> &cc){
	T tmp=min(cc);
	return in1<tmp?in1:tmp;
}

template<class T, int N>
inline T min(const coord<T, N> &cc, T in1 ){
	T tmp=min(cc);
	return in1<tmp?in1:tmp;
}
//******************Max***********************

template<class T1, int N>
inline T1 ApMax(const coord<T1, N>& a){
	T1 ss=a[0];
	for(int i=1;i<N;i++){
		if(ss<=a[i]) ss=a[i];
	}
	return ss;
}

template<class T1, int N>
inline T1 max(const coord<T1,N>& a){
	return ApMax(a);
}

template<class T, int N>
inline T max(T in1, const coord<T, N> &cc){
	T tmp=max(cc);
	return in1>tmp?in1:tmp;
}

template<class T, int N>
inline T max(const coord<T, N> &cc, T in1){
	T tmp=max(cc);
	return in1>tmp?in1:tmp;
}

//***********************dot *********************
template<class T1, class T2>
inline OutType(FloatType(T1),FloatType(T2)) ApDot(const coord<T1, 1>& a, const coord<T2, 1> &b)
{
	return a[0]*b[0];
}

template<class T1, class T2>
inline OutType(FloatType(T1),FloatType(T2)) ApDot(const coord<T1, 2>& a, const coord<T2, 2> &b)
{
	return a[0]*b[0]+a[1]*b[1];
}

template<class T1, class T2>
inline OutType(FloatType(T1),FloatType(T2)) ApDot(const coord<T1, 3>& a, const coord<T2, 3> &b)
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

template<class T1, class T2>
inline OutType(FloatType(T1),FloatType(T2)) ApDot(const coord<T1, 4>& a, const coord<T2, 4> &b)
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
}

template<class T1, class T2, int N>
inline OutType(FloatType(T1),FloatType(T2)) ApDot(const coord<T1, N>& a, const coord<T2, N> &b)
{
	OutType(FloatType(T1),FloatType(T2)) ss=a[0]*b[0];
	for(int i=1;i<N;i++) ss+=a[i]*b[i];
	return ss;
}

template<class T1, class T2, int N>
inline OutType(FloatType(T1),FloatType(T2)) ApDot(const coord<T1, N>& a, const Vector<T2> &b)
{
	OutType(FloatType(T1),FloatType(T2)) ss=a[0]*b[0];
	for(int i=1;i<N;i++) ss+=a[i]*b[i];
	return ss;
}

template<class T1, class T2, int N>
inline OutType(FloatType(T1),FloatType(T2)) ApDot(const Vector<T1>& a, const coord<T2, N> &b)
{
	OutType(FloatType(T1),FloatType(T2)) ss=a[0]*b[0];
	for(int i=1;i<N;i++) ss+=a[i]*b[i];
	return ss;
}

template<class T1, class expr, int N>
inline OutType(FloatType(T1),FloatType(typename VExpr<expr>::numtype))
ApDot(const coord<T1, N>& a, const VExpr<expr> &b)
{
	OutType(FloatType(T1),FloatType(typename VExpr<expr>::numtype))  ss=a[0]*b[0];
	for(int i=1;i<N;i++) ss+=a[i]*b[i];
	return ss;
}

template<class expr, class T2, int N>
inline OutType(FloatType(typename VExpr<expr>::numtype),FloatType(T2))
ApDot(const VExpr<expr>& a, const coord<T2, N> &b)
{
	OutType(FloatType(typename VExpr<expr>::numtype),FloatType(T2)) ss=a[0]*b[0];
	for(int i=1;i<N;i++) ss+=a[i]*b[i];
	return ss;
}

template<class T1, class T2, int N>
inline OutType(FloatType(T1),FloatType(T2))
dot(const coord<T1, N>& a, const coord<T2, N> &b)
{
	return ApDot(a,b);
}

template<class expr, class T2, int N>
inline OutType(FloatType(typename VExpr<expr>::numtype),FloatType(T2))
dot(const VExpr<expr>& a, const coord<T2, N> &b)
{
	return ApDot(a,b);
}

template<class T1, class expr, int N>
inline OutType(FloatType(T1),FloatType(typename VExpr<expr>::numtype))
dot(const coord<T1, N>& a, const VExpr<expr> &b)
{
	return ApDot(a,b);
}

template<class T1, class T2, int N>
inline OutType(FloatType(T1),FloatType(T2))
dot(const Vector<T1>& a, const coord<T2, N> &b)
{
	return ApDot(a,b);
}

template<class T1, class T2, int N>
inline OutType(FloatType(T1),FloatType(T2))
dot(const coord<T1,N>& a, const Vector<T2> &b)
{
	return ApDot(a,b);
}
//*********************** cross *********************
template<class T1>
inline coord<T1,1>
	ApCross(const coord<T1, 1>& a)
{
	return a;
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)),1>
	ApCross(const coord<T1, 1>& a, const coord<T2, 1> &b) 
{
	BLEXCEPTION(std::string(" cross product (asymetric multiply) is NOT allowed for ") + 
		std::string("\n items of vector length '1' "))
	return ZeroType<coord<OutType(FloatType(T1),FloatType(T2)),1> > ::zero();
}

template<class T1>
inline coord<FloatType(T1),2>
	ApCross(const coord<T1, 2>& a)
{
	typedef coord<FloatType(T1),2> OutT;
	return OutT(a[1], -a[0]);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)),2>
	ApCross(const coord<T1, 2>& a, const coord<T2, 2> &b) 
{
	BLEXCEPTION(std::string(" cross product (asymetric multiply) is NOT allowed for ") + 
		std::string("\n items of vector length '2' "))
	return ZeroType<coord<OutType(FloatType(T1),FloatType(T2)),2> > ::zero();
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)), 3>
ApCross(const coord<T1, 3>& a, const coord<T2, 3> &b)
{
	typedef coord<OutType(FloatType(T1),FloatType(T2)),3> OutT;
	return OutT(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}


template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)), 3>
ApCross(const coord<T1, 3>& a, const Vector<T2> &b)
{
	typedef coord<OutType(FloatType(T1),FloatType(T2)),3> OutT;
	return OutT(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)), 3>
ApCross(const Vector<T1>& a, const coord<T1, 3> &b)
{
	typedef coord<OutType(FloatType(T1),FloatType(T2)),3> OutT;
	return OutT(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)), 3>
ApCross(const coord<T1, 3>& a, const VExpr<T2> &b)
{
	typedef coord<FloatType(T1),3> OutT;
	return OutT(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)), 3>
ApCross(const VExpr<T1>& a, const coord<T1, 3> &b)
{
	typedef coord<FloatType(T2),3> OutT;
	return OutT(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
ApCross(const coord<T1, 4>& a, const coord<T2, 4> &b,const coord<T3, 4> &c)
{

	typedef coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}


template<class T1, class T2, class Expr>
inline coord<OutType(FloatType(Expr),OutType(FloatType(T1),FloatType(T2))), 4> ApCross(const coord<T1, 4>& a, const coord<T2, 4> &b,const VExpr<Expr> &c)
{

	typedef coord<OutType(FloatType(Expr),OutType(FloatType(T1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}

template<class T1, class Expr1, class Expr2>
inline coord<OutType(FloatType(Expr2),OutType(FloatType(T1),FloatType(Expr1))), 4>
  ApCross(const coord<T1, 4>& a, const VExpr<Expr1> &b,const VExpr<Expr2> &c)
{

	typedef coord<OutType(FloatType(Expr2),OutType(FloatType(T1),FloatType(Expr1))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}


template<class Expr1, class T2, class Expr2>
inline coord<OutType(FloatType(Expr2),OutType(FloatType(Expr1),FloatType(T2))), 4>
  ApCross( const VExpr<Expr1> &a, const coord<T2, 4>& b,const VExpr<Expr2> &c)
{

	typedef coord<OutType(FloatType(Expr2),OutType(FloatType(Expr1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4> ApCross(const coord<T1, 4>& a, const VExpr<T2> &b,const coord<T3, 4> &c)
{

	typedef coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4> ApCross(const VExpr<T1>& a, const coord<T2, 4> &b,const coord<T3, 4> &c)
{

	typedef coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}


template<class T1, class T2, class Expr>
inline coord<OutType(FloatType(Expr),OutType(FloatType(T1),FloatType(T2))), 4>
ApCross(const coord<T1, 4>& a, const coord<T2, 4> &b,const Vector<Expr> &c)
{

	typedef coord<OutType(FloatType(Expr),OutType(FloatType(T1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}

template<class T1, class Expr1, class Expr2>
inline coord<OutType(FloatType(Expr2),OutType(FloatType(T1),FloatType(Expr1))), 4> ApCross(const coord<T1, 4>& a, const Vector<Expr1> &b,const Vector<Expr2> &c)
{

	typedef coord<OutType(FloatType(Expr2),OutType(FloatType(T1),FloatType(Expr1))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}


template<class Expr1, class T2, class Expr2>
inline coord<OutType(FloatType(Expr2),OutType(FloatType(Expr1),FloatType(T2))), 4>
ApCross( const Vector<Expr1> &a, const coord<T2, 4>& b,const Vector<Expr2> &c)
{

	typedef coord<OutType(FloatType(Expr2),OutType(FloatType(Expr1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
ApCross(const coord<T1, 4>& a, const Vector<T2> &b,const coord<T3, 4> &c)
{

	typedef coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4> ApCross(const Vector<T1>& a, const coord<T2, 4> &b,const coord<T3, 4> &c)
{

	typedef coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4> OutT;
	return OutT(a(3)*(-(b(2)*c(1)) + b(1)*c(2)) + a(2)*(b(3)*c(1) - b(1)*c(3)) +
	    a(1)*(-(b(3)*c(2)) + b(2)*c(3)),
	   a(3)*(b(2)*c(0) - b(0)*c(2)) + a(2)*(-(b(3)*c(0)) + b(0)*c(3)) +
	    a(0)*(b(3)*c(2) - b(2)*c(3)),
	   a(3)*(-(b(1)*c(0)) + b(0)*c(1)) + a(1)*(b(3)*c(0) - b(0)*c(3)) +
	    a(0)*(-(b(3)*c(1)) + b(1)*c(3)),
	   a(2)*(b(1)*c(0) - b(0)*c(1)) + a(1)*(-(b(2)*c(0)) + b(0)*c(2)) +
    a(0)*(b(2)*c(1) - b(1)*c(2)));
}


template<class T1>
inline coord<FloatType(T1),1>
cross(const coord<T1, 1>& a)
{
	return ApCross(a);
}

template<class T1>
inline coord<FloatType(T1),2>
cross(const coord<T1, 2>& a)
{
	return ApCross(a);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)),3>
cross(const coord<T1, 3>& a, const coord<T2, 3> &b)
{
	return ApCross(a,b);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)),3>
cross(const VExpr<T1>& a, const coord<T2, 3> &b)
{
	return ApCross(a,b);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)),3>
cross(const coord<T1, 3>& a, const VExpr<T2> &b)
{
	return ApCross(a,b);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)),3>
cross(const Vector<T1>& a, const coord<T2, 3> &b)
{
	return ApCross(a,b);
}

template<class T1, class T2>
inline coord<OutType(FloatType(T1),FloatType(T2)),3>
cross(const coord<T1, 3>& a, const Vector<T2> &b)
{
	return ApCross(a,b);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const coord<T1, 4>& a, const coord<T2, 4> &b, const coord<T3, 4> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const coord<T1, 4>& a, const coord<T2, 4> &b, const VExpr<T3> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const coord<T1, 4>& a, const VExpr<T2> &b, const VExpr<T3> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const VExpr<T1>& a, const coord<T2, 4> &b, const VExpr<T3> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const VExpr<T1>& a, const coord<T2, 4> &b, const coord<T2, 4> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const VExpr<T1>& a, const VExpr<T2> &b, const coord<T3, 4> &c)
{
	return ApCross(a,b,c);
}


template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const coord<T1, 4>& a, const coord<T2, 4> &b, const Vector<T3> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const coord<T1, 4>& a, const Vector<T2> &b, const Vector<T3> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const Vector<T1>& a, const coord<T2, 4> &b, const Vector<T3> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const Vector<T1>& a, const coord<T2, 4> &b, const coord<T2, 4> &c)
{
	return ApCross(a,b,c);
}

template<class T1, class T2, class T3>
inline coord<OutType(FloatType(T3),OutType(FloatType(T1),FloatType(T2))), 4>
cross(const Vector<T1>& a, const Vector<T2> &b, const coord<T2, 4> &c)
{
	return ApCross(a,b,c);
}


/*********** Products ******************/

template<class T1>
inline FloatType(T1) ApProd(const coord<T1, 1>& a)
{
	return a[0];
}

template<class T1>
inline FloatType(T1) ApProd(const coord<T1, 2>& a)
{
	return a[0]*a[1];
}

template<class T1>
inline FloatType(T1) ApProd(const coord<T1, 3>& a)
{
	return a[0]*a[1]*a[2];
}

template<class T1>
inline FloatType(T1) ApProd(const coord<T1, 4>& a)
{
	return a[0]*a[1]*a[2]*a[3];
}

template<class T1, int N>
inline FloatType(T1) ApProd(const coord<T1, N>& a)
{
	FloatType(T1) ss=a[0];
	for(int i=1;i<N;i++) ss*=a[i];
	return ss;
}


template<class T1,int N>
inline FloatType(T1)
prod(const coord<T1, N>& a)
{
	return ApProd(a);
}

/*********** Sums ******************/

template<class T1>
inline SumType(T1) ApSum(const coord<T1, 1>& a)
{
	return a[0];
}

template<class T1>
inline SumType(T1) ApSum(const coord<T1, 2>& a)
{
	return a[0]+a[1];
}

template<class T1>
inline SumType(T1) ApSum(const coord<T1, 3>& a)
{
	return a[0]+a[1]+a[2];
}

template<class T1>
inline SumType(T1) ApSum(const coord<T1, 4>& a)
{
	return a[0]+a[1]+a[2]+a[3];
}

template<class T1, int N>
inline  SumType(T1) ApSum(const coord<T1, N>& a)
{
	SumType(T1) ss=a[0];
	for(int i=1;i<N;i++) ss*=a[i];
	return ss;
}


template<class T1,int N>
inline FloatType(T1)
sum(const coord<T1, N>& a)
{
	return ApSum(a);
}

/******************** ROTATIONS  EXTERNAL to CLASS *********/
template<class T, int N>
coord<T,N> Rotate3D(const coord<T,N> &in, double theta, double phi, double psi=0.0)
{
	T tx=in.x(), ty=in.y(), tz=in.z();
	coord<T,N> tmp;
	double cth=cos(theta), sth=sin(theta), cps=cos(psi), sps=sin(psi), cph=cos(phi), sph=sin(phi);
	tmp.put(0,tx*(cps*cph-cth*sph*sps)+
		ty*(cps*sph+cth*cph*sps)+
		tz*(sps*sth));
	tmp.put(1,tx*(-sps*cph-cth*sph*cps)+
		ty*(-sps*sph+cth*cph*cps)+
		tz*(cps*sth));
	tmp.put(2,tx*(sth*sph)-
		ty*(sth*cph)+
		tz*cth);
	return tmp;
}
template<class T, int N>
coord<T,N> Rotate3D(const coord<T,N> &in,const coord<double,N> &rot)
{
	T tx=in.x(), ty=in.y(), tz=in.z();
	coord<T,N> tmp;
	double cth=cos(rot.x()), sth=sin(rot.x()), cps=cos(rot.z()), sps=sin(rot.z()), cph=cos(rot.y()), sph=sin(rot.y());
	tmp.put(0,tx*(cps*cph-cth*sph*sps)+
		ty*(cps*sph+cth*cph*sps)+
		tz*(sps*sth));
	tmp.put(1,tx*(-sps*cph-cth*sph*cps)+
		ty*(-sps*sph+cth*cph*cps)+
		tz*(cps*sth));
	tmp.put(2,tx*(sth*sph)-
		ty*(sth*cph)+
		tz*cth);
	return tmp;
}

/**************** VECTOR COORD FUNCTIONS  ***********/

//****************************************************************************
#define VecFromCoord(op)                                            \
	template<class T> template<class T1, int N>                     \
	inline Vector<T>& Vector<T>::operator op (const coord<T1,N> &Texpr){     \
		for (int i=0; i < length_; ++i) {                              \
			(*this)[i] op Texpr;                                      \
		}                                                                   \
    	return *this;				\
	}

VecFromCoord(+=)
VecFromCoord(-=)
VecFromCoord(*=)
VecFromCoord(/=)


/*template<class T, int N>
inline coord<T,N> chop(const coord<T, N> &cc){
	return chop(cc, 1.e-12);
}

template<class T, int N>
inline coord<T,N> chop(const coord<T, N> &cc, double lim){
	return chop(cc, lim);
}*/

//***********************************************
//**************IO*******************
template<class T, int N>
std::ostream& operator<<(std::ostream& otr, const coord<T,N>& x){
	for (int i=0; i < N; ++i)	otr <<x[i]<<" ";
    return otr;
}

END_BL_NAMESPACE



#endif



