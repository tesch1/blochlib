/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 04.12.02
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

/* Vectoraux.h -- contains basic FFT on a vector */

#ifndef _Vectoraux_h_
#define _Vectoraux_h_ 1


//Functions specific to a particular Vector Type (i.e. initiats an instance of Vector<T>
#include "container/Vector/Vector.h"

//Need to see if FFTW is obtained...
#ifndef ON_WINDOWS
 #include "blochconfig.h"
#endif


BEGIN_BL_NAMESPACE


//The FFT SHIFT operation...useful for 'normal' viewing of
// FFT data going from -w...0...w rather then the
// algorithmic view of 0..w..-w

template<class T>
Vector<T> rotateRight(const Vector<T> &in, int sh)
{
	//RunTimeAssert(sh<in.size() && sh>=0);
	Vector<T> out(in.size());
	int i,j;
	for(i=sh, j=0; i<in.size();++i, ++j)	out[j]=in[i];
	for(i=0, j=sh;i<sh;++i, ++j) out[j]=in[i];
	return out;
}

template<class T>
Vector<T> fftshift(const Vector<T> &in)
{
    return rotateRight(in, int(std::ceil(double(in.size())/2.0)));
}

template<class T>
Vector<T> ifftshift(const Vector<T> &in)
{
    return rotateRight(in, int(std::floor(double(in.size())/2.0)));
}

//wrappers for 'lowercase' functioncalling
template<class T>
Vector<complex> ifft(const Vector<T> &in)
{	return IFFT(in);	}

template<class T>
Vector<complex> fft(const Vector<T> &in)
{	return FFT(in);	}

template<class T>
Vector<complex> fft(const Vector<T> &in, int ff)
{	return FFT(in, ff);	}

/** Vector Expression **/
//wrappers for 'lowercase' functioncalling
template<class T>
Vector<complex> ifft(const VExpr<T> &in)
{	return IFFT(Vector<typename VExpr<T>::numtype>(in));	}

template<class T>
Vector<complex> fft(const VExpr<T> &in)
{	return FFT(Vector<typename VExpr<T>::numtype>(in));	}

template<class T>
Vector<complex> fft(const VExpr<T> &in, int ff)
{	return FFT(Vector<typename VExpr<T>::numtype>(in), ff);	}

/** Vector Expression **/
//wrappers for 'Uppercase' functioncalling
template<class T>
Vector<complex> IFFT(const VExpr<T> &in)
{	return IFFT(Vector<typename VExpr<T>::numtype>(in));	}

template<class T>
Vector<complex> FFT(const VExpr<T> &in)
{	return FFT(Vector<typename VExpr<T>::numtype>(in));	}

template<class T>
Vector<complex> FFT(const VExpr<T> &in, int ff)
{	return FFT(Vector<typename VExpr<T>::numtype>(in), ff);	}


//Algorithm for NO FFTW library
#ifndef HAVE_FFTW
template<class T1>
Vector<complex> FFT_(const Vector<T1> &in, int isign)//FFT OUTPUT WILL BE COMPEX!!!
{
	Vector<complex> dat_;
	dat_=in;
	int size=in.size();
	if(size==0) return dat_;
	complex w, wp, temp;
	double theta, sto2;
	int i,j,m,mmax,istep,ii;

	#ifndef PI
	#define PI 3.14159265358979323846
	#endif

//	if(isign==-1){
//		for(i=0,j=size/2; j<size; i++,j++)	swap_(dat_[i], dat_[j]);
//	}
	for(j=0,i=0; i<size; i++){
		if(j>i )	swap_(dat_[i],dat_[j]);
		m = size/2;
		while((m>1) && (j>=m)){
			j -=m;
			m /= 2;
		}
		j+=m;
	}
	mmax=1;
	while (size>mmax){
		istep = 2*mmax;
		theta = PI/(isign*mmax);
		sto2 = sin(theta*.5);
		wp = complex(-2.*sto2*sto2, sin(theta));
		w=1;
		for(ii=0; ii<mmax; ii++){
			for(i=ii; i<size; i+=istep){
				j = i+mmax;			// Danielson Lanczos
				temp = w*dat_[j];
				dat_[j] = dat_[i] - temp;
				dat_[i] += temp;
			}
			w+=w*wp;
		}
		mmax*=2;
	}
	return dat_;
}

template< class T >
Vector<complex> IFFT(const Vector<T> &in)
{
	return FFT_(in, -1);
}

template< class T >
Vector<complex> FFT(const Vector<T> &in)
{
	return FFT_(in, 1);
}

#endif

#ifdef HAVE_FFTW



#ifdef HAVE_FFTW_II
//Algorithm for FFTW VERSION 2 library

#include "fftw.h"
//we need to create a workspace becuase the data format
// for the FFTW is slightly differenent then our format
// basically their complex is "re, im" and ours is "real, imag"

template<class T>
void vec_pack_fftw(T *in, fftw_complex *out, int len)
{
	for(int i=0;i<len;++i){
		c_re(out[i])=Re(in[i]);
		c_im(out[i])=Im(in[i]);
	}
}

template<class T>
void vec_unpack_fftw(T *out, fftw_complex *in, int len)
{
	for(int i=0;i<len;++i){
		Re(out[i],c_re(in[i]));
		Im(out[i],c_im(in[i]));
	}
}


template<class T>
Vector<complex> FFT(const Vector<T> &in)
{
	fftw_plan fplan;
	fplan=fftw_create_plan(in.size(), FFTW_FORWARD,FFTW_ESTIMATE);

	Vector<complex> out(in.size());

	fftw_complex *fftw_dat, *fftw_outdat;
	fftw_dat= new fftw_complex[in.size()];
	fftw_outdat=new fftw_complex[in.size()];

	vec_pack_fftw(in.data(), fftw_dat,in.size());

	fftw_one(fplan, fftw_dat, fftw_outdat);
	vec_unpack_fftw(out.data(), fftw_outdat, in.size());
	fftw_destroy_plan(fplan);
	delete [] fftw_dat;
	delete [] fftw_outdat;

	return out;
}

template<class T>
Vector<complex> IFFT(const Vector<T> &in)
{
	fftw_plan fplan;
	fplan=fftw_create_plan(in.size(), FFTW_BACKWARD,FFTW_ESTIMATE);

	Vector<complex> out(in.size());

	fftw_complex *fftw_dat, *fftw_outdat;
	fftw_dat= new fftw_complex[in.size()];
	fftw_outdat=new fftw_complex[in.size()];

	vec_pack_fftw(in.data(), fftw_dat,in.size());

	fftw_one(fplan, fftw_dat, fftw_outdat);
	vec_unpack_fftw(out.data(), fftw_outdat, in.size());
	fftw_destroy_plan(fplan);
	delete [] fftw_dat;
	delete [] fftw_outdat;
	return out/in.size();
}

#endif

/** VERSION 3 OF FFTW ***/

#ifdef HAVE_FFTW_III
#include <fftw3.h>
//we need to create a workspace becuase the data format
// for the FFTW is slightly differenent then our format
// basically their complex is "re, im" and ours is "real, imag"

template<class T>
void vec_pack_fftw(T *in, fftw_complex *out, int len)
{
	for(int i=0;i<len;++i){
		out[i][0]=Re(in[i]);
		out[i][1]=Im(in[i]);
	}
}



template<class T>
void vec_unpack_fftw(T *out, fftw_complex *in, int len)
{
	for(int i=0;i<len;++i){
		Re(out[i],(in[i][0]));
		Im(out[i],(in[i][1]));
	}
}

template<class T>
Vector<complex> FFT(const Vector<T> &in)
{
	fftw_plan fplan;
	Vector<complex> out(in.size());

	fftw_complex *fftw_dat, *fftw_outdat;
	fftw_dat= new fftw_complex[in.size()];
	fftw_outdat=  new fftw_complex[in.size()];

	vec_pack_fftw(in.data(), fftw_dat,in.size());

	fplan = fftw_plan_dft_1d(in.size(), fftw_dat, fftw_outdat, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fplan);
	vec_unpack_fftw(out.data(), fftw_outdat, in.size());
	
	fftw_destroy_plan(fplan);

	delete [] fftw_dat;
	delete [] fftw_outdat;
	
	return out;
}

template<class T>
Vector<complex> IFFT(const Vector<T> &in)
{
	fftw_plan fplan;
	Vector<complex> out(in.size());

	fftw_complex *fftw_dat, *fftw_outdat;
	fftw_dat= new fftw_complex[in.size()];
	fftw_outdat=  new fftw_complex[in.size()];

	vec_pack_fftw(in.data(), fftw_dat,in.size());
	
	fplan = fftw_plan_dft_1d(in.size(), fftw_dat, fftw_outdat, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(fplan);
	
	vec_unpack_fftw(out.data(), fftw_outdat, in.size());
	fftw_destroy_plan(fplan);

	delete [] fftw_dat;
	delete [] fftw_outdat;
	return out/in.size();
}

#endif //end HAVE_FFTW_III


template<class T>
Vector<complex> FFT(const Vector<T> &in, int ff)
{
	return (ff==1)?FFT(in):IFFT(in);
}


#endif


END_BL_NAMESPACE





#endif
