/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 04-11-02
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

/**************************************************
 matfft.h

 1) contains the code that does Fourier Transforms on Matrices

 it will use the FFTW library (http://www.fftw.org) if the lib
 has been configured for it..

 **************************************************/


#ifndef Mat_FFT_H_
#define Mat_FFT_H_ 1

//#include "math.h"
#include "container/Vector/Vector.h"
#include "container/rankType.h"
#include "container/matrix/_matrix.h"
#include "utils/utils.h"

//Need to see if FFTW is obtained...
#ifndef ON_WINDOWS
 #include "blochconfig.h"
#endif

BEGIN_BL_NAMESPACE


//The FFT SHIFT operation...useful for 'normal' viewing of
// FFT data going from -w...0...w rather then the
// algorithmic view of 0..w..-w

template<class nt1, class st1>
_matrix<nt1, st1>  rotateRight(const _matrix<nt1, st1> &in, int rSh, int cSh)
{
	//RunTimeAssert(sh<in.size() && sh>=0);
	_matrix<nt1, st1> outt(in.rows(), in.cols());
	int i,j;
	Vector<int> Rows(in.rows()), Cols(in.cols());

	Rows(Range(0, in.rows()-rSh-1))=Range(rSh, in.rows()-1);
	Rows(Range(in.rows()-rSh, in.rows()))=Range(0, rSh);

	Cols(Range(0, in.cols()-cSh-1))=Range(cSh, in.cols()-1);
	Cols(Range(in.cols()-cSh, in.cols()))=Range(0, cSh);

	for(i=0; i<outt.rows();++i){
		for(j=0; j<outt.cols();++j){
			outt(i,j)=in(Rows(i), Cols(j));
		}
	}
	return outt;
}

template<class nt1, class st1>
_matrix<nt1, st1>  fftshift(const _matrix<nt1, st1> &in)
{
	return rotateRight(in, ceil(double(in.rows())/2.0), ceil(double(in.cols())/2.0));
}

template<class nt1, class st1>
_matrix<nt1, st1>  ifftshift(const _matrix<nt1, st1> &in)
{
	return rotateRight(in, floor(double(in.rows())/2.0), floor(double(in.cols())/2.0));
}

/*** Matrix Expression ***/
template<class expr>
_matrix<
   typename _matrixExpr<expr>::numtype,
   typename _matrixExpr<expr>::structure>
fftshift(const _matrixExpr<expr> &in)
{
	typedef typename _matrixExpr<expr>::numtype nt1;
	typedef typename _matrixExpr<expr>::structure st1;
	return fftshift( _matrix<nt1, st1>(in));
}

template<class expr>
_matrix<
   typename _matrixExpr<expr>::numtype,
   typename _matrixExpr<expr>::structure>
ifftshift(const _matrixExpr<expr> &in)
{
	typedef typename _matrixExpr<expr>::numtype nt1;
	typedef typename _matrixExpr<expr>::structure st1;

	return ifftshift( _matrix<nt1, st1>(in));
}

//various 'auxilliarly names' for the FFT
template<class nt1, class st1>
_matrix<complex, FullMatrix> fft(const _matrix<nt1, st1> &in)
{	return FFT(in);	}

template<class nt1, class st1>
_matrix<complex, FullMatrix> fft2(const _matrix<nt1, st1> &in)
{	return FFT(in);	}

template<class nt1, class st1>
_matrix<complex, FullMatrix> FFT2(const _matrix<nt1, st1> &in)
{	return FFT(in);	}

template<class nt1, class st1>
_matrix<complex, FullMatrix> ifft(const _matrix<nt1, st1> &in)
{	return IFFT(in);	}

template<class nt1, class st1>
_matrix<complex, FullMatrix> ifft2(const _matrix<nt1, st1> &in)
{	return IFFT(in);	}

template<class nt1, class st1>
_matrix<complex, FullMatrix> IFFT2(const _matrix<nt1, st1> &in)
{	return IFFT(in);	}

template<class nt1, class st1>
_matrix<complex, FullMatrix> FFT(const _matrix<nt1, st1> &in, int dir)
{
	return (in==1)?FFT(in):IFFT(in);
}

template<class nt1, class st1>
_matrix<complex, FullMatrix> fft(const _matrix<nt1, st1> &in, int dir)
{
	return (in==1)?FFT(in):IFFT(in);
}

/*** Matrix Expression ***/
//various 'auxilliarly names' for the FFT Matrix Expressions
template<class expr1>
_matrix<complex, FullMatrix> fft(const _matrixExpr<expr1> &in)
{	return FFT(_matrix<complex, FullMatrix>(in));	}

template<class expr1>
_matrix<complex, FullMatrix> fft2(const _matrixExpr<expr1> &in)
{	return FFT(_matrix<complex, FullMatrix>(in));	}

template<class expr1>
_matrix<complex, FullMatrix> FFT2(const _matrixExpr<expr1> &in)
{	return FFT(_matrix<complex, FullMatrix>(in));	}

template<class expr1>
_matrix<complex, FullMatrix> ifft(const _matrixExpr<expr1> &in)
{	return IFFT(_matrix<complex, FullMatrix>(in));	}

template<class expr1>
_matrix<complex, FullMatrix> ifft2(const _matrixExpr<expr1> &in)
{	return IFFT(_matrix<complex, FullMatrix>(in));	}

template<class expr1>
_matrix<complex, FullMatrix> IFFT2(const _matrixExpr<expr1> &in)
{	return IFFT(_matrix<complex, FullMatrix>(in));	}

template<class expr1>
_matrix<complex, FullMatrix> FFT(const _matrixExpr<expr1> &in, int dir)
{
	return (in==1)?FFT(_matrix<complex, FullMatrix>(in)):IFFT(_matrix<complex, FullMatrix>(in));
}

template<class expr1>
_matrix<complex, FullMatrix> fft(const _matrixExpr<expr1> &in, int dir)
{
	return (in==1)?FFT(_matrix<complex, FullMatrix>(in)):IFFT(_matrix<complex, FullMatrix>(in));
}


/** Algorithm for NO FFTW library **/
#ifndef HAVE_FFTW

//bad algorithm, we need to go row by row, then col, by col
template<class nt1, class st1>
_matrix<complex, FullMatrix> FFT(const _matrix<nt1, st1> &in)
{
	int i;
	_matrix<complex, FullMatrix> out(in);
	for(i=0;i<in.rows();++i)
		out.putRow(i, FFT(out.row(i)));

	for(i=0;i<in.cols();++i)
		out.putCol(i, FFT(out.col(i)));

	return out;
}

//bad algorithm, we need to go row by row, then col, by col
template<class nt1, class st1>
_matrix<complex, FullMatrix> IFFT(const _matrix<nt1, st1> &in)
{
	int i;
	_matrix<complex, FullMatrix> out(in);
	for(i=0;i<in.rows();++i)
		out.putRow(i, IFFT(out.row(i)));

	for(i=0;i<in.cols();++i)
		out.putCol(i, IFFT(out.col(i)));

	return out;
}


/** USE FFTW If we HAVE it***/
#endif

#ifdef HAVE_FFTW


#ifdef HAVE_FFTW_II

#include "fftw.h"
//we need to create a workspace becuase the data format
// for the FFTW is slightly differenent then our format
// basically their complex is "re, im" and ours is "real, imag"

template<class nt1, class st1>
void mat_pack_fftw(const _matrix<nt1, st1> &in, fftw_complex *out, int r, int c)
{
	int i;
	const nt1 *dd=in.data();
	for(i=0;i<r*c;++i){
		c_re(out[i])=Re(dd[i]);
		c_im(out[i])=Im(dd[i]);
	}
}

template<class nt1, class st1>
void mat_unpack_fftw(_matrix<nt1, st1> &out, fftw_complex *in, int r, int c)
{
	int i;
	nt1 *dd=out.data();
	for(i=0;i<r*c;++i){
		Re(dd[i],c_re(in[i]));
		Im(dd[i],c_im(in[i]));
	}
}

template<class nt1, class st1>
_matrix<complex, FullMatrix> FFT(const _matrix<nt1, st1> &in)
{
	fftwnd_plan fplan;
	fplan=fftw2d_create_plan(in.rows(),in.cols(), FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_complex *fftw_dat,*fftw_outdat;

	fftw_dat= new fftw_complex[in.rows()*in.cols()];
	fftw_outdat=new fftw_complex[in.rows()*in.cols()];

//pack it into fftw format
	mat_pack_fftw(in, fftw_dat, in.rows(), in.cols());

//do the FFT
	fftwnd_one(fplan, fftw_dat, fftw_outdat);

//unpack the fftw format
	_matrix<complex, FullMatrix> out(in.rows(), in.cols());
	mat_unpack_fftw(out, fftw_outdat, in.rows(), in.cols());
	fftwnd_destroy_plan(fplan);
	delete [] fftw_dat;
	delete [] fftw_outdat;

	return out;
}


template<class nt1, class st1>
_matrix<complex, FullMatrix> IFFT(const _matrix<nt1, st1> &in)
{
	fftwnd_plan fplan;
	fplan=fftw2d_create_plan(in.rows(),in.cols(), FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_complex *fftw_dat,*fftw_outdat;

	fftw_dat= new fftw_complex[in.rows()*in.cols()];
	fftw_outdat=new fftw_complex[in.rows()*in.cols()];

//pack it into fftw format
	mat_pack_fftw(in, fftw_dat, in.rows(), in.cols());

//do the iFFT
	fftwnd_one(fplan, fftw_dat, fftw_outdat);

//unpack the fftw format
	_matrix<complex, FullMatrix> out(in.rows(), in.cols());
	mat_unpack_fftw(out, fftw_outdat, in.rows(), in.cols());
	fftwnd_destroy_plan(fplan);
	delete [] fftw_dat;
	delete [] fftw_outdat;

//return normalized fft
	return out/(in.rows()*in.cols());
}


#endif

#ifdef HAVE_FFTW_III
#include <fftw3.h>

//we need to create a workspace becuase the data format
// for the FFTW is slightly differenent then our format
// basically their complex is "re, im" and ours is "real, imag"

template<class nt1, class st1>
void mat_pack_fftw(const _matrix<nt1, st1> &in, fftw_complex *out, int r, int c)
{
	int i;
	const nt1 *dd=in.data();
	for(i=0;i<r*c;++i){
		out[i][0]=Re(dd[i]);
		out[i][1]=Im(dd[i]);
	}
}

template<class nt1, class st1>
void mat_unpack_fftw(_matrix<nt1, st1> &out, fftw_complex *in, int r, int c)
{
	int i;
	nt1 *dd=out.data();
	for(i=0;i<r*c;++i){
		Re(dd[i],in[i][0]);
		Im(dd[i],in[i][1]);
	}
}

template<class nt1, class st1>
_matrix<complex, FullMatrix> FFT(const _matrix<nt1, st1> &in)
{
	fftw_plan fplan;
	fftw_complex *fftw_dat,*fftw_outdat;

	fftw_dat= new fftw_complex[in.rows()*in.cols()];
	fftw_outdat=new fftw_complex[in.rows()*in.cols()];

	//pack it into fftw format
	mat_pack_fftw(in, fftw_dat, in.rows(), in.cols());

	//do the FFT
	fplan = fftw_plan_dft_2d( in.rows(),  in.cols(),
						  fftw_dat, fftw_outdat,
						  FFTW_FORWARD,FFTW_ESTIMATE);
						  
	fftw_execute(fplan);

	//unpack the fftw format
	_matrix<complex, FullMatrix> out(in.rows(), in.cols());
	mat_unpack_fftw(out, fftw_outdat, in.rows(), in.cols());
	fftw_destroy_plan(fplan);
	delete [] fftw_dat;
	delete [] fftw_outdat;

	return out;
}


template<class nt1, class st1>
_matrix<complex, FullMatrix> IFFT(const _matrix<nt1, st1> &in)
{
	fftw_plan fplan;
	
	fftw_complex *fftw_dat,*fftw_outdat;

	fftw_dat= new fftw_complex[in.rows()*in.cols()];
	fftw_outdat=new fftw_complex[in.rows()*in.cols()];

	//pack it into fftw format
	mat_pack_fftw(in, fftw_dat, in.rows(), in.cols());

	fplan = fftw_plan_dft_2d( in.rows(),  in.cols(),
						   fftw_dat, fftw_outdat,
						   FFTW_BACKWARDS,FFTW_ESTIMATE);

	fftw_execute(fplan);
	
	//unpack the fftw format
	_matrix<complex, FullMatrix> out(in.rows(), in.cols());
	mat_unpack_fftw(out, fftw_outdat, in.rows(), in.cols());
	fftw_destroy_plan(fplan);
	delete [] fftw_dat;
	delete [] fftw_outdat;

	//return normalized fft
	return out/(in.rows()*in.cols());
}

#endif
#endif

END_BL_NAMESPACE


#endif

