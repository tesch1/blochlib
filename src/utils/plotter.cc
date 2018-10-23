/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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
	powell.h-->METHODS for fitting of 1-D spectra..several basic minimazation
	and bracketing algorithms are used here

	main methods are found in 'powellmethods.h'


*/

#ifndef _plotter_cc_
#define _plotter_cc_

#include "container/Vector/Vector.h"
#include "utils/plotter.h"
#include "utils/utils.h"
#include <string>

#ifndef ON_WINDOWS
 #include "blochconfig.h"
#endif


BEGIN_BL_NAMESPACE



std::string plotterFID(Vector<Complex<double> > &fid, std::string datasave, double dt, double appo, bool dccorrect, bool zfill, bool toplot){
	int npts=fid.size();
	double dw=1./dt/npts;
	double min_w=-1./dt/2.;
	Vector<complex> tmfid=fid;
	double tmt=0;

#ifndef HAVE_FFTW

	//check the power of 2
	double csize=tmfid.size();
	int mypow; //hold the 2^mypow
	double f=frexp(csize,&mypow);
	if(f==0.5) mypow=mypow-1; //already a power of 2

	int newsize=int(pow(2.0, double(mypow)));
	if(npts!=newsize){
		std::cerr<<"Warning: VEctor size not a power of 2...currently, FFTs "<<std::endl;
		std::cerr<<" can only be done with powers of 2..." <<std::endl;
		std::cerr<<" zero filling to the next power of 2"<<std::endl;
		tmfid.resizeAndPreserve(newsize, ZeroType<complex>::zero());
		npts=tmfid.size();
		dw=1./dt/npts;
		min_w=-1./dt/2.;
	}
#endif

	if(dccorrect){
		tmfid-=sum(tmfid)/tmfid.size();
		std::cout<<sum(tmfid)/tmfid.size()<<std::endl;
	}
	if(appo!=0){
		for(int i=0;i<fid.size();i++){
			tmfid(i)*=exp(-abs(appo)*tmt);
			tmt+=dt;
		}
	}
	int tmpows=0;
	if(zfill){

#ifndef HAVE_FFTW
		tmpows=int(pow(2.0, double(mypow+1)));
#else
		tmpows=2*tmfid.size();
#endif
		tmfid.resizeAndPreserve(tmpows, ZeroType<complex>::zero());
		dw=1./dt/tmfid.size();
	}


	Vector<complex> tmfft=FFT(tmfid);
    tmfft=fftshift(tmfft);
    std::string dat_file=datasave;
	std::ifstream testdat(dat_file.c_str());
	int i=0, failed=false;
	while(!testdat.fail()){
		testdat.close();
		dat_file=datasave+itost(i);
		testdat.open(dat_file.c_str());
		i++;
		failed=true;
	}
	if(failed){
		std::cout<<std::endl<<"data file "<<datasave<<std::endl;
		std::cout<<" already exsists...changing name to "<<dat_file<<std::endl;
	}

	std::ofstream oo(dat_file.c_str());
	if(oo.fail()){
		std::cerr<<std::endl<<"Error: write data..."<<std::endl;
		std::cerr<<" Could not open the file to write data......"<<std::endl;
		std::cerr<<" you will have to figure out the problem "<<std::endl;
		std::cerr<<" (i.e. directory permissions, etc) and run the simulation again"<<std::endl;
		std::cerr<<" good bye..."<<std::endl;
		return dat_file;
	}
	tmt=0;
	for(i=0; i<tmfid.size();i++){
		oo<<tmt<<" "<<min_w<<" "<<Re(tmfft(i))<<" "<<Im(tmfft(i))<<" "
		<<Re(tmfid(i))<<" "<<Im(tmfid(i))<<" "<<norm(tmfft(i))<<std::endl;
		min_w+=dw;
		tmt+=dt;
	}
	oo.close();
	return dat_file;

}



std::string plotterFID(Vector<Complex<float> > &fid, std::string datasave, double dt, double appo, bool dccorrect, bool zfill, bool toplot){
	int npts=fid.size();
	float dw=1./dt/npts;
	float min_w=-1./dt/2.;
	Vector<Complex<float> > tmfid=fid;
	float tmt=0;

#ifndef HAVE_FFTW

	//check the power of 2
	float csize=tmfid.size();
	int mypow; //hold the 2^mypow
	double f=frexp(csize,&mypow);
	if(f==0.5) mypow=mypow-1; //already a power of 2

	int newsize=int(pow(2.0, double(mypow)));
	if(npts!=newsize){
		std::cerr<<"Warning: VEctor size not a power of 2...currently, FFTs "<<std::endl;
		std::cerr<<" can only be done with powers of 2..." <<std::endl;
		std::cerr<<" zero filling to the next power of 2"<<std::endl;
		tmfid.resizeAndPreserve(newsize, ZeroType<complex>::zero());
		npts=tmfid.size();
		dw=1./dt/npts;
		min_w=-1./dt/2.;
	}
#endif

	if(dccorrect){
		tmfid-=sum(tmfid)/tmfid.size();
		std::cout<<sum(tmfid)/tmfid.size()<<std::endl;
	}
	if(appo!=0){
		for(int i=0;i<fid.size();i++){
			tmfid(i)*=exp(-abs(appo)*tmt);
			tmt+=dt;
		}
	}
	int tmpows=0;
	if(zfill){

#ifndef HAVE_FFTW
		tmpows=int(pow(2.0, double(mypow+1)));
#else
		tmpows=2*tmfid.size();
#endif
		tmfid.resizeAndPreserve(tmpows, ZeroType<Complex<float> >::zero());
		dw=1./dt/tmfid.size();
	}


	Vector<Complex<float> > tmfft=FFT(tmfid);
    tmfft=fftshift(tmfft);
    std::string dat_file=datasave;
	std::ifstream testdat(dat_file.c_str());
	int i=0, failed=false;
	while(!testdat.fail()){
		testdat.close();
		dat_file=datasave+itost(i);
		testdat.open(dat_file.c_str());
		i++;
		failed=true;
	}
	if(failed){
		std::cout<<std::endl<<"data file "<<datasave<<std::endl;
		std::cout<<" already exsists...changing name to "<<dat_file<<std::endl;
	}

	std::ofstream oo(dat_file.c_str());
	if(oo.fail()){
		std::cerr<<std::endl<<"Error: write data..."<<std::endl;
		std::cerr<<" Could not open the file to write data......"<<std::endl;
		std::cerr<<" you will have to figure out the problem "<<std::endl;
		std::cerr<<" (i.e. directory permissions, etc) and run the simulation again"<<std::endl;
		std::cerr<<" good bye..."<<std::endl;
		return dat_file;
	}
	tmt=0;
	for(i=0; i<tmfid.size();i++){
		oo<<tmt<<" "<<min_w<<" "<<Re(tmfft(i))<<" "<<Im(tmfft(i))<<" "
		<<Re(tmfid(i))<<" "<<Im(tmfid(i))<<" "<<norm(tmfft(i))<<std::endl;
		min_w+=dw;
		tmt+=dt;
	}
	oo.close();
	return dat_file;

}
//writes a gnuplot script to plot the data from the plotter function above
// if 'toplot' is true it makes a system call to run the script file
void WriteGnuplotFID(std::string datasave, bool toplot)
{
	std::ofstream gnplot("datplot");
		std::string mess="set data style lines\n";
		mess+="set origin 0,0\n";
		mess+="set size 1,1\n";
		mess+="set title \"Bloch Eq Sim\"\n";
		mess+="set nokey\n";
	mess+="\n";
		mess+="set multiplot\n";
		mess+="set border 31\n";
	mess+="\n";
		mess+="set ylabel \"\"\n";
		mess+="set xlabel \"t\", 1\n";
	mess+="\n";
		mess+="set origin 0,0\n";
		mess+="set size .3,.48\n";
		mess+="set title \"imag FID\" ,-1\n";
		mess+="plot \""+datasave+"\" u 1:6\n";
	mess+="\n";
		mess+="set origin 0,.52\n";
		mess+="set title \"real FID\" ,-1\n";
		mess+="plot \""+datasave+"\" u 1:5\n";
	mess+="\n";
		mess+="set xlabel \"Hz\"	,1\n";
	mess+="\n";
		mess+="set origin .3,.52\n";
		mess+="set title \"real fft\" ,-1\n";
		mess+="plot \""+datasave+"\" u 2:3\n";

		mess+="set origin .3,.0\n";
		mess+="set title \"imag fft\" ,-1\n";
		mess+="plot \""+datasave+"\" u 2:4\n";

		mess+="set size .3,1\n";
		mess+="set origin .6,0\n";
		mess+="set title \"power fft\" ,-1\n";
		mess+="plot \""+datasave+"\" u 2:7\n";

		mess+="set nomultiplot\n";
		mess+="set term postscript\n";
		mess+="set output 'printme.ps'\n";
		mess+="set multiplot\n";
		mess+="set origin 0,0\n";
		mess+="set size 1,1\n";
		mess+="set title \"Bloch Eq Sim\"\n";
		mess+="set nokey\n";
	mess+="\n";
		mess+="set multiplot\n";
		mess+="set border 31\n";
	mess+="\n";
		mess+="set ylabel \"\"\n";
		mess+="set xlabel \"t\", 1\n";
	mess+="\n";
		mess+="set origin 0,0\n";
		mess+="set size .3,.48\n";
		mess+="set title \"imag FID\" ,-1\n";
		mess+="plot \""+datasave+"\" u 1:6\n";
	mess+="\n";
		mess+="set origin 0,.52\n";
		mess+="set title \"real FID\" ,-1\n";
		mess+="plot \""+datasave+"\" u 1:5\n";
	mess+="\n";
		mess+="set xlabel \"Hz\"	,1\n";
	mess+="\n";
		mess+="set origin .3,.52\n";
		mess+="set title \"real fft\" ,-1\n";
		mess+="plot \""+datasave+"\" u 2:3\n";

		mess+="set origin .3,.0\n";
		mess+="set title \"imag fft\" ,-1\n";
		mess+="plot \""+datasave+"\" u 2:4\n";

		mess+="set size .3,1\n";
		mess+="set origin .6,0\n";
		mess+="set title \"power fft\" ,-1\n";
		mess+="plot \""+datasave+"\" u 2:7\n";

		mess+="set nomultiplot\n";
		mess+="set border 31\n";
		mess+="set xlabel \"\"\n";
		mess+="set zlabel \"\"\n";
		mess+="set ylabel \"\"\n";
		mess+="set title \"\"\n";
		mess+="set size 1,1\n";
		mess+="set origin 0,0\n";
		mess+="set nolabel\n";
	#ifdef __CYGWIN__
		mess+="set term windows\n";
	#else
		mess+="set term x11\n";
	#endif

	if(gnplot.fail()){
		std::cerr<<std::endl<<"Error: write plotter file..."<<std::endl;
		std::cerr<<" Could not open the file to write plotting script......"<<std::endl;
		std::cerr<<" you will have to figure out the problem "<<std::endl;
		std::cerr<<" (i.e. directory permissions, etc) "<<std::endl<<std::endl;
		return;
	}
	gnplot<<mess<<std::endl;
	if(toplot){
		std::string sysca="gnuplot datplot";
		system(sysca.c_str());
	}
}


//writes a gnuplot scripto to plot the data from
// if 'toplot' is true it makes a system call to run the script file
// the 'cols' is how many columns in the saved file are present....
void WriteGnuplotLyp(std::string datasave, int cols, bool toplot)
{
	std::ofstream gnplot("lyplot");
	std::string mess="set data style lines\n";
	mess+="set title \"Blochlib Lyapunov plot\"\n";
	mess+="set nokey\n";
	if(cols>0) mess+="plot ";
	for(int i=2;i<cols;i++)
	{
		if(i==2) mess += "'"+datasave+"' u 1:"+itost(i)+" ";
		else mess += ",'"+datasave+"' u 1:"+itost(i)+" ";
	}
	if(gnplot.fail()){
		std::cerr<<std::endl<<"Error: write plotter file..."<<std::endl;
		std::cerr<<" Could not open the file to write plotting script......"<<std::endl;
		std::cerr<<" you will have to figure out the problem "<<std::endl;
		std::cerr<<" (i.e. directory permissions, etc) "<<std::endl<<std::endl;
		return;
	}
	gnplot<<mess<<std::endl;
	if(toplot)
	{
		std::string sysca="gnuplot lypplot";
		system(sysca.c_str());
	}
}


END_BL_NAMESPACE



#endif


