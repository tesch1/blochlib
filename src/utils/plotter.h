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
	plotter.h-->a single easy to use plotter function
	that writes the data out in a 6 column format

	  frequncey time RealFID ImagFID RealFFT ImagFFT MagnitudeFFT

	also writes a gnuplot (http://gnuplot.org) script to plot
	the data nicely...

*/

#ifndef _plotter_h_
#define _plotter_h_

#include "container/Vector/Vector.h"
#include <string>

BEGIN_BL_NAMESPACE


//dumps data in 'fid' to a file called 'datasave' and performs various
//dataprocs if desired
// returns the 'true' file name...it does NOT over write any data file
// that happens t have the same file name

std::string plotterFID(Vector<Complex<double> > &fid, std::string datasave="fid", double dt=1, double appo=0, bool dccorrect=0, bool zfill=0, bool toplot=0);
std::string plotterFID(Vector<Complex<float> > &fid, std::string datasave="fid", double dt=1, double appo=0, bool dccorrect=0, bool zfill=0, bool toplot=0);
//std::string plotter(Vector<complex> &fid, std::string datasave="fid", double dt=1, double appo=0, bool dccorrect=0, bool zfill=0, bool toplot=0)
//{	return std::plotterFID(fid,datasave, dt, appo, dccorrect, zfill, toplot);	}

//writes a gnuplot script to plot the data from the plotter function above
// if 'toplot' is true it makes a system call to run the script file
void WriteGnuplotFID(std::string datasave, bool toplot=false);


//writes a gnuplot scripto to plot the data from
// if 'toplot' is true it makes a system call to run the script file
// the 'cols' is how many columns in the saved file are present....
void WriteGnuplotLyp(std::string datasave, int cols, bool toplot=false);

END_BL_NAMESPACE


#endif


