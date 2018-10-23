/* Isotope.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-4-01
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
	Isotope.cc--> atomic isotopes data storage container....

I'd like to thank to the makers of Gamma

  S.A. Smith, T.O. Levante, B.H. Meier and R.R. Ernst
  Computer Simulations in Magnetic Resonance:
  An Object-Oriented Programming Approach
  J. Magn. Reson., 106a, 75-105, (1994S, Series A)

http://gamma.magnet.fsu.edu/

for the hard data

*/

#ifndef _Isotope_cc_
#define _Isotope_cc_ 1

#include "bloch/Isotope.h"
#include "utils/utils.h"
#include "utils/constants.h"


BEGIN_BL_NAMESPACE


const int NIso=131;

const int IsoSpins[NIso] = { 2, 3, 2, 2, 1, 1, 3, 4, 5, 4, 7, 4, 1, 2, 3, 2, 1,
				 6, 1, 2, 1, 4, 1, 4, 6, 6, 2, 2, 4, 4, 4, 4, 4, 8,
				 8, 6, 8, 13, 8, 4, 6, 2, 8, 4, 4, 4, 6, 4, 4, 10, 4,
				 2, 4, 4, 6, 4, 10, 2, 6, 10, 6, 6, 10, 6, 6, 2, 6, 2,
				 2, 2, 2, 10, 10, 2, 2, 6, 8, 2, 2, 6, 2, 4, 8, 4, 4,
				 11, 8, 6, 8, 8, 8, 8, 6, 6, 4, 4, 4, 6, 6, 8, 8, 2,
				 2, 6, 8, 15, 8, 10, 8, 2, 6, 6, 2, 4, 4, 4, 2, 4, 2,
				 4, 2, 2, 2, 10, 10, 4, 6, 4, 4, 8, 2 };

/*          An array of known isotope atomic numbers (1=H, 6=C, ...)        */

const int IsoNumbers[NIso] = { 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 5, 6, 6,
				   7, 7, 8, 8, 8, 9, 10, 10, 10, 11, 12, 13, 14, 15,
				   16, 17, 17, 19, 19, 20, 21, 22, 22, 23, 23, 24, 25, 26,
				   27, 28, 29, 29, 30, 31, 31, 31, 32, 33, 34, 35, 37, 37,
				   38, 39, 40, 41, 42, 42, 43, 44, 44, 45, 46, 47, 47, 48,
				   48, 49, 49, 50, 50, 51, 51, 52, 52, 53, 54, 54, 55, 56,
				   56, 57, 57, 59, 60, 60, 62, 62, 63, 63, 64, 64, 65, 66,
				   66, 67, 68, 69, 70, 70, 71, 71, 72, 72, 73, 74, 75, 75,
				   76, 76, 77, 77, 78, 79, 80, 80, 81, 81, 82, 83, 85, 89,
				   90, 91, 91, 92, 0 };

/*     An array of known isotope atomic masses  (1=1H, 2=2H, 13=13C, ...)   */

const int IsoMasses[NIso] = { 1, 2, 3, 3, 4, 6, 6, 7, 8, 9, 10,
				  11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
				  21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41,
				  43, 45, 47, 49, 50, 51, 53, 55, 57, 59, 61,
				  63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 85,
				  87, 87, 89, 91, 93, 95, 97, 99, 99, 101, 103,
				  105, 107, 109, 111, 113, 113, 115, 117, 119, 121, 123,
				  123, 125, 127, 129, 131, 133, 135, 137, 138, 139, 141,
				  143, 145, 147, 149, 151, 153, 155, 157, 159, 161, 163,
				  165, 167, 169, 171, 173, 175, 176, 177, 179, 181, 183,
				  185, 187, 187, 189, 191, 193, 195, 197, 199, 201, 203,
				  205, 207, 209, 211, 227, 229, 231, 233, 235, 0 };

/*   An array of known isotope atomic weights (1H = 1.008665 g/mole, ...)   */

const double IsoWeights[NIso] = { 1.008665, 2.0140, 3.01605, 3.01603, 4.00260,
					  6.01888, 6.01512, 7.01600, -1.1e6, 9.01218,
					  10.0129, 11.00931, 12.00000, 13.00335, 14.00307,
					  15.00011, 15.99491, -1.1e6, -1.1e6, 18.99840,
					  19.99244, 20.99395, 21.99138, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6,
					  0.000001 };

/*        An array of known isotope receptivities (based on 13C = 1)        */

const double IsoRecepts[NIso] = { 5.68E+3, 8.21E-3, -1.1e6, 3.26E-3, -1.1e6,
					  -1.1e6, 3.58, 1.54E+3, -1.1e6, 7.88E+1,
					  2.21E+1, 7.54E+2, -1.1e6, 1.00, 5.69,
					  2.19E-2, -1.1e6, 6.11E-2, -1.1e6, 4.73E+3,
					  -1.1e6, 3.59E-2, -1.1e6, 5.25E+2, 1.54,
					  1.17E+3, 2.09, 3.77E+2, 9.73E-2, 2.02E+1,
					  3.77, 2.69, 3.28E-2, 5.27E-2, 1.71E+3,
					  8.64E-1, 1.18, 7.55E-1, 2.16E+3, 4.90E-1,
					  9.94E+2, 4.19E-3, 1.57E+3, 2.41E-1, 3.65E+2,
					  2.01E+2, 6.65E-1, 2.37E+2, 3.19E+2, 6.17E-1,
					  1.43E+2, 2.98, 2.26E+2, 2.77E+2, 4.34E+1,
					  2.77E+2, 1.07, 6.68E-1, 6.04, 2.74E+3,
					  2.88, 1.84, 2.13E+3, 8.30E-1, 1.56,
					  1.77E-1, 1.41, 1.95E-1, 2.76E-1, 6.93,
					  7.6, 8.38E+1, 1.89E+3, 1.95E+1, 2.52E+1,
					  5.20E+2, 1.11E+2, 8.90E-1, 1.25E+1, 5.30E+2,
					  3.18E+1, 3.31, 2.69E+2, 1.83, 4.41,
					  4.30E-1, 3.36E+2, 1.65E+3, 2.34, 3.70E-1,
					  1.26, 0.59, 4.83E+2, 4.53E+1, 2.33E-1,
					  4.84E-1, 3.31E+2, 4.47E-1, 1.58, 1.02E+3,
					  6.59E-1, 3.21, 4.44, 1.22, 1.72E+2,
					  5.47, 6.70E-1, 1.69E-1, 2.04E+2, 5.89E-2,
					  2.80E+2, 4.90E+2, 1.14E-3, 2.13, 5.36E-2,
					  1.16E-1, 1.91E+1, 1.43E-1, 5.42, 1.08,
					  2.89E+2, 7.69E+2, 1.18E+1, 7.77E+2, -1.1e6,
					  -1.1e6, -1.1e6, -1.1e6, -1.1e6, 4.95E-3,
					  -1.1e6 };

/*  An array of known isotope relative frequencies (from 1H @ 400.130 MHz)  */

const double IsoRelFreqs[NIso] = { 400.130, 61.424, 426.791, 304.811, -1.1e6, -1.1e6,
					   58.883, 155.503, -1.1e6, 56.226, 42.990, 128.330,
					   -1.1e6, 100.613, 28.913, -40.561, -1.1e6, 54.242,
					   -1.1e6, 376.498, -1.1e6, 31.586, -1.1e6, 105.842,
					   24.496, 104.262, 79.494, 161.977, 30.714, 39.205,
					   32.635, 18.670, 10.247, 26.925, 97.200, 22.563,
					   22.557, 39.893, 105.246, 22.615, 99.012, 12.956,
					   94.939, 35.756, 106.146, 113.705, 25.036, 96.035,
					   122.028, 13.957, 68.514, 76.313, 100.249, 108.063,
					   38.633, 130.927, 17.342, 19.606, 37.196, 97.936,
					   26.076, 26.625, 90.061, 18.462, 20.691, 12.744,
					   18.310, 16.197, 18.622, 84.832, 88.741, 87.492,
					   87.680, 142.574, 149.212, 95.755, 51.853, 104.714,
					   126.241, 80.062, 110.668, 32.811, 52.481, 39.749,
					   44.466, 52.793, 56.522, 117.202, 21.755, 13.384,
					   16.517, 13.160, 99.236, 43.818, 15.281, 19.102,
					   90.741, 13.180, 18.338, 82.079, 11.564, 33.095,
					   70.475, 19.414, 45.643, 31.722, 12.484, 7.478,
					   47.976, 16.669, 90.129, 91.038, 9.131, 31.070,
					   6.874, 7.486, 85.876, 6.918, 71.667, 26.457,
					   228.970, 231.219, 83.459, 64.297, -1.1e6, -1.1e6,
					   -1.1e6, -1.1e6, -1.1e6, 7.162, -263376.3044 };

/*                       An array of known isotope names                     */

const std::string IsoNames[87] = { "Hydrogen", "Deuterium", "Tritium", "Helium",
					"Lithium", "Beryllium", "Boron", "Carbon",
					"Nitrogen", "Oxygen", "Fluorine", "Neon",
					"Sodium", "Magnesium", "Aluminum", "Silicon",
					"Phosphorus", "Sulfur", "Chlorine", "Potassium",
					"Calcium", "Scandium", "Titanium", "Vanadium",
					"Chromium", "Manganese", "Iron", "Cobalt",
					"Nickel", "Copper", "Zinc", "Gallium",
					"Germanium", "Arsenic", "Selenium", "Bromine",
					"Rubidium", "Strontium", "Yttrium", "Zirconium",
					"Niobium", "Molybdenum", "Technetium", "Ruthenium",
					"Rhodium", "Palladium", "Silver", "Cadmium",
					"Indium", "Tin", "Antimony", "Tellurium",
					"Iodine", "Xenon", "Cesium", "Barium",
					"Lanthanum", "Praseodymium", "Neodymium", "Samarium",
					"Europium", "Gadolinium", "Terbium", "Dysprosium",
					"Holmium", "Erbium", "Thulium", "Ytterbium",
					"Lutetium", "Hafnium", "Tantalum", "Tungsten",
					"Rhenium", "Osmium", "Iridium", "Platinum",
					"Gold", "Mercury", "Thallium", "Lead",
					"Bismuth", "Astatine", "Actinium", "Thorium",
					"Protactinium", "Uranium", "Electron" };

/*                       Known Isotope Indices For Previous Names Array   */

const int IsoNamesIdx[NIso] = { 0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 6,
					6, 7, 7, 8, 8, 9, 9, 9, 10, 11, 11,
					11, 12, 13, 14, 15, 16, 17, 18, 18, 19, 19,
					20, 21, 22, 22, 23, 23, 24, 25, 26, 27, 28,
					29, 29, 30, 31, 31, 32, 33, 34, 35, 35, 36,
					36, 37, 38, 39, 40, 41, 41, 42, 43, 43, 44,
					45, 46, 46, 47, 47, 48, 48, 49, 49, 50, 50,
					51, 51, 52, 53, 53, 54, 55, 55, 56, 56, 57,
					58, 58, 59, 59, 60, 60, 61, 61, 62, 63, 63,
					64, 65, 66, 67, 67, 68, 68, 69, 69, 70, 71,
					72, 72, 73, 73, 74, 74, 75, 76, 77, 77, 78,
					78, 79, 80, 81, 82, 83, 84, 84, 85, 86 };

const std::string IsoElements[85] = { "H", "He", "Li", "Be", "B", "C", "N", "O",
									   "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
									   "Cl", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
									   "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
									   "Se", "Br", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
									   "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
									   "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Pr",
									   "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
									   "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
									   "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "At",
									   "Ac", "Th", "Pa", "U", "e" };

const int IsoIndexElements[NIso] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 4,
									 4, 5, 5, 6, 6, 7, 7, 7, 8, 9, 9,
					 				 9, 10, 11, 12, 13, 14, 15, 16, 16, 17, 17,
									 18, 19, 20, 20, 21, 21, 22, 23, 24, 25, 26,
									 27, 27, 28, 29, 29, 30, 31, 32, 33, 33, 34,
									 34, 35, 36, 37, 38, 39, 39, 40, 41, 41, 42,
									 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 48,
									 49, 49, 50, 51, 51, 52, 53, 53, 54, 54, 55,
									 56, 56, 57, 57, 58, 58, 59, 59, 60, 61, 61,
									 62, 63, 64, 65, 65, 66, 66, 67, 67, 68, 69,
									 70, 70, 71, 71, 72, 72, 73, 74, 75, 75, 76,
									 76, 77, 78, 79, 80, 81, 82, 82, 83, 84 };

const std::string IsoSymbol[NIso] = {"1H","2H","3H","3He","4He","6He","6Li","7Li",
									"8Li","9Be","10B","11B","12C","13C","14N","15N",
									"16O","17O","18O","19F","20Ne","21Ne","21Ne","23Na",
									"25Mg","27Al","29Si","31P","33S","35Cl","37Cl","39K",
									"41K","43Ca","45Sc","47Ti","49Ti","50V","51V","53Cr",
									"55Mn","57Fe","59Co","61Ni","63Cu","65Cu","67Zn","69Ga",
									"71Ga","73Ge","75As","77Se","79Br","81Br","85Rb","87Rb",
									"87Sr","89Y","91Zr","93Nb","95Mo","97Mo","99Tc","99Ru",
									"101Ru","103Rh","105Pd","107Ag","109Ag","111Cd","113Cd","113In",
									"115In","117Sn","119Sn","121Sb","123Sb","123Te","125Te","127I",
									"129Xe","131Xe","133Cs","135Ba","137Ba","138La","139La","141Pr",
									"143Nd","145Nd","147Sm","149Sm","151Eu","153Eu","155Gd","157Gd",
									"159Tb","161Dy","163Dy","165Ho","167Er","169Tm","171Yb","173Yb",
									"175Lu","176Lu","177Hf","179Hf","181Ta","183W","185Re","187Re",
									"187Os","189Os","191Ir","193Ir","195Pt","197Au","199Hg","201Hg",
									"203Tl","205Tl","207Pb","209Bi","211At","227Ac","229Th","231Pa",
									"233Pa","235U","0e"};




const double RelRF1H=400.13;
const double gammTesla=GAMMA1H/RelRF1H;
const double gammGauss=GAMMA1H/RelRF1H/10000.;


const double IsoGammaPerTesla[NIso]={IsoRelFreqs[0]*gammTesla,IsoRelFreqs[1]*gammTesla,IsoRelFreqs[2]*gammTesla,IsoRelFreqs[3]*gammTesla,IsoRelFreqs[4]*gammTesla,
IsoRelFreqs[5]*gammTesla,IsoRelFreqs[6]*gammTesla,IsoRelFreqs[7]*gammTesla,IsoRelFreqs[8]*gammTesla,
IsoRelFreqs[9]*gammTesla,IsoRelFreqs[10]*gammTesla,IsoRelFreqs[11]*gammTesla,IsoRelFreqs[12]*gammTesla,
IsoRelFreqs[13]*gammTesla,IsoRelFreqs[14]*gammTesla,IsoRelFreqs[15]*gammTesla,IsoRelFreqs[16]*gammTesla,
IsoRelFreqs[17]*gammTesla,IsoRelFreqs[18]*gammTesla,IsoRelFreqs[19]*gammTesla,IsoRelFreqs[20]*gammTesla,
IsoRelFreqs[21]*gammTesla,IsoRelFreqs[22]*gammTesla,IsoRelFreqs[23]*gammTesla,IsoRelFreqs[24]*gammTesla,
IsoRelFreqs[25]*gammTesla,IsoRelFreqs[26]*gammTesla,IsoRelFreqs[27]*gammTesla,IsoRelFreqs[28]*gammTesla,
IsoRelFreqs[29]*gammTesla,IsoRelFreqs[30]*gammTesla,IsoRelFreqs[31]*gammTesla,IsoRelFreqs[32]*gammTesla,
IsoRelFreqs[33]*gammTesla,IsoRelFreqs[34]*gammTesla,IsoRelFreqs[35]*gammTesla,IsoRelFreqs[36]*gammTesla,
IsoRelFreqs[37]*gammTesla,IsoRelFreqs[38]*gammTesla,IsoRelFreqs[39]*gammTesla,IsoRelFreqs[40]*gammTesla,
IsoRelFreqs[41]*gammTesla,IsoRelFreqs[42]*gammTesla,IsoRelFreqs[43]*gammTesla,IsoRelFreqs[44]*gammTesla,
IsoRelFreqs[45]*gammTesla,IsoRelFreqs[46]*gammTesla,IsoRelFreqs[47]*gammTesla,IsoRelFreqs[48]*gammTesla,
IsoRelFreqs[49]*gammTesla,IsoRelFreqs[50]*gammTesla,IsoRelFreqs[51]*gammTesla,IsoRelFreqs[52]*gammTesla,
IsoRelFreqs[53]*gammTesla,IsoRelFreqs[54]*gammTesla,IsoRelFreqs[55]*gammTesla,IsoRelFreqs[56]*gammTesla,
IsoRelFreqs[57]*gammTesla,IsoRelFreqs[58]*gammTesla,IsoRelFreqs[59]*gammTesla,IsoRelFreqs[60]*gammTesla,
IsoRelFreqs[61]*gammTesla,IsoRelFreqs[62]*gammTesla,IsoRelFreqs[63]*gammTesla,IsoRelFreqs[64]*gammTesla,
IsoRelFreqs[65]*gammTesla,IsoRelFreqs[66]*gammTesla,IsoRelFreqs[67]*gammTesla,IsoRelFreqs[68]*gammTesla,
IsoRelFreqs[69]*gammTesla,IsoRelFreqs[70]*gammTesla,IsoRelFreqs[71]*gammTesla,IsoRelFreqs[72]*gammTesla,
IsoRelFreqs[73]*gammTesla,IsoRelFreqs[74]*gammTesla,IsoRelFreqs[75]*gammTesla,IsoRelFreqs[76]*gammTesla,
IsoRelFreqs[77]*gammTesla,IsoRelFreqs[78]*gammTesla,IsoRelFreqs[79]*gammTesla,IsoRelFreqs[80]*gammTesla,
IsoRelFreqs[81]*gammTesla,IsoRelFreqs[82]*gammTesla,IsoRelFreqs[83]*gammTesla,IsoRelFreqs[84]*gammTesla,
IsoRelFreqs[85]*gammTesla,IsoRelFreqs[86]*gammTesla,IsoRelFreqs[87]*gammTesla,IsoRelFreqs[88]*gammTesla,
IsoRelFreqs[89]*gammTesla,IsoRelFreqs[90]*gammTesla,IsoRelFreqs[91]*gammTesla,IsoRelFreqs[92]*gammTesla,
IsoRelFreqs[93]*gammTesla,IsoRelFreqs[94]*gammTesla,IsoRelFreqs[95]*gammTesla,IsoRelFreqs[96]*gammTesla,
IsoRelFreqs[97]*gammTesla,IsoRelFreqs[98]*gammTesla,IsoRelFreqs[99]*gammTesla,IsoRelFreqs[100]*gammTesla,
IsoRelFreqs[101]*gammTesla,IsoRelFreqs[102]*gammTesla,IsoRelFreqs[103]*gammTesla,IsoRelFreqs[104]*gammTesla,
IsoRelFreqs[105]*gammTesla,IsoRelFreqs[106]*gammTesla,IsoRelFreqs[107]*gammTesla,IsoRelFreqs[108]*gammTesla,
IsoRelFreqs[109]*gammTesla,IsoRelFreqs[110]*gammTesla,IsoRelFreqs[111]*gammTesla,IsoRelFreqs[112]*gammTesla,
IsoRelFreqs[113]*gammTesla,IsoRelFreqs[114]*gammTesla,IsoRelFreqs[115]*gammTesla,IsoRelFreqs[116]*gammTesla,
IsoRelFreqs[117]*gammTesla,IsoRelFreqs[118]*gammTesla,IsoRelFreqs[119]*gammTesla,IsoRelFreqs[120]*gammTesla,
IsoRelFreqs[121]*gammTesla,IsoRelFreqs[122]*gammTesla,IsoRelFreqs[123]*gammTesla,IsoRelFreqs[124]*gammTesla,
IsoRelFreqs[125]*gammTesla,IsoRelFreqs[126]*gammTesla,IsoRelFreqs[127]*gammTesla,IsoRelFreqs[128]*gammTesla,
IsoRelFreqs[129]*gammTesla,IsoRelFreqs[130]*gammTesla};

const double IsoGammaPerGauss[NIso]={IsoRelFreqs[0]*gammGauss,IsoRelFreqs[1]*gammGauss,IsoRelFreqs[2]*gammGauss,IsoRelFreqs[3]*gammGauss,IsoRelFreqs[4]*gammGauss,
IsoRelFreqs[5]*gammGauss,IsoRelFreqs[6]*gammGauss,IsoRelFreqs[7]*gammGauss,IsoRelFreqs[8]*gammGauss,
IsoRelFreqs[9]*gammGauss,IsoRelFreqs[10]*gammGauss,IsoRelFreqs[11]*gammGauss,IsoRelFreqs[12]*gammGauss,
IsoRelFreqs[13]*gammGauss,IsoRelFreqs[14]*gammGauss,IsoRelFreqs[15]*gammGauss,IsoRelFreqs[16]*gammGauss,
IsoRelFreqs[17]*gammGauss,IsoRelFreqs[18]*gammGauss,IsoRelFreqs[19]*gammGauss,IsoRelFreqs[20]*gammGauss,
IsoRelFreqs[21]*gammGauss,IsoRelFreqs[22]*gammGauss,IsoRelFreqs[23]*gammGauss,IsoRelFreqs[24]*gammGauss,
IsoRelFreqs[25]*gammGauss,IsoRelFreqs[26]*gammGauss,IsoRelFreqs[27]*gammGauss,IsoRelFreqs[28]*gammGauss,
IsoRelFreqs[29]*gammGauss,IsoRelFreqs[30]*gammGauss,IsoRelFreqs[31]*gammGauss,IsoRelFreqs[32]*gammGauss,
IsoRelFreqs[33]*gammGauss,IsoRelFreqs[34]*gammGauss,IsoRelFreqs[35]*gammGauss,IsoRelFreqs[36]*gammGauss,
IsoRelFreqs[37]*gammGauss,IsoRelFreqs[38]*gammGauss,IsoRelFreqs[39]*gammGauss,IsoRelFreqs[40]*gammGauss,
IsoRelFreqs[41]*gammGauss,IsoRelFreqs[42]*gammGauss,IsoRelFreqs[43]*gammGauss,IsoRelFreqs[44]*gammGauss,
IsoRelFreqs[45]*gammGauss,IsoRelFreqs[46]*gammGauss,IsoRelFreqs[47]*gammGauss,IsoRelFreqs[48]*gammGauss,
IsoRelFreqs[49]*gammGauss,IsoRelFreqs[50]*gammGauss,IsoRelFreqs[51]*gammGauss,IsoRelFreqs[52]*gammGauss,
IsoRelFreqs[53]*gammGauss,IsoRelFreqs[54]*gammGauss,IsoRelFreqs[55]*gammGauss,IsoRelFreqs[56]*gammGauss,
IsoRelFreqs[57]*gammGauss,IsoRelFreqs[58]*gammGauss,IsoRelFreqs[59]*gammGauss,IsoRelFreqs[60]*gammGauss,
IsoRelFreqs[61]*gammGauss,IsoRelFreqs[62]*gammGauss,IsoRelFreqs[63]*gammGauss,IsoRelFreqs[64]*gammGauss,
IsoRelFreqs[65]*gammGauss,IsoRelFreqs[66]*gammGauss,IsoRelFreqs[67]*gammGauss,IsoRelFreqs[68]*gammGauss,
IsoRelFreqs[69]*gammGauss,IsoRelFreqs[70]*gammGauss,IsoRelFreqs[71]*gammGauss,IsoRelFreqs[72]*gammGauss,
IsoRelFreqs[73]*gammGauss,IsoRelFreqs[74]*gammGauss,IsoRelFreqs[75]*gammGauss,IsoRelFreqs[76]*gammGauss,
IsoRelFreqs[77]*gammGauss,IsoRelFreqs[78]*gammGauss,IsoRelFreqs[79]*gammGauss,IsoRelFreqs[80]*gammGauss,
IsoRelFreqs[81]*gammGauss,IsoRelFreqs[82]*gammGauss,IsoRelFreqs[83]*gammGauss,IsoRelFreqs[84]*gammGauss,
IsoRelFreqs[85]*gammGauss,IsoRelFreqs[86]*gammGauss,IsoRelFreqs[87]*gammGauss,IsoRelFreqs[88]*gammGauss,
IsoRelFreqs[89]*gammGauss,IsoRelFreqs[90]*gammGauss,IsoRelFreqs[91]*gammGauss,IsoRelFreqs[92]*gammGauss,
IsoRelFreqs[93]*gammGauss,IsoRelFreqs[94]*gammGauss,IsoRelFreqs[95]*gammGauss,IsoRelFreqs[96]*gammGauss,
IsoRelFreqs[97]*gammGauss,IsoRelFreqs[98]*gammGauss,IsoRelFreqs[99]*gammGauss,IsoRelFreqs[100]*gammGauss,
IsoRelFreqs[101]*gammGauss,IsoRelFreqs[102]*gammGauss,IsoRelFreqs[103]*gammGauss,IsoRelFreqs[104]*gammGauss,
IsoRelFreqs[105]*gammGauss,IsoRelFreqs[106]*gammGauss,IsoRelFreqs[107]*gammGauss,IsoRelFreqs[108]*gammGauss,
IsoRelFreqs[109]*gammGauss,IsoRelFreqs[110]*gammGauss,IsoRelFreqs[111]*gammGauss,IsoRelFreqs[112]*gammGauss,
IsoRelFreqs[113]*gammGauss,IsoRelFreqs[114]*gammGauss,IsoRelFreqs[115]*gammGauss,IsoRelFreqs[116]*gammGauss,
IsoRelFreqs[117]*gammGauss,IsoRelFreqs[118]*gammGauss,IsoRelFreqs[119]*gammGauss,IsoRelFreqs[120]*gammGauss,
IsoRelFreqs[121]*gammGauss,IsoRelFreqs[122]*gammGauss,IsoRelFreqs[123]*gammGauss,IsoRelFreqs[124]*gammGauss,
IsoRelFreqs[125]*gammGauss,IsoRelFreqs[126]*gammGauss,IsoRelFreqs[127]*gammGauss,IsoRelFreqs[128]*gammGauss,
IsoRelFreqs[129]*gammGauss,IsoRelFreqs[130]*gammGauss};


double Isotope::RelativeRF_1H = RelRF1H;

Isotope::Isotope(){	Index=0; }

Isotope::Isotope(const std::string &ss){
	Index=getCorrectIndex(ss);
}


Isotope::Isotope(const  Isotope &rhs){
	Index=rhs.Index;
}

Isotope &Isotope::operator=(const Isotope &rhs){
	if(&rhs==this) return *this;
	Index=rhs.Index;
	return *this;
}

Isotope &Isotope::operator=(const std::string &rhs){
	Index=getCorrectIndex(rhs);
	return *this;
}


//hilbert space size
int Isotope::HS()const{	return IsoSpins[Index];	}

double Isotope::spin()const{ return (IsoSpins[Index])?(IsoSpins[Index]-1)/2.:0.;	}		//spin quantum number

double Isotope::qn()const{ return spin();	}		//spin quantum number

int Isotope::number()const{ return  IsoNumbers[Index];	}		//atomic number
double Isotope::mass()const{ return IsoMasses[Index];	}		//amu mass
double Isotope::weight()const{ return IsoWeights[Index];	}	//g/mole weight
double Isotope::receptivity()const{ return IsoRecepts[Index];	}	//relative sensitivity
std::string Isotope::name()const{ return IsoNames[IsoNamesIdx[Index]];	}	//"Hydrogen", "Carbon", etc
std::string Isotope::element()const{ return IsoElements[IsoIndexElements[Index]];	}	//"H", "C", etc
std::string Isotope::symbol() const{ return IsoSymbol[Index];	}	//"1H", "2H", etc
double Isotope::gamma()const{ return IsoGammaPerTesla[Index]; 	}		//gamma factor*2Pi (i.e.in rads/Telsa)
double Isotope::gammaGauss() const { return IsoGammaPerGauss[Index]; 	}	//gamma factor*2*pi (i.e. rad/GAUSS)
double Isotope::relativeFrequency()const{ return IsoRelFreqs[Index]/RelRF1H; }	//frequency relative to "1H" frequency

std::string Isotope::momentum() const {
  if(!HS()) return std::string("0");
  else if(HS()%2) return itost((HS()-1)/2);
  else return std::string(itost(HS()-1))+"/2";
}

void Isotope::SpinType(char *sname){
	setCorrectIndex(std::string(sname));
}

void Isotope::SpinType(std::string sname){
	setCorrectIndex(sname);
}

void Isotope::setCorrectIndex(std::string sname){
	Index=getCorrectIndex(sname);
}

void Isotope::Spin(std::string in)
{
	Index=getCorrectIndex(in);
}


int Isotope::getCorrectIndex(std::string ss){
	std::string tm;
	for(int i=0;i<NIso;i++){
		tm=(itost(IsoMasses[i])+IsoElements[IsoIndexElements[i]]);
		if(ss==tm){
			return i;
		}
	}
	std::cerr<<std::endl<<"Error::Isotope(std::string)"<<std::endl;
	std::cerr<<" Isotope \""<<ss<<"\" Not found"<<std::endl;
	std::cerr<<" Setting to \"1H\"..."<<std::endl;
	return 0;
}

void Isotope::print(std::ostream &oo){
	oo<<(*this);
}

std::ostream &operator<<(std::ostream &oo, const Isotope &out){
	oo<<"Isotope: "<<out.symbol()<<"(momentum="<<out.momentum();
	oo<<", gamma="<<out.gamma()<<") ";
	return oo;
}

END_BL_NAMESPACE



#endif

