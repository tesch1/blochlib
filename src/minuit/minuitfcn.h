/* C interface to minuit function */
/* Andr\'as Major 1997 */

#ifndef _MINUITFCN_H_WRAP__
#define _MINUITFCN_H_WRAP__ 1

//#ifndef ON_WINDOWS
//	#include "blochconfig.h"
//#endif

#include "cfortran.h"


void fcn(int npar, double* grad, double* fcnval, double* xval, int iflag, void* futil);

#if defined(CFORTRAN_HACK)
void minuitfcn(int npar, double* grad, double* fcnval, double* xval, int iflag, void* futil);
FCALLSCSUB6(minuitfcn,MINUITFCN,minuitfcn,INT,DOUBLEV,PDOUBLE,DOUBLEV,INT,PINT);
void minuitfcn(int npar, double* grad, double* fcnval, double* xval, int iflag, void* futil)
{
  fcn(npar, grad, fcnval, xval, iflag, futil);
}
#else
void minuitfcn(int* npar, double* grad, double* fcnval, double* xval, int* iflag, void* futil);
FCALLSCSUB6(minuitfcn,MINUITFCN,minuitfcn,PINT,DOUBLEV,PDOUBLE,DOUBLEV,PINT,PINT)
void minuitfcn(int* npar, double* grad, double* fcnval, double* xval, int* iflag, void* futil)
{
  fcn(*npar, grad, fcnval, xval, *iflag, futil);
}
#endif

int intrac_();

#endif
