/*
 *  Interface to minuit using cfortran.h
 *
 *  Edit history:
 *  G.Folger  12-Dec-94  change some to use ROUTINE for passing fcn/futil

 *  Edit history:
 *  Bo Blanton 01-30-02 added the PROTOCALLSUBs to allow compilation
 *			   with C++
 */
#ifndef _MINUIT_H_WRAP__
#define _MINUIT_H_WRAP__ 1


#include "cfortran.h"
#include "minuitfcn.h"


/*------------------------------------------------------------------
fortran filename   : bl_init.f (blochlib initializartion call for mninit(int, int,int)
------------------------------------------------------------------*/

PROTOCCALLSFSUB1(BL_MNINIT,bl_mninit,STRING)
#define BL_MNINIT(A1)  CCALLSFSUB2(BL_MNINIT,bl_mninit,STRING, A1)


/*------------------------------------------------------------------
fortran filename   : minuit.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MINUIT,minuit,ROUTINE,ROUTINE)
#define MINUIT(A1,A2)  CCALLSFSUB2(MINUIT,minuit,ROUTINE,ROUTINE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnamin.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNAMIN,mnamin,DOUBLE,DOUBLE)
#define MNAMIN(A1,A2)  CCALLSFSUB2(MNAMIN,mnamin,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnbins.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB7(MNBINS,mnbins,DOUBLE,DOUBLE,INT,PDOUBLE,PDOUBLE,PINT,PDOUBLE)
#define MNBINS(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(MNBINS,mnbins,DOUBLE,DOUBLE,INT,PDOUBLE,PDOUBLE,PINT,PDOUBLE,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : mncalf.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB4(MNCALF,mncalf,DOUBLE,DOUBLEV,PDOUBLE,DOUBLE)
#define MNCALF(A1,A2,A3,A4)  CCALLSFSUB4(MNCALF,mncalf,DOUBLE,DOUBLEV,PDOUBLE,DOUBLE,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : mncler.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB0(MNCLER,mncler)
#define MNCLER() CCALLSFSUB0(MNCLER,mncler)

/*------------------------------------------------------------------
fortran filename   : mncntr.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB5(MNCNTR,mncntr,DOUBLE,INT,INT,PINT,DOUBLE)
#define MNCNTR(A1,A2,A3,A4,A5)  CCALLSFSUB5(MNCNTR,mncntr,DOUBLE,INT,INT,PINT,DOUBLE,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : mncomd.f
------------------------------------------------------------------*/
/*
#define mncomd_ELEMS_2          ZTRINGV_NUM(1)
#define mncomd_ELEMLEN_2        ZTRINGV_NUM(255)
*/

PROTOCCALLSFSUB4(MNCOMD,mncomd,ROUTINE,STRING,PINT,ROUTINE)
#define MNCOMD(A1,A2,A3,A4)  CCALLSFSUB4(MNCOMD,mncomd,ROUTINE,STRING,PINT,ROUTINE,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : mncont.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB8(MNCONT,mncont,ROUTINE,INT,INT,INT,PDOUBLE,PDOUBLE,PINT,ROUTINE)
#define MNCONT(A1,A2,A3,A4,A5,A6,A7,A8)  CCALLSFSUB8(MNCONT,mncont,ROUTINE,INT,INT,INT,PDOUBLE,PDOUBLE,PINT,ROUTINE,A1,A2,A3,A4,A5,A6,A7,A8)

/*-----------------------------------------------------------------
fortran filename   : mncrck.f
------------------------------------------------------------------*/
/*
#define mncrck_ELEMS_1          ZTRINGV_NUM(1)
#define mncrck_ELEMLEN_1        ZTRINGV_NUM(255)
#define mncrck_ELEMS_3          ZTRINGV_NUM(1)
#define mncrck_ELEMLEN_3        ZTRINGV_NUM(255)
*/

PROTOCCALLSFSUB9(MNCRCK,mncrck,STRING,INT,PSTRING,PINT,INT,PDOUBLE,PINT,PINT,INT)
#define MNCRCK(A1,A2,A3,A4,A5,A6,A7,A8,A9)  CCALLSFSUB9(MNCRCK,mncrck,STRING,INT,PSTRING,PINT,INT,PDOUBLE,PINT,PINT,INT,A1,A2,A3,A4,A5,A6,A7,A8,A9)

/*------------------------------------------------------------------
fortran filename   : mncros.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB4(MNCROS,mncros,DOUBLE,PDOUBLE,PINT,DOUBLE)
#define MNCROS(A1,A2,A3,A4)  CCALLSFSUB4(MNCROS,mncros,DOUBLE,PDOUBLE,PINT,DOUBLE,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : mncuve.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNCUVE,mncuve,DOUBLE,DOUBLE)
#define MNCUVE(A1,A2)  CCALLSFSUB2(MNCUVE,mncuve,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnderi.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNDERI,mnderi,DOUBLE,DOUBLE)
#define MNDERI(A1,A2)  CCALLSFSUB2(MNDERI,mnderi,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mndxdi.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB3(MNDXDI,mndxdi,DOUBLE,INT,PDOUBLE)
#define MNDXDI(A1,A2,A3)  CCALLSFSUB3(MNDXDI,mndxdi,DOUBLE,INT,PDOUBLE,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : mneig.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB7(MNEIG,mneig,PDOUBLE,INT,INT,INT,PDOUBLE,DOUBLE,PINT)
#define MNEIG(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(MNEIG,mneig,PDOUBLE,INT,INT,INT,PDOUBLE,DOUBLE,PINT,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : mnemat.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNEMAT,mnemat,PDOUBLE,INT)
#define MNEMAT(A1,A2)  CCALLSFSUB2(MNEMAT,mnemat,PDOUBLE,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnerrs.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB5(MNERRS,mnerrs,INT,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define MNERRS(A1,A2,A3,A4,A5)  CCALLSFSUB5(MNERRS,mnerrs,INT,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : mneval.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB5(MNEVAL,mneval,DOUBLE,DOUBLE,PDOUBLE,PINT,DOUBLE)
#define MNEVAL(A1,A2,A3,A4,A5)  CCALLSFSUB5(MNEVAL,mneval,DOUBLE,DOUBLE,PDOUBLE,PINT,DOUBLE,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : mnexcm.f
------------------------------------------------------------------*/
/*
#define mnexcm_ELEMS_2          ZTRINGV_NUM(1)
#define mnexcm_ELEMLEN_2        ZTRINGV_NUM(255)
*/

PROTOCCALLSFSUB6(MNEXCM,mnexcm,ROUTINE,STRING,DOUBLEV,INT,PINT,ROUTINE)
#define MNEXCM(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(MNEXCM,mnexcm,ROUTINE,STRING,DOUBLEV,INT,PINT,ROUTINE,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : mnexin.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB1(MNEXIN,mnexin,PDOUBLE)
#define MNEXIN(A1)  CCALLSFSUB1(MNEXIN,mnexin,PDOUBLE,A1)

/*------------------------------------------------------------------
fortran filename   : mnfixp.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNFIXP,mnfixp,INT,PINT)
#define MNFIXP(A1,A2)  CCALLSFSUB2(MNFIXP,mnfixp,INT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnfree.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB1(MNFREE,mnfree,INT)
#define MNFREE(A1)  CCALLSFSUB1(MNFREE,mnfree,INT,A1)

/*------------------------------------------------------------------
fortran filename   : mngrad.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNGRAD,mngrad,DOUBLE,DOUBLE)
#define MNGRAD(A1,A2)  CCALLSFSUB2(MNGRAD,mngrad,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnhelp.f
------------------------------------------------------------------*/
/*
#define mnhelp_ELEMS_1          ZTRINGV_NUM(1)
#define mnhelp_ELEMLEN_1        ZTRINGV_NUM(255)
*/

PROTOCCALLSFSUB2(MNHELP,mnhelp,STRING,INT)
#define MNHELP(A1,A2)  CCALLSFSUB2(MNHELP,mnhelp,STRING,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnhes1.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNHES1,mnhes1,DOUBLE,DOUBLE)
#define MNHES1(A1,A2)  CCALLSFSUB2(MNHES1,mnhes1,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnhess.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNHESS,mnhess,DOUBLE,DOUBLE)
#define MNHESS(A1,A2)  CCALLSFSUB2(MNHESS,mnhess,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnimpr.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNIMPR,mnimpr,DOUBLE,DOUBLE)
#define MNIMPR(A1,A2)  CCALLSFSUB2(MNIMPR,mnimpr,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mninex.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB1(MNINEX,mninex,DOUBLEV)
#define MNINEX(A1)  CCALLSFSUB1(MNINEX,mninex,DOUBLEV,A1)

/*------------------------------------------------------------------
fortran filename   : mninit.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB3(MNINIT,mninit,INT,INT,INT)
#define MNINIT(A1,A2,A3)  CCALLSFSUB3(MNINIT,mninit,INT,INT,INT,A1,A2,A3)
/*------------------------------------------------------------------
fortran filename   : mninpu.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNINPU,mninpu,INT,PINT)
#define MNINPU(A1,A2)  CCALLSFSUB2(MNINPU,mninpu,INT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnintr.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNINTR,mnintr,ROUTINE,ROUTINE)
#define MNINTR(A1,A2)  CCALLSFSUB2(MNINTR,mnintr,ROUTINE,ROUTINE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnlims.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNLIMS,mnlims,DOUBLE,DOUBLE)
#define MNLIMS(A1,A2)  CCALLSFSUB2(MNLIMS,mnlims,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnline.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB7(MNLINE,mnline,DOUBLE,DOUBLEV,DOUBLE,DOUBLEV,DOUBLE,DOUBLE,DOUBLE)
#define MNLINE(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(MNLINE,mnline,DOUBLE,DOUBLEV,DOUBLE,DOUBLEV,DOUBLE,DOUBLE,DOUBLE,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : mnmatu.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB1(MNMATU,mnmatu,INT)
#define MNMATU(A1)  CCALLSFSUB1(MNMATU,mnmatu,INT,A1)

/*------------------------------------------------------------------
fortran filename   : mnmigr.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNMIGR,mnmigr,DOUBLE,DOUBLE)
#define MNMIGR(A1,A2)  CCALLSFSUB2(MNMIGR,mnmigr,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnmnos.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNMNOS,mnmnos,DOUBLE,DOUBLE)
#define MNMNOS(A1,A2)  CCALLSFSUB2(MNMNOS,mnmnos,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnmnot.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB6(MNMNOT,mnmnot,DOUBLE,INT,INT,DOUBLE,DOUBLE,DOUBLE)
#define MNMNOT(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(MNMNOT,mnmnot,DOUBLE,INT,INT,DOUBLE,DOUBLE,DOUBLE,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : mnparm.f
------------------------------------------------------------------*/
/*
#define mnparm_ELEMS_2          ZTRINGV_NUM(1)
#define mnparm_ELEMLEN_2        ZTRINGV_NUM(255)
*/

PROTOCCALLSFSUB7(MNPARM,mnparm,INT,STRING,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PINT)
#define MNPARM(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(MNPARM,mnparm,INT,STRING,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PINT,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : mnpars.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNPARS,mnpars,STRING,PINT)
#define MNPARS(A1,A2)  CCALLSFSUB2(MNPARS,mnpars,STRING,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnpfit.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB5(MNPFIT,mnpfit,DOUBLEV,DOUBLEV,INT,PDOUBLE,PDOUBLE)
#define MNPFIT(A1,A2,A3,A4,A5)  CCALLSFSUB5(MNPFIT,mnpfit,DOUBLEV,DOUBLEV,INT,PDOUBLE,PDOUBLE,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : mnpint.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB3(MNPINT,mnpint,PDOUBLE,INT,PDOUBLE)
#define MNPINT(A1,A2,A3)  CCALLSFSUB3(MNPINT,mnpint,PDOUBLE,INT,PDOUBLE,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : mnplot.f
------------------------------------------------------------------*/
/*
#define mnplot_ELEMS_3          ZTRINGV_NUM(#?ð#)
#define mnplot_ELEMLEN_3        ZTRINGV_NUM(1)
*/

PROTOCCALLSFSUB7(MNPLOT,mnplot,PDOUBLE,PDOUBLE,PSTRINGV,INT,INT,INT,INT)
#define MNPLOT(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(MNPLOT,mnplot,PDOUBLE,PDOUBLE,PSTRINGV,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : mnpout.f
------------------------------------------------------------------*/
/*
#define mnpout_ELEMS_2          ZTRINGV_NUM(1)
#define mnpout_ELEMLEN_2        ZTRINGV_NUM(255)
*/

PROTOCCALLSFSUB7(MNPOUT,mnpout,INT,PSTRING,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PINT)
#define MNPOUT(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(MNPOUT,mnpout,INT,PSTRING,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PINT,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : mnprin.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNPRIN,mnprin,INT,DOUBLE)
#define MNPRIN(A1,A2)  CCALLSFSUB2(MNPRIN,mnprin,INT,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnpsdf.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB0(MNPSDF,mnpsdf)
#define MNPSDF() CCALLSFSUB0(MNPSDF,mnpsdf)

/*------------------------------------------------------------------
fortran filename   : mnrazz.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB5(MNRAZZ,mnrazz,DOUBLE,DOUBLEV,PDOUBLE,PINT,PINT)
#define MNRAZZ(A1,A2,A3,A4,A5)  CCALLSFSUB5(MNRAZZ,mnrazz,DOUBLE,DOUBLEV,PDOUBLE,PINT,PINT,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : mnread.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB4(MNREAD,mnread,DOUBLE,INT,PINT,DOUBLE)
#define MNREAD(A1,A2,A3,A4)  CCALLSFSUB4(MNREAD,mnread,DOUBLE,INT,PINT,DOUBLE,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : mnrn15.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNRN15,mnrn15,PDOUBLE,PINT)
#define MNRN15(A1,A2)  CCALLSFSUB2(MNRN15,mnrn15,PDOUBLE,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnrset.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB1(MNRSET,mnrset,INT)
#define MNRSET(A1)  CCALLSFSUB1(MNRSET,mnrset,INT,A1)

/*------------------------------------------------------------------
fortran filename   : mnsave.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB0(MNSAVE,mnsave)
#define MNSAVE() CCALLSFSUB0(MNSAVE,mnsave)

/*------------------------------------------------------------------
fortran filename   : mnscan.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNSCAN,mnscan,DOUBLE,DOUBLE)
#define MNSCAN(A1,A2)  CCALLSFSUB2(MNSCAN,mnscan,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnseek.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNSEEK,mnseek,DOUBLE,DOUBLE)
#define MNSEEK(A1,A2)  CCALLSFSUB2(MNSEEK,mnseek,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnset.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNSET,mnset,DOUBLE,DOUBLE)
#define MNSET(A1,A2)  CCALLSFSUB2(MNSET,mnset,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnseti.f
------------------------------------------------------------------*/
/*
#define mnseti_ELEMS_1          ZTRINGV_NUM(1)
#define mnseti_ELEMLEN_1        ZTRINGV_NUM(255)
*/

PROTOCCALLSFSUB1(MNSETI,mnseti,STRING)
#define MNSETI(A1)  CCALLSFSUB1(MNSETI,mnseti,STRING,A1)

/*------------------------------------------------------------------
fortran filename   : mnsimp.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNSIMP,mnsimp,DOUBLE,DOUBLE)
#define MNSIMP(A1,A2)  CCALLSFSUB2(MNSIMP,mnsimp,DOUBLE,DOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnstat.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB6(MNSTAT,mnstat,PDOUBLE,PDOUBLE,PDOUBLE,PINT,PINT,PINT)
#define MNSTAT(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(MNSTAT,mnstat,PDOUBLE,PDOUBLE,PDOUBLE,PINT,PINT,PINT,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : mnstin.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNSTIN,mnstin,BYTE,PINT)
#define MNSTIN(A1,A2)  CCALLSFSUB2(MNSTIN,mnstin,BYTE,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mntiny.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB2(MNTINY,mntiny,DOUBLE,PDOUBLE)
#define MNTINY(A1,A2)  CCALLSFSUB2(MNTINY,mntiny,DOUBLE,PDOUBLE,A1,A2)

/*------------------------------------------------------------------
fortran filename   : mnunpt.f
------------------------------------------------------------------*/

PROTOCCALLSFFUN1(LOGICAL,MNUNPT,mnunpt,BYTE)
//PROTOCCALLSFSUB1(MNUNPT,mnunpt,BYTE)
#define MNUNPT(A2)  CCALLSFFUN1(MNUNPT,mnunpt,BYTE,A2)

/*------------------------------------------------------------------
fortran filename   : mnvert.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB5(MNVERT,mnvert,PDOUBLE,INT,INT,INT,PINT)
#define MNVERT(A1,A2,A3,A4,A5)  CCALLSFSUB5(MNVERT,mnvert,PDOUBLE,INT,INT,INT,PINT,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : mnwarn.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB3(MNWARN,mnwarn,BYTE,BYTE,BYTE)
#define MNWARN(A1,A2,A3)  CCALLSFSUB3(MNWARN,mnwarn,BYTE,BYTE,BYTE,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : mnwerr.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB0(MNWERR,mnwerr)
#define MNWERR() CCALLSFSUB0(MNWERR,mnwerr)

/*------------------------------------------------------------------
fortran filename   : stand.f
------------------------------------------------------------------*/

PROTOCCALLSFSUB0(STAND,stand)
#define STAND() CCALLSFSUB0(STAND,stand)


#endif

