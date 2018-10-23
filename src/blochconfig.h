/* src/blochconfig.h.  Generated from blochconfig.h.in by configure.  */
/* config.h.in.  Generated automatically from configure.in by autoheader.  */

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define as __inline if that's what the C compiler calls it.  */
/* #undef inline */

/* Define if you need to in order for stat and other things to work.  */
/* #undef _POSIX_SOURCE */

/* Define to `unsigned' if <sys/types.h> doesn't define.  */
/* #undef size_t */

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if you can safely include both <sys/time.h> and <time.h>.  */
#define TIME_WITH_SYS_TIME 1

/* Define if you have the <sys/time.h> header file.  */
#define HAVE_SYS_TIME_H 1

/* Define if you have the <unistd.h> header file.  */
#define HAVE_UNISTD_H 1

/* Define if you have the <climits> header file.  */
#define HAVE_CLIMITS 1

/* Define if you have the <limits.h> header file.  */
#define HAVE_LIMITS_H 1

/* Define if you have the function isnan.  */
#define HAVE_ISNAN 1

/* Define if you have the function finite.  */
#define HAVE_FINITE 1

/* Name of package */
#define PACKAGE "blochlib"

/* Version number of package */
#define VERSION "1.2"

/* Useing minuit or not */
#define USE_MINUIT 1

/* Useing f2c fortran or not */
/* #undef f2cFortran */

/* do you have the 'getline' function */
#define HAVE_GETLINE 1

/* do you have MPI library */
/* #undef HAVE_MPI */

/* using an auxilliary BLAS lib for dgemm */
#define AUX_BLAS_LIB 1

/* do we have the max(thing) function */
#define HAVE_MAX1 1

/* do we have the max(thing, thing) function */
/* #undef HAVE_MAX2 */

/* do we have the min(thing) function */
#define HAVE_MIN1 1

/* do we have the min(thing, thing) function */
/* #undef HAVE_MIN2 */

/* do we have the sign(thing) function */
#define HAVE_SIGN 1

/* Use the FFTW fast Fourier Transform Package at all */
/* #undef HAVE_FFTW */

/* Use the FFTW fast Fourier Transform Package Version 2 */
/* #undef HAVE_FFTW_II */

/* Use the FFTW v 3 fast Fourier Transform Package */
/* #undef HAVE_FFTW_III */

/* Are we going to use exceptions? */
/* #undef USE_EXCEPTIONS */


/* Do We have signal.h header */
#define HAVE_SIGNAL_H 1

/* the namespace name */
#define BL_NAMESPACE BlochLib

/* if namespaces are allowed */
#define BEGIN_BL_NAMESPACE namespace BL_NAMESPACE {

/* if namespaces are allowed */
#define END_BL_NAMESPACE }

