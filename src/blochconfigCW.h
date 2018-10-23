
//Config file for CodeWarrior ide.....

#define ON_WINDOWS 1

#ifndef __MINGW32__

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if you can safely include both <sys/time.h> and <time.h>.  */
#define TIME_WITH_SYS_TIME 0

/* Define if you have the <sys/time.h> header file.  */
#define HAVE_SYS_TIME_H 0

/* Define if you have the <unistd.h> header file.  */
#define HAVE_UNISTD_H 0

/* Name of package */
#define PACKAGE "blochlib"

/* Version number of package */
#define VERSION "0.9"

/* do we have the max(thing) function */
#define HAVE_MAX1 1

/* do we have the max(thing, thing) function */
/* #undef HAVE_MAX2 */

/* do we have the min(thing) function */
#define HAVE_MIN1 1

/* do we have the min(thing, thing) function */
/* #undef HAVE_MIN2 */

/* have the isnan function */
#undef HAVE_ISNAN

#endif
