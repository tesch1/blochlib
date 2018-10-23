


###################### BLOCHLIB FUNCTIONS ###########


AC_DEFUN(ACX_CHECK_CXX_FLAGS,
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(whether ${CXX-cc} accepts $1, ac_$2,
[echo 'void f(){}' > conftest.c
if test -z "`${CXX-cc} $1 -c conftest.c 2>&1`"; then
	ac_$2=yes
else
	ac_$2=no
fi
rm -f conftest*
])
if test "$ac_$2" = yes; then
	:
	$3
else
	:
	$4
fi
])

AC_DEFUN(ACX_CHECK_INLINE,
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(for inline , ok,
[echo 'inline int f(){}' > conftest.c
if test -z "`${CXX-cc} $1 -c conftest.c 2>&1`"; then
	ok=yes
else
		
  echo "**************************************************************"
  echo "* ERROR: The "'``'"inline'' keyword does not appear to work     *"
  echo "* properly on this system.  THIS WILL NOT COMPILE.             *"
  echo "* (It's >= 2001.  Get an ANSI C++ compiler.)                   *"
  case "${host_cpu}-${host_os}" in
    sparc-solaris2*) 
      echo "* Maybe you are using the wrong C++ compiler (try --with-gcc instead).  *"
      ;;
  esac
  echo "**************************************************************"
fi
rm -f conftest*
])
if test "ok" = yes; then
	:
	$3
else
	:
	$4
fi
])

AC_DEFUN(ACX_CHECK_CONST,
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(for const, ok,
[echo 'const int i=9;' > conftest.c
if test -z "`${CXX-cc} $1 -c conftest.c 2>&1`"; then
 echo "**************************************************************"
  echo "* ERROR: The "'``'"const'' keyword does not appear to work     *"
  echo "* properly on this system.  THIS WILL NOT COMPILE.             *"
  echo "* (It's >= 2001.  Get an ANSI C++ compiler.)                   *"
  case "${host_cpu}-${host_os}" in
    sparc-solaris2*) 
      echo "* Maybe you are using the wrong C++ compiler (try --with-gcc instead).  *"
      ;;
  esac
  echo "**************************************************************"
else
	ok=no
fi
rm -f conftest*
])
if test "$ok" = yes; then
	:
	$3
else
	:
	$4
fi
])


AC_DEFUN(ACX_TEST_GCC,
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(whether we are using gcc $1.$2 or later, ac_cv_prog_gcc_$1_$2,
[
dnl The semicolon after "yes" below is to pacify NeXT's syntax-checking cpp.
cat > conftest.c <<EOF
#ifdef __GNUC__
	yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX-cc} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  GCXX="yes"
else
  GCXX="no"
fi
])
if test "$GCXX" = yes; then
	:
	$3
else
	:
	$4
fi
])


AC_DEFUN(ACX_PROG_GCXX_VERSION,
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(whether we are using gcc $1.$2 or later, ac_cv_prog_gcc_$1_$2,
[
dnl The semicolon after "yes" below is to pacify NeXT's syntax-checking cpp.
cat > conftest.c <<EOF
#ifdef __GNUC__
#  if (__GNUC__ > $1) || (__GNUC__ == $1 && __GNUC_MINOR__ >= $2)
     yes;
#  endif
#endif
EOF
if AC_TRY_COMMAND(${CXX-cc} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  ac_cv_prog_gcc_$1_$2=yes
else
  ac_cv_prog_gcc_$1_$2=no
fi
])
if test "$ac_cv_prog_gcc_$1_$2" = yes; then
	:
	$3
else
	:
	$4
fi
])

AC_DEFUN(ACX_PROG_CXX_EGCS,
[ACX_PROG_GCXX_VERSION(2,95,acx_prog_egcs=yes,acx_prog_egcs=no)])

# Check to see if we are using a version of gcc that aligns the stack
# (true in gcc-2.95+, which have the -mpreferred-stack-boundary flag).
# Also check for stack alignment bug in gcc-2.95.x
# (see http://egcs.cygnus.com/ml/gcc-bugs/1999-11/msg00259.html), and
# whether main() is correctly aligned by the OS/libc/loader.
AC_DEFUN(ACX_GCXX_ALIGNS_STACK,
[
AC_REQUIRE([AC_PROG_CXX])
acx_gcc_aligns_stack=no
if test "$GCXX" = "yes"; then
ACX_CHECK_CXX_FLAGS(-mpreferred-stack-boundary=4, m_pref_stack_boundary_4)
if test "$ac_m_pref_stack_boundary_4" = "yes"; then
	AC_MSG_CHECKING([whether the stack is correctly aligned by gcc])
	save_CXXFLAGS="$CXXFLAGS"
	CXXFLAGS="-O -malign-double"
	AC_TRY_RUN([#include <stdlib.h>
#       include <stdio.h>
	struct yuck { int blechh; };
	int one(void) { return 1; }
	struct yuck ick(void) { struct yuck y; y.blechh = 3; return y; }
#       define CHK_ALIGN(x) if ((((long) &(x)) & 0x7)) { fprintf(stderr, "bad alignment of " #x "\n"); exit(1); }
	void blah(int foo) { double foobar; CHK_ALIGN(foobar); }
	int main(void) { double ok1; struct yuck y; double ok2; CHK_ALIGN(ok1);
                         CHK_ALIGN(ok2); y = ick(); blah(one()); return 0; }
	], [acx_gcc_aligns_stack=yes; acx_gcc_stack_align_bug=no], 
	acx_gcc_stack_align_bug=yes, acx_gcc_stack_align_bug=yes)
	CXXFLAGS="$save_CXXFLAGS"
	AC_MSG_RESULT($acx_gcc_aligns_stack)
fi
fi
if test "$acx_gcc_aligns_stack" = yes; then
	:
	$1
else
	:
	$2
fi
])


AC_DEFUN(ACX_PROG_CXX_MAXOPT,
[
AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([ACX_PROG_CXX_EGCS])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine "good" native compiler flags if none specified on command
# line
#if test "$ac_test_CXXFLAGS" != "set"; then
  CXXFLAGS=""
  case "${host_cpu}-${host_os}" in

  *linux*)
	echo "*       Congratulations! You are running linux.       *"; GCXX="yes"; 
	;;
  sparc-solaris2*) if test "$CXX" = cc; then
                    CXXFLAGS="-native -fast -xO5 -dalign"
                 fi;;

  alpha*-osf*)  if test "$CXX" = cc; then
                    CXXFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host -arch host -std1"
                fi;;

  hppa*-hpux*)  if test "$CXX" = cc; then
                    CXXFLAGS="-Ae +O3 +Oall"
                fi;;

   rs6000*-aix*)  if test "$CXX" = cc -o "$CXX" = xlc; then
                    CXXFLAGS="-O3 -qarch=pwrx -qtune=pwrx -qansialias -w"
                fi;;
   powerpc*-aix*)
	if test "$CXX" = cc -o "$CXX" = xlc; then
        	CXXFLAGS="-O3 -qarch=ppc -qansialias -w"
		echo "*******************************************************"
		echo "*  You have AIX on an unknown powerpc system.  It is  *"
		echo "*  recommended that you use                           *"
		echo "*                                                     *"
		echo "*   CXXFLAGS=-O3 -qarch=ppc -qtune=xxx -qansialias -w *"
		echo "*                                 ^^^                 *"
		echo "*  where xxx is 601, 603, 604, or whatever kind of    *"
                echo "*  PowerPC CPU you have.   For more info, man cc.     *"
		echo "*******************************************************"
        fi;;
  esac

  # use default flags for gcc on all systems
  if test $ac_cv_prog_gcc = yes; then
     CXXFLAGS="-O3 -fomit-frame-pointer -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
  fi

  # the egcs scheduler is too smart and destroys our own schedule.
  # Disable the first instruction scheduling pass.  The second
  # scheduling pass (after register reload) is ok.

#  if test "$acx_prog_egcs" = yes; then
#     CXXFLAGS="$CXXFLAGS -fno-schedule-insns -fschedule-insns2"
#  fi

  # test for gcc-specific flags:
  if test $ac_cv_prog_gcc = yes; then
    # -malign-double for x86 systems
    #ACX_CHECK_CXX_FLAGS(-malign-double,align_double,
#	CXXFLAGS="$CXXFLAGS -malign-double")
    # -fstrict-aliasing for gcc-2.95+
    ACX_CHECK_CXX_FLAGS(-fstrict-aliasing,fstrict_aliasing,
	CXXFLAGS="$CXXFLAGS -fstrict-aliasing")
 
    # loop unrolling
    ACX_CHECK_CXX_FLAGS(-funroll-loops,funrol_loops,
	CXXFLAGS="$CXXFLAGS -funroll-loops")
    # inlineing unrolling
    ACX_CHECK_CXX_FLAGS(-finline-functions,finline_functions,
	CXXFLAGS="$CXXFLAGS -finline-functions")
    # open up the template depth
    ACX_CHECK_CXX_FLAGS(-ftemplate-depth-40,ftemplate_depth_40,
	CXXFLAGS="$CXXFLAGS -ftemplate-depth-40")


  fi

  CPU_FLAGS=""
  if test "$GCXX" = "yes"; then
	  dnl try to guess correct CPU flags, at least for linux
	  case "${host_cpu}" in
	  i586*)  ACX_CHECK_CXX_FLAGS(-mcpu=pentium,cpu_pentium,
			[CPU_FLAGS=-mcpu=pentium],
			[ACX_CHECK_CXX_FLAGS(-mpentium,pentium,
				[CPU_FLAGS=-mpentium])])
		  ;;
	  i686*)  ACX_CHECK_CXX_FLAGS(-mcpu=pentiumpro,cpu_pentiumpro,
			[CPU_FLAGS=-mcpu=pentiumpro],
			[ACX_CHECK_CXX_FLAGS(-mpentiumpro,pentiumpro,
				[CPU_FLAGS=-mpentiumpro])])
		  ;;
	  powerpc*)
		cputype=`(grep cpu /proc/cpuinfo | head -1 | cut -d: -f2 | sed 's/ //g') 2> /dev/null`
		is60x=`echo $cputype | egrep "^60[0-9]e?$"`
		if test -n "$is60x"; then
			ACX_CHECK_CXX_FLAGS(-mcpu=$cputype,m_cpu_60x,
				CPU_FLAGS=-mcpu=$cputype)
		elif test "$cputype" = 750; then
                        ACX_PROG_GCXX_VERSION(2,95,
                                ACX_CHECK_CXX_FLAGS(-mcpu=750,m_cpu_750,
					CPU_FLAGS=-mcpu=750))
		fi
		if test -z "$CPU_FLAGS"; then
		        ACX_CHECK_CXX_FLAGS(-mcpu=powerpc,m_cpu_powerpc,
				CPU_FLAGS=-mcpu=powerpc)
		fi
		if test -z "$CPU_FLAGS"; then
			ACX_CHECK_CXX_FLAGS(-mpowerpc,m_powerpc,
				CPU_FLAGS=-mpowerpc)
		fi
	  esac
  fi

  if test -n "$CPU_FLAGS"; then
        CXXFLAGS="$CXXFLAGS $CPU_FLAGS"
  fi

  if test -z "$CXXFLAGS"; then
	echo ""
	echo "************************************************************"
        echo "* WARNING: Don't know the best CXXFLAGS for this system    *"
        echo "* Use  make CXXFLAGS=..., or edit the top level Makefile   *"
	echo "* (otherwise, a default of CXXFLAGS=-O3 will be used)      *"
	echo "************************************************************"
	echo ""
        CXXFLAGS="-O3"
  fi

  ACX_CHECK_CXX_FLAGS(${CXXFLAGS}, guessed_cflags, , [
	echo ""
        echo "**********************************************************"
        echo "* WARNING: The guessed CXXFLAGS don't seem to work with  *"
        echo "* your compiler.                                         *"
        echo "* Use  make CXXFLAGS=..., or edit the top level Makefile *"
        echo "**********************************************************"
        echo ""
        CXXFLAGS=""
  ])

#fi
])


dnl like AC_SUBST, but replace XXX_variable_XXX instead of @variable@
dnl This macro protects VARIABLE from being diverted twice
dnl if this macro is called twice for it.
dnl AC_SUBST(VARIABLE)
define(ACX_SUBST_XXX,
[ifdef([ACX_SUBST_XXX_$1], ,
[define([ACX_SUBST_XXX_$1], )dnl
AC_DIVERT_PUSH(AC_DIVERSION_SED)dnl
s=XXX_$1_XXX=[$]$1=g
AC_DIVERT_POP()dnl
])])



