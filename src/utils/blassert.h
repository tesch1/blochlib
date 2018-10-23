

/* blassert.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-25-01
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
 	blassert.h-->runtime assertion, and compile time assertions

 	much thanks to 'Cheeta' (http://www.acl.lanl.gov/pooma/) for some helpful
 	hints here....
 */

#ifndef _bl_assert_h_
#define _bl_assert_h_ 1

#include <iostream>
#include <string>

#include "blochconfig.h"

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif


BEGIN_BL_NAMESPACE


//The exception Class
class BL_exception
{
	private:
		std::string mess;
		std::string itost_bl(int i);
	public:
		BL_exception(){}
		BL_exception(std::string inmes):mess(inmes){}
		BL_exception(const char *file, const char *function, int line, std::string message="");
		void print(std::ostream &out=std::cout);
};

// the macro that decides if we are to use exceptions or not
#ifdef USE_EXCEPTIONS 
	#define BLEXCEPTION(C) throw BL_exception(__FILE__, __FUNCTION__, __LINE__ ,C);
#else
	#ifdef HAVE_SIGNAL_H
		#define BLEXCEPTION(C) \
			std::cerr<<std::endl<<" Error: in file: "<<__FILE__<<"\n at line: "<<__LINE__<<"\n  function: "<<__FUNCTION__<<std::endl; \
			std::cerr<<" C "<<std::endl; \
			raise(SIGINT);
	#else
		#define BLEXCEPTION(C) \
			std::cerr<<std::endl<<" Error: in file: "<<__FILE__<<"\n at line: "<<__LINE__<<"\n  function: "<<__FUNCTION__<<std::endl;  \
			std::cerr<<" C"<<std::endl; \
			exit(0);
	#endif
#endif

/*
 Classes:
 RunTimeAssert macro, RunTimeInsist macro, CompTimeAssert struct
*/

/*
 RunTimeAssert(bool c) is a Run Time assertion macro. (use these for testing phases checking)
 RunTimeInsist(bool c) is a run-time assertion macro. (use this for things the MUST happen)
 CompTimeAssert(bool c, char *m) use this for compile time checking
*/

#ifndef ON_WINDOWS
 #include "blochconfig.h"
#endif

//This class here is quite standard on most libs i have ever seen...

/* : std::runtime_error */
class Assertion
{
  char *msg_m;
  char *file_m;
  int line_m;
public:
  Assertion(const char *msg, const char *file, int line);
  Assertion(const Assertion &a);
  ~Assertion() { delete[] msg_m; delete [] file_m; }
  Assertion &operator=(const Assertion &a);
  const char *what() const { return msg_m; }
  const char *file() const { return file_m; }
  int line() const { return line_m; }

  template<class OS>
  void print(OS &os) const
  {
    os << "*** BlochLib Assertion Failure ***\n";
    os << "*** " << what() << "\n";
    os << "*** File " << file() << "; Line " << line();
  }
};

//The main death function...'barf'( i.e. throw up, i.e. toss thy cookies, i.e. die of infulenza)

void barf(const char *msg, const char *file, int line ...);

/*
 CompTimeAssert: compile time assert.

 It tests the condition at compile time and if it is false
 you get a compile error that it can't find ComTimeAssert<false>::test().

 If NOCTASSERT is defined, CompTimeAssert will revert to the equivalent of
 RunTimeAssert. To turn off the test completely, define NORTASSERT as well.
*/

#if defined(NOCTASSERT)
 #if defined(NORTASSERT)
  #define CompTimeAssert(c)
 #else
  #define CompTimeAssert(c) if (c) {} else barf(#c, __FILE__, __LINE__)
 #endif
#else
 #define CompTimeAssert(c) CompTimeAssertStruct<(c)>::test()
  template<bool B> struct CompTimeAssertStruct {};
  template<> struct CompTimeAssertStruct<true> { static void test() {} };
#endif


/*
 RunTimeAssert: Run-time assertion


 The RunTimeAssert macro is intended to be used for validating preconditions
 which SHOULD be true in order for following code to be correct, etc.

 Things like bounds checking on vector lengths, and such can take lots
 of processor time, and are not nessesary if the the code has been
 written properly, so this macro should be used for debugging and
 testing phases

 or NORTASSERT is defined, RunTimeAssert will just be an empty macro.

*/

#if defined NORTASSERT
 #define RunTimeAssert(c)
#else
 #define RunTimeAssert(c) if (c) {} else barf(#c, __FILE__, __LINE__)
#endif


/*
 RunTimeInsist:

 this macro is for ALL errors that MUST happen regardless of debugging phases
 of testing...things like file opening MUST happen before anything else can continue
 THere is no possible way to turn this off

*/

#define RunTimeInsist(c,m) if (c) {} else ::barf(m, __FILE__, __LINE__)
#define RunTimeInsist1(c,m,a1)                                                    \
  if (c) {} else barf(m,__FILE__,__LINE__,a1)
#define RunTimeInsist2(c,m,a1,a2)                                                 \
  if (c) {} else barf(m,__FILE__,__LINE__,a1,a2)
#define RunTimeInsist3(c,m,a1,a2,a3)                                              \
  if (c) {} else barf(m,__FILE__,__LINE__,a1,a2,a3)
#define RunTimeInsist4(c,m,a1,a2,a3,a4)                                           \
  if (c) {} else barf(m,__FILE__,__LINE__,a1,a2,a3,a4)
#define RunTimeInsist5(c,m,a1,a2,a3,a4,a5)                                        \
  if (c) {} else barf(m,__FILE__,__LINE__,a1,a2,a3,a4,a5)
#define RunTimeInsist6(c,m,a1,a2,a3,a4,a5,a6)                                     \
  if (c) {} else barf(m,__FILE__,__LINE__,a1,a2,a3,a4,a5,a6)
#define RunTimeInsist7(c,m,a1,a2,a3,a4,a5,a6,a7)                                  \
  if (c) {} else barf(m,__FILE__,__LINE__,a1,a2,a3,a4,a5,a6,a7)



END_BL_NAMESPACE


#endif

