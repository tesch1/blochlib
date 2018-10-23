
/* blassert.cc ********/


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
 	blassert.cc-->Assertion Class definitions

 	much thanks to 'Cheeta' (http://www.acl.lanl.gov/pooma/) for some helpful
 	hints here....
 */

#ifndef _bl_assert_cc_
#define _bl_assert_cc_ 1

#include "utils/blassert.h"


#ifndef ON_WINDOWS
 #include "blochconfig.h"
#endif

#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>

BEGIN_BL_NAMESPACE


std::string BL_exception::itost_bl(int i){
	char buffer[80];
	sprintf(buffer, "%d", i);
	return std::string(buffer);
} //std::string from an int

BL_exception::BL_exception(const char *file, const char * function, int line, std::string message)
{
	mess = 
	std::string("*** file: \"")+std::string(file)+
	std::string("\"\n*** line: ")+std::string(itost_bl(line))+
	std::string("\n*** function: \"")+ std::string(function)+"\"";
	if(message !="") mess+="\n"+message;
}


 void BL_exception::print(std::ostream &out)
{	out<<std::endl<<"Program Failure...quitting"<<std::endl<<mess<<std::endl<<std::endl;	}

/* Assertion class defs */

Assertion::Assertion(const char *msg, const char *file, int line)
{
  msg_m = new char[strlen(msg) + 1];
  strcpy(msg_m, msg);
  file_m = new char[strlen(file) + 1];
  strcpy(file_m, file);
  line_m = line;
}

Assertion::Assertion(const Assertion &a)
{
  msg_m = new char[strlen(a.what())+1];
  strcpy(msg_m, a.what());
  file_m = new char[strlen(a.file())+1];
  strcpy(file_m, a.file());
  line_m = a.line();
}

Assertion &Assertion::operator=(const Assertion &a)
{
  msg_m = new char[strlen(a.what())+1];
  strcpy(msg_m, a.what());
  file_m = new char[strlen(a.file())+1];
  strcpy(file_m, a.file());
  line_m = a.line();

  return *this;
}

/* the mighty barf function...*/

void barf(const char *msg, const char *file, int line ...)
{
  va_list ap;
  va_start(ap, line);
  char buf[256];
  vsprintf(buf, msg, ap);
  va_end(ap);

  Assertion a(buf, file, line);

  std::cerr<<"*** BlochLib Assertion Failure **"<<std::endl;
  std::cerr<<"*** "<< a.what() <<std::endl;
  std::cerr<<"*** File: "<<a.file() <<" Line: "<< a.line()<<std::endl;
#ifdef USE_EXCEPTIONS
	throw a;
#else
	#ifdef HAVE_SIGNAL_H
		raise(SIGINT);
	#else
		exit(0);
	#endif
#endif
}

END_BL_NAMESPACE


#endif


