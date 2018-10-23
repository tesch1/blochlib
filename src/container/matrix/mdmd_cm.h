

#ifndef _BLPHPIP_ZEMM_H_
#define _BLPHPIP_ZEMM_H_ 1


BEGIN_BL_NAMESPACE
#include "container/complex.h"

void mdmd_a1bc_md_cm(const int,const  int,const  int,const  complex*,const  complex*, complex*,
			const  int,const  int,const  int,const  complex);
void mdmd_acbc_md_cm(const int,const  int, const int,const  complex*,const  complex*, complex*,
			 const int,const  int,const int,const  complex,const  complex);

void mdmdt_a1bc_md_cm(const int, const int, const int,const  complex*, const complex*, complex*,
			  const int, const int, const int,const  complex);
void mdmdt_acbc_md_cm(const int,const  int,const  int,const  complex*, const complex*, complex*,
			 const  int,const  int,const  int, const complex, const complex);

void mdtmd_a1bc_md_cm(const int,const  int,const  int,const  complex*,const  complex*, complex*,
			  const int,const  int, const int, complex);
void mdtmd_acbc_md_cm(const int,const  int,const  int,const  complex*,const  complex*, complex*,
			  const int,const  int,const  int,const  complex,const  complex);

void mdtmdt_a1bc_md_cm(const int,const  int,const  int, const complex*,const  complex*, complex*,
			   const int,const  int,const  int, complex);
void mdtmdt_acbc_md_cm(const int,const  int, const int,const  complex*,const  complex*, complex*,
			  const  int,const  int,const  int,const  complex,const  complex);

END_BL_NAMESPACE


#endif
