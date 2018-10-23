/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-25-01
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
	space_ten.cc--> spacial spherical tensors up to rank 2
*/

#ifndef _space_ten_cc_
#define _space_ten_cc_

#include "QMspins/space_ten.h"
#include "container/rankType.h"
#include "utils/blassert.h"

BEGIN_BL_NAMESPACE


complex A1(double theta, double phi, int m){
  if(m!=0&&m!=1&&m!=-1){
    BLEXCEPTION(std::string("rank one SPACE tensors cannot have m,")+itost(m)+", like that!")
    return 0;
  }else{
    if(m==1){
      return -pow(3./(8.),0.5)*sin(theta)*exp(complexi*phi);
    }else if(m==0){
      return pow(3./(4.), 0.5)*cos(theta);
    }else if(m==-1.){
      return pow(3./(8.),0.5)*sin(theta)*exp(-complexi*phi);
    }
  }
  return 0;
}

complex A2(double theta, double phi, int m){
	if(m!=0&&m!=1&&m!=-1&&m!=-2&&m!=2){
		BLEXCEPTION(std::string("rank two SPACE tensors cannot have m,")+itost(m)+", like that!")
    return 0;
	}else{
		if(m==-2){
			return pow(3./8.,0.5)*sin(theta)*sin(theta)*exp(2.*complexi*phi);
		}else if(m==-1){
			return pow(3./8.,0.5)*sin(2.*theta)*exp(complexi*phi);
		}else if(m==0){
			return 1./2.*(3.*cos(theta)*cos(theta)-1.);
		}else if(m==1){
			return -pow(3./8.,0.5)*sin(2.*theta)*exp(-complexi*phi);
		}else if(m==2){
			return pow(3./8.,0.5)*sin(theta)*sin(theta)*exp(-2.*complexi*phi);
		}
  }
  return 0;
}

complex wigner1(double alpha, double beta,double gamma,int m, int mp){
	if(m!=0&&m!=1&&m!=-1&&mp!=0&&mp!=1&&mp!=-1){
		BLEXCEPTION(std::string("rank one WIGNERS Rotations tensors cannot have m,")+itost(m)+", like that!")
		return 0;
	}else{
		if(mp==-1){
			if(m==1){
				return pow(sin(beta/2.), 2.)*exp(-complexi*alpha)*exp(complexi*gamma);
			}else if(m==0){
				return -pow(1./2.,0.5)*sin(beta)*exp(complexi*gamma);
			}else if(m==-1.){
				return pow(cos(beta/2.), 2.)*exp(complexi*alpha)*exp(complexi*gamma);
			}
		}else if(mp==0){
			if(m==1.){
				return -pow(1./2.,0.5)*sin(beta)*exp(-complexi*alpha);
			}else if(m==0){
				return cos(beta);
			}else if(m==-1){
				return pow(1./2.,0.5)*sin(beta)*exp(complexi*alpha);
			}
		}else if(mp==1){
			if(m==1){
				return pow(cos(beta/2.), 2.)*exp(-complexi*alpha)*exp(-complexi*gamma);
			}else if(m==0){
				return pow(1./2.,0.5)*sin(beta)*exp(-complexi*gamma);
			}else if(m==-1){
				return pow(sin(beta/2.), 2.)*exp(complexi*alpha)*exp(-complexi*gamma);
			}
		}
	}
  return 0;
}

complex wigner2(double alpha, double beta, double gamma,int m, int n){
  if(m!=0&&m!=1&&m!=-1&&m!=-2&&m!=2&&n!=0&&n!=1&&n!=-1&&n!=-2&&n!=2){
    BLEXCEPTION(std::string("rank two WIGNERS Rotations tensors cannot have m,")+itost(m)+", like that!")		
    return 0;
  }else{
	double cb=cos(beta);
	double sb=sin(beta);
    if(n==2){
      if(m==-2){
				return 1./4.*(1.-cb)*(1.-cb)*exp(-complexi*(-2.*alpha+2.*gamma));
      }else if(m==-1){
				return 1./2.*(1.-cb)*sb*exp(-complexi*(-alpha+2.*gamma));
      }else if(m==0){
				return pow(3./8., 0.5)*sb*sb*exp(-2.*complexi*gamma);
      }else if(m==1){
				return 1./2.*(1.+cb)*sb*exp(-complexi*(alpha+2.*gamma));
      }else if(m==2){
				return 1./4.*(1.+cb)*(1.+cb)*exp(-complexi*(2.*alpha+2.*gamma));
      }
    }else if(n==1){
      if(m==-2){
				return 1./2.*sb*(1-cb)*exp(-complexi*(-2.*alpha+gamma));
      }else if(m==-1.){
				return (1./2.*(1.+cb)-cb*cb)*exp(-complexi*(-alpha+gamma));
      }else if(m==0){
				return pow(3./8., 0.5)*sin(2.*beta)*exp(-complexi*gamma);
      }else if(m==1){
				return (1./2.*(cb-1.)+cb*cb)*exp(-complexi*(alpha+gamma));
      }else if(m==2){
				return -1./2.*sb*(1+cb)*exp(-complexi*(2.*alpha+gamma));
      }
    }else if(n==0){
      if(m==-2){
				return pow(3./8., 0.5)*sb*sb*exp(2.*complexi*alpha);
      }else if(m==-1){
				return pow(3./8., 0.5)*sin(2.*beta)*exp(complexi*alpha);
      }else if(m==0){
				return 1./2.*(3.*cb*cb-1.);
      }else if(m==1){
				return -pow(3./8., 0.5)*sin(2.*beta)*exp(-complexi*alpha);
      }else if(m==2){
				return pow(3./8., 0.5)*sb*sb*exp(-2.*complexi*alpha);
      }
    }else if(n==-1){
      if(m==-2){
				return 1./2.*(1+cb)*sb*exp(-complexi*(-2.*alpha-gamma));
      }else if(m==-1.){
				return (1./2.*(cb-1)+cb*cb)*exp(-complexi*(-alpha-gamma));
      }else if(m==0){
				return -pow(3./8., 0.5)*sin(2.*beta)*exp(complexi*gamma);
      }else if(m==1){
				return (1./2.*(cb+1)-cb*cb)*exp(-complexi*(alpha-gamma));
      }else if(m==2){
				return -1./2.*(1-cb)*sb*exp(-complexi*(2.*alpha-gamma));
      }
    }else if(n==-2){
      if(m==-2){
				return 1./4.*(1.+cb)*(1.+cb)*exp(-complexi*(-2.*alpha-2.*gamma));
      }else if(m==-1){
				return -1./2.*(1.+cb)*sb*exp(-complexi*(-alpha-2.*gamma));
      }else if(m==0){
				return pow(3./8., 0.5)*sb*sb*exp(2.*complexi*gamma);
      }else if(m==1){
				return -1./2.*(1.-cb)*sb*exp(-complexi*(alpha-2.*gamma));
      }else if(m==2){
				return 1./4.*(1.-cb)*(1.-cb)*exp(-complexi*(2.*alpha-2.*gamma));
      }
    }
  }
  return 0;
}


double dnm2(double beta, int n, int m){
   switch(n){
   case -2:
      switch(m){
      case -2:
         return  pow(cos(beta/2.0), 4.0);
      case -1:
         return  (1.0+cos(beta))*sin(beta)/2.0;
      case 0:
         return  pow(3./8., 0.5)*pow(sin(beta), 2.0);
      case 1:
         return  2.0*cos(beta/2.0)*pow(sin(beta/2.0), 3.0);
      case 2:
         return  pow(sin(beta/2.0), 4.0);
      default:

       BLEXCEPTION(std::string(" invalid m")+itost(m))
      }
   case -1:
      switch(m){
      case -2:
         return  -2.0*pow(cos(beta/2.0), 3.0)*sin(beta/2.0);
      case -1:
         return  (cos(beta)+cos(2.0*beta))/2.0;
      case 0:
         return  pow(3.0/2.0, 0.5)*cos(beta)*sin(beta);
      case 1:
         return  (cos(beta)-cos(2.0*beta))/2.0;
      case 2:
         return  2.0*cos(beta/2.0)*pow(sin(beta/2.0), 3.0);
      default:
		BLEXCEPTION(std::string("invalid m ")+itost(m))
      }
   case 0:
      switch(m){
      case -2:
         return  pow(3./8., 0.5)*pow(sin(beta), 2.0);
      case -1:
         return  -pow(3.0/2.0, 0.5)*cos(beta)*sin(beta);
      case 0:
         return  (1.0+3.0*cos(2.0*beta))/4.0;
      case 1:
         return  pow(3.0/2., 0.5)*cos(beta)*sin(beta);
      case 2:
         return  pow(3./8., 0.5)*pow(sin(beta), 2.0);
      default:

		BLEXCEPTION(std::string("invalid m ")+itost(m))

      }
   case 1:
      switch(m){
      case -2:
         return  (sin(2.0*beta)-2.0*sin(beta))/4.0;
      case -1:
         return  (cos(beta)-cos(2.0*beta))/2.0;
      case 0:
         return  -pow(3.0/2., 0.5)*cos(beta)*sin(beta);
      case 1:
         return  (cos(beta)+cos(2.0*beta))/2.0;
      case 2:
         return  (1.0+cos(beta))*sin(beta)/2.0;
      default:

 		BLEXCEPTION(std::string("invalid m ")+itost(m))
     }
   case 2:
      switch(m){
      case -2:
         return  pow(sin(beta/2.0), 4.0);
      case -1:
         return  (sin(2.0*beta)-2.0*sin(beta))/4.0;
      case 0:
         return  pow(3./8., 0.5)*pow(sin(beta), 2.0);
      case 1:
         return  -2.0*pow(cos(beta/2.0), 3.0)*sin(beta/2.0);
      case 2:
         return  pow(cos(beta/2.0), 4.0);
      default:

		BLEXCEPTION(std::string("invalid m ")+itost(m))
      }
   default:

		BLEXCEPTION("invalid n "+itost(n))
    	return 0;
   }
}

double dnm4(double beta, int n, int m){
	   switch(n){
	   case -4:
	      switch(m){
	      case -4:
	         return  pow(cos(beta/2.0),8.0);
	      case -3:
	         return  pow(8.0, 0.5)*pow(cos(beta/2.0), 7.0)*sin(beta/2.0);
	      case -2:
	         return  pow(28.0, 0.5)*pow(cos(beta/2.0), 6.0)*pow(sin(beta/2.0), 2.0);
	      case -1:
	         return  pow(56.0, 0.5)*pow(cos(beta/2.0),5.0)*pow(sin(beta/2.0), 3.0);
	      case 0:
	         return  pow(35.0/128.0, 0.5)*pow(sin(beta), 4.0);
	      case 1:
	         return  pow(56.0, 0.5)*pow(cos(beta/2.0), 3.0)*pow(sin(beta/2.0), 5.0);
	      case 2:
	         return  pow(28.0, 0.5)*pow(cos(beta/2.0), 2.0)*pow(sin(beta/2.0), 6.0);
	      case 3:
	         return  pow(8.0, 0.5)*cos(beta/2.0)*pow(sin(beta/2.0), 7.0);
	      case 4:
	         return  pow(sin(beta/2.0), 8.0);
	      default:

			BLEXCEPTION(std::string("invalid m ")+itost(m))
	         return 0.0;
	      }
	   case -3:
	      switch(m){
	      case -4:
	         return  -pow(8.0,0.5)*pow(cos(beta/2.0), 7.0)*sin(beta/2.0);
	      case -3:
	         return  pow(cos(beta/2.0), 6.0)*(4.0*cos(beta)-3.0);
	      case -2:
	         return  pow(14.0,0.5)*pow(cos(beta/2.0), 5.0)*(sin(3.0*beta/2.0)-2.0*sin(beta/2.0));
	      case -1:
	         return  pow(7.0,0.5)*pow(cos(beta/2.0), 4.0)*(4.0*cos(beta)-1.0)*pow(sin(beta/2.0), 2.0);
	      case 0:
	         return  pow(35.0/16.0,0.5)*cos(beta)*pow(sin(beta), 3.0);
	      case 1:
	         return  pow(7.0,0.5)*pow(cos(beta/2.0), 2.0)*(4.0*cos(beta)+1.0)*pow(sin(beta/2.0), 4.0);
	      case 2:
	         return  pow(14.0,0.5)*(2.0*cos(beta/2.0)+cos(3.0*beta/2.0))*pow(sin(beta/2.0), 5.0);
	      case 3:
	         return  (3.0+4.0*cos(beta))*pow(sin(beta/2.0), 6.0);
	      case 4:
	         return  pow(8.0,0.5)*cos(beta/2.0)*pow(sin(beta/2.0), 7.0);
	      default:

	         BLEXCEPTION(std::string("invalid m ")+itost(m))
			 return 0.0;
	      }
	   case -2:
	      switch(m){
	      case -4:
	         return  pow(28.0, 0.5)*pow(cos(beta/2.0), 6.0)*pow(sin(beta/2.0), 2.0);
	      case -3:
	         return  -pow(14.0, 0.5)*pow(cos(beta/2.0), 5.0)*(sin(3.0*beta/2.0)-2.0*sin(beta/2.0));
	      case -2:
	         return  pow(cos(beta/2.0),4.0)*(9-14.0*cos(beta)+7.0*cos(2.0*beta))/2.0;
	      case -1:
	         return  (3.0*sin(beta)-2.0*sin(2.0*beta)+7.0*(sin(3.0*beta)+sin(4.0*beta)))/pow(512.0,0.5);
	      case 0:
	         return  pow(5.0/128.0, 0.5)*(5+7.0*cos(2.0*beta))*pow(sin(beta), 2.0);
	      case 1:
	         return  (3.0*sin(beta)+2.0*sin(2.0*beta)+7.0*(sin(3.0*beta)-sin(4.0*beta)))/pow(512.0, 0.5);
	      case 2:
	         return  pow(sin(beta/2.0), 4.0)*(9.0+14.0*cos(beta)+7.0*cos(2.0*beta))/2.0;
	      case 3:
	         return  pow(14.0, 0.5)*(2.0*cos(beta/2.0)+cos(3.0*beta/2.0))*pow(sin(beta/2.0), 5.0);
	      case 4:
	         return  pow(28.0, 0.5)*pow(cos(beta/2.0), 2.0)*pow(sin(beta/2.0),6.0);
	      default:

	         BLEXCEPTION(std::string("invalid m ")+itost(m))
			 return 0.0;
	      }
	   case -1:
	      switch(m){
	      case -4:
	         return  -pow(56.0, 5.0)*pow(cos(beta/2.0),5.0)*pow(sin(beta/2.0), 3.0);
	      case -3:
	         return  pow(7.0, 0.5)*pow(cos(beta/2.0),4.0)*(4.0*cos(beta)-1)*pow(sin(beta/2.0), 2.0);
	      case -2:
	         return  (-3.0*sin(beta)+2.0*sin(2.0*beta)-7.0*(sin(3.0*beta)+sin(4.0*beta)))/pow(512.0, 0.5);
	      case -1:
	         return  (9.0*cos(beta)+2.0*cos(2.0*beta)+7.0*(cos(3.0*beta)+2.0*cos(4.0*beta)))/32.0;
	      case 0:
	         return  pow(5.0/1024.0, 0.5)*(2.0*sin(2.0*beta)+7.0*sin(4.0*beta));
	      case 1:
	         return  (9.0*cos(beta)-2.0*cos(2.0*beta)+7.0*(cos(3.0*beta)-2.0*cos(4.0*beta)))/32.0;
	      case 2:
	         return  (3.0*sin(beta)+2.0*sin(2.0*beta)+7.0*(sin(3.0*beta)-sin(4.0*beta)))/pow(512.0, 0.5);
	      case 3:
	         return  pow(7.0, 0.5)*pow(cos(beta/2.0), 2.0)*(4.0*cos(beta)+1)*pow(sin(beta/2.0), 4.0);
	      case 4:
	         return  pow(56.0, 0.5)*pow(cos(beta/2.0),3.0)*pow(sin(beta/2.0), 5.0);
	      default:

	         BLEXCEPTION(std::string("invalid m ")+itost(m))
			 return 0.0;
	      }
	   case 0:
	      switch(m){
	      case -4:
	         return  pow(35./128., 0.5)*pow(sin(beta),4.0);
	      case -3:
	         return  -pow(35./16., 0.5)*cos(beta)*pow(sin(beta),3.0);
	      case -2:
	         return  pow(5./128., 0.5)*(5.0+7.0*cos(2.0*beta))*pow(sin(beta),2.0);
	      case -1:
	         return  -pow(5./1024., 0.5)*(2.0*sin(2.0*beta)+7.0*sin(4.0*beta));
	      case 0:
	         return  (9.0+20.0*cos(2.0*beta)+35*cos(4.0*beta))/64.0;
	      case 1:
	         return  pow(5./1024., 0.5)*(2.0*sin(2.0*beta)+7.0*sin(4.0*beta));
	      case 2:
	         return  pow(5./128., 0.5)*(5.0+7.0*cos(2.0*beta))*pow(sin(beta),2.0);
	      case 3:
	         return  pow(35./16., 0.5)*cos(beta)*pow(sin(beta), 3.0);
	      case 4:
	         return  pow(35./128., 0.5)*pow(sin(beta),4.0);
	      default:

			 BLEXCEPTION(std::string("invalid m ")+itost(m))
			 return 0.0;
	      }
	   case 1:
	      switch(m){
	      case -4:
	         return  -pow(56.0, 0.5)*pow(cos(beta/2.0),3.0)*pow(sin(beta/2.0),5.0);
	      case -3:
	         return  pow(7.0, 0.5)*pow(cos(beta/2.0),2.0)*(4.0*cos(beta)+1)*pow(sin(beta/2.0),4.0);
	      case -2:
	         return  -(3.0*sin(beta)+2.0*sin(2.0*beta)+7.0*(sin(3.0*beta)-sin(4.0*beta)))/pow(512.0, 0.5);
	      case -1:
	         return  (9.0*cos(beta)-2.0*cos(2.0*beta)+7.0*(cos(3.0*beta)-2.0*cos(4.0*beta)))/32.0;
	      case 0:
	         return  -pow(5.0/1024.0, 0.5)*(2.0*sin(2.0*beta)+7.0*sin(4.0*beta));
	      case 1:
	         return  (9.0*cos(beta)+2.0*cos(2.0*beta)+7.0*(cos(3.0*beta)+2.0*cos(4.0*beta)))/32.0;
	      case 2:
	         return  (3.0*sin(beta)-2.0*sin(2.0*beta)+7.0*(sin(3.0*beta)+sin(4.0*beta)))/pow(512.0, 0.5);
	      case 3:
	         return  pow(7.0, 0.5)*pow(cos(beta/2.0),4.0)*(4.0*cos(beta)-1)*sin(beta/2.0)*sin(beta/2.0);
	      case 4:
	         return  pow(56.0, 0.5)*pow(cos(beta/2.0),5.0)*pow(sin(beta/2.0),3.0);
	      default:

	         BLEXCEPTION(std::string("invalid m ")+itost(m))
			 return 0.0;
	      }
	   case 2:
	      switch(m){
	      case -4:
	         return  pow(28.0,0.5)*pow(cos(beta/2.0),2.0)*pow(sin(beta/2.0), 6.0);
	      case -3:
	         return  -pow(14.0,0.5)*(2.0*cos(beta/2.0)+cos(3.0*beta/2.0))*pow(sin(beta/2.0),5.0);
	      case -2:
	         return  pow(sin(beta/2.0),4.0)*(9+14.0*cos(beta)+7.0*cos(2.0*beta))/2.0;
	      case -1:
	         return  -(3.0*sin(beta)+2.0*sin(2.0*beta)+7.0*(sin(3.0*beta)-sin(4.0*beta)))/pow(512.0,0.5);
	      case 0:
	         return  pow(5.0/128.0, 0.5)*(5.0+7.0*cos(2.0*beta))*pow(sin(beta), 2.0);
	      case 1:
	         return  (-3.0*sin(beta)+2.0*sin(2.0*beta)-7.0*(sin(3.0*beta)+sin(4.0*beta)))/pow(512.0, 0.5);
	      case 2:
	         return  pow(cos(beta/2.0), 4.0)*(9-14.0*cos(beta)+7.0*cos(2.0*beta))/2.0;
	      case 3:
	         return  pow(14.0, 0.5)*pow(cos(beta/2.0), 5.0)*(sin(3.0*beta/2.0)-2.0*sin(beta/2.0));
	      case 4:
	         return  pow(28.0, 0.5)*pow(cos(beta/2.0), 6.0)*pow(sin(beta/2.0),2.0);
	      default:
	         BLEXCEPTION(std::string("invalid m ")+itost(m))
			 return 0.0;
	      }
	   case 3:
	      switch(m){
	      case -4:
	         return  -pow(8.0, 0.5)*cos(beta/2.0)*pow(sin(beta/2.0), 7.0);
	      case -3:
	         return  pow(sin(beta/2.0), 6.0)*(4.0*cos(beta)+3.0);
	      case -2:
	         return  -pow(14.0, 0.5)*(cos(3.0*beta/2.0)+2.0*cos(beta/2.0))*pow(sin(beta/2.0), 5.0);
	      case -1:
	         return  pow(7.0, 0.5)*pow(cos(beta/2.0), 2.0)*(4.0*cos(beta)+1)*pow(sin(beta/2.0),4.0);
	      case 0:
	         return  -pow(35.0/16.0, 0.5)*cos(beta)*pow(sin(beta), 3.0);
	      case 1:
	         return  pow(7.0, 0.5)*pow(cos(beta/2.0), 4.0)*(4.0*cos(beta)-1)*pow(sin(beta/2.0), 2.0);
	      case 2:
	         return  -pow(14.0, 0.5)*(sin(3.0*beta/2.0)-2.0*sin(beta/2.0))*pow(cos(beta/2.0), 5.0);
	      case 3:
	         return  (4.0*cos(beta)-3.0)*pow(cos(beta/2.0), 6.0);
	      case 4:
	         return  pow(8.0,0.5)*pow(cos(beta/2.0), 7.0)*sin(beta/2.0);
	      default:
	         BLEXCEPTION(std::string("invalid m ")+itost(m))
			 return 0.0;
	      }
	   case 4:
	      switch(m){
	      case -4:
	         return  pow(sin(beta/2.0), 8.0);
	      case -3:
	         return  -pow(8.0,0.5)*cos(beta/2.0)*pow(sin(beta/2.0), 7.0);
	      case -2:
	         return  pow(28.0,0.5)*pow(cos(beta/2.0), 2.0)*pow(sin(beta/2.0), 6.0);
	      case -1:
	         return  -pow(56.0,0.5)*pow(cos(beta/2.0), 3.0)*pow(sin(beta/2.0),5.0);
	      case 0:
	         return  pow(35.0/128.0,0.5)*pow(sin(beta), 4.0);
	      case 1:
	         return  -pow(56.0,0.5)*pow(cos(beta/2.0), 5.0)*pow(sin(beta/2.0), 3.0);
	      case 2:
	         return  pow(28.0,0.5)*pow(cos(beta/2.0),6.0)*pow(sin(beta/2.0), 2.0);
	      case 3:
	         return  -pow(8.0,0.5)*pow(cos(beta/2.0), 7.0)*sin(beta/2.0);
	      case 4:
	         return  pow(cos(beta/2.0),8.0);
	      default:
	         BLEXCEPTION(std::string("invalid m ")+itost(m))
			 return 0.0;
	      }
	   default:
	      BLEXCEPTION(std::string("invalid n ")+itost(n))
		  return 0.0;
   }
}

double dnml(double beta, int n, int m, int l){
	switch(l){
		case 2: return dnm2(beta, n, m);
		case 4: return dnm4(beta, n, m);
		default:
	      BLEXCEPTION(std::string("invalid l (so far only 2 and 4 allowed) ")+itost(l))
		  return 0.0;
	}
}


//gives the values of the Quadropolar hamiltonian
//in the PAS frame: only rank 2 matters
//note: the basic coupling is added during the hamiltonian
//set up as it requires the spin quantum number and
//that is not know until we get our data from the file
//it just multiplies the entire thing by a constant
//         [ 5 ]1/2                       [  5  ]1/2
//   A   = |---|    A    =  0      A     =|-----]*eta
//    2,0  [4Pi]     2,+/-1         2,+/-2[24*pi]
//
//
double A_Q_pas(double eta, int m){
  if(m!=0&&m!=1&&m!=-1&&m!=-2&&m!=2){
    BLEXCEPTION(std::string("rank two QUADRAPOLE Rotations tensors cannot m's, ")+itost(m)+", like that!")
    return 0;
  }else{
    if(m==2||m==-2){
      return 1./2.*eta;
    }else if(m==1||m==-1){
      return 0.;
    }else if(m==0){
      return 1./3.;
    }
  }
  return 0;
}

//CSA spacial tensors in the PAS frame...
//i am negelcting the usual Rank 1 compenents as they
//are asymetric and introduce little after powder averaging
//
//           1/2
//  A  =  -[3]  iso
//   0,0
//       [ 3 ]1/2
//  A   =|---| *del     A =  0     A  =   1/2 del*eta
//   2,0 [ 2 ]           2,+/-1     2,+/-2
//
//


double A_CSA_pas(double iso, double del, double eta,int l, int m){
  if(l==0&&m==0){
    return -pow(3., 0.5)*iso;
  }else if(l==1){
    return 0;
  }else if(l==2){
    if(m==-2||m==2){
      return 1./2.*del*eta;
    }else if(m==-1||m==1){
      return 0;
    }else if(m==0){
      return del;
    }
  }
  return 0;
}

//Dipole spacial tensors in the PAS frame
//all are rank 2
double A_Dip_pas(double dip, int m){
  if(m==0){
    return (-1.)*pow(6.,0.5)*dip;
  }else{
    return 0;
  }
  return 0;
}

rmatrix rotationMatrix3D(double a, double b, double g)
{
	rmatrix rot(3,3);
	double cg=cos(g), cb=cos(b), ca=cos(a), sg=sin(g), sb=sin(b), sa=sin(a);
	/*rot(0,0)=cos(a)*cos(b)*cos(g)-sin(a)*sin(g);
	rot(0,1)=sin(a)*cos(b)*cos(g)+cos(a)*sin(g);
	rot(0,2)=-sin(b)*cos(g);

	rot(1,0)=-cos(a)*cos(b)*sin(g)-sin(a)*cos(g);
	rot(1,1)=-sin(a)*cos(b)*sin(g)+cos(a)*cos(g);
	rot(1,2)=sin(b)*sin(g);

	rot(2,0)= cos(a)*sin(b);
	rot(2,1)= sin(a)*sin(b);
	rot(2,2)=  cos(b);*/

	rot(0,0)=ca*cb*cg-sa*sg;
	rot(0,1)=sa*cb*cg+ca*sg;
	rot(0,2)=-sb*cg;

	rot(1,0)=-ca*cb*sg-sa*cg;
	rot(1,1)=-sa*cb*sg+ca*cg;
	rot(1,2)=sb*sg;

	rot(2,0)= ca*sb;
	rot(2,1)= sa*sb;
	rot(2,2)=  cb;
	/*rot(0,0)=ca*cg-cb*sa*sg;
	rot(0,1)=cb*cg*sa+ca*sg;
	rot(0,2)=sa*sb;

	rot(1,0)=-cg*sa-ca*cb*sg;
	rot(1,1)=ca*cb*cg-sa*sg;
	rot(1,2)=ca*sb;

	rot(2,0)= sb*sg;
	rot(2,1)= -cg*sb;
	rot(2,2)=  cb;*/

	return rot;
}

rmatrix rotationMatrix3D(double phi, double theta)
{
	rmatrix rot(3,3);

	double ct=cos(theta), cp=cos(phi), sp=sin(phi), st=sin(theta);
	/*rot(0,0)=cp;
	rot(0,1)=ct*sp;
	rot(0,2)=sp*st;


	rot(1,0)=-sp;
	rot(1,1)=cp*ct;
	rot(1,2)=cp*st;

	rot(2,0)=0.0;
	rot(2,1)=-st;
	rot(2,2)=ct;*/

	rot(0,0)=cp;
	rot(0,1)=sp;
	rot(0,2)=0.0;


	rot(1,0)=-ct*sp;
	rot(1,1)=cp*ct;
	rot(1,2)=st;

	rot(2,0)=sp*st;
	rot(2,1)=-cp*st;
	rot(2,2)=ct;

	return rot;
}

Rotations::Rotations(){
	alpha=-11415; theta=-11415; beta=-11415; phi=-11415; gamma=-11415, chi=-11415;
	powderAs.resize(5,5,ZeroType<complex>::zero());
	spinnerAs.resize(5,0);
	spinnerAsMat.resize(5,5,0);
	cartPowder.resize(3,3);
	cartSpin.resize(3,3);
	RotationType=All;
}

Rotations::Rotations(int inty){
	alpha=-11415; theta=-11415; beta=-11415; phi=-11415; gamma=-11415, chi=-11415;
	powderAs.resize(5,5,ZeroType<complex>::zero());
	spinnerAs.resize(5,0);
	spinnerAsMat.resize(5,5,0);
	cartPowder.resize(3,3);
	cartSpin.resize(3,3);
	RotationType=inty;
	rotate=false;
}

Rotations::Rotations(matrix &inP, Vector<complex> &inA){
	powderAs=inP;
	spinnerAs=inA;
	RotationType=Spherical;
	rotate=false;
}

Rotations::Rotations(const Rotations &rhs)
{
	alpha=rhs.alpha;
	beta=rhs.beta;
	gamma=rhs.gamma;
	theta=rhs.theta;
	phi=rhs.phi;
	chi=rhs.chi;
	spinnerAs=rhs.spinnerAs;
	powderAs=rhs.powderAs;
	spinnerAsMat=rhs.spinnerAsMat;
	cartPowder=rhs.cartPowder;
	cartSpin=rhs.cartSpin;
	RotationType=rhs.RotationType;
	rotate=rhs.rotate;
}

Rotations::Rotations(double inp, double inth, double ing, double ina, double inb)
{
	setPowderAngles(inp, inth, ing);
	setRotorAngles(ina, inb);
}

complex Rotations::getPow(int i, int j){
	return powderAs(i,j);
}

complex Rotations::getSpn(int i){
	return spinnerAs[i];
}


void Rotations::setPowderAs(double pin, double thin, double gin)
{	setPowderAngles( pin,  thin,  gin);	}

void Rotations::setPowderAngles(double pin, double thin, double gin){
	if(pin != phi || thin != theta || gin !=gamma)
	{
		phi=pin;
		theta=thin;
		gamma=gin;
		rotate=true;
		int m=-2, mp=-2,i=0,j=0;
		switch(RotationType){
			case Spherical:
				for(i=0,m=-2;m<=2;m++,i++){
					for(j=0,mp=-2;mp<=2;mp++,j++){
						//std::cout<<i<<" "<<j<<" "<<wigner2(phi, theta, gamma, m,mp)<<std::endl;
						powderAs(j,i)=wigner2(phi, theta, gamma, m,mp);
					}
				}
				break;
			case Cartesian:
				cartPowder=rotationMatrix3D(phi, theta, gamma);
				break;
			default:
				for(i=0,m=-2;m<=2;m++,i++){
					for(j=0,mp=-2;mp<=2;mp++,j++){
						//std::cout<<i<<" "<<j<<" "<<wigner2(phi, theta, gamma, m,mp)<<std::endl;
						powderAs(j,i)=wigner2(phi, theta, gamma, m,mp);
					}
				}
				cartPowder=rotationMatrix3D(phi, theta, gamma);
		}
	}
	rotate=false;
}



void Rotations::setSpinnerAs(double ain, double bin, double gin )
{	setRotorAngles(ain, bin, gin );	}


void Rotations::setRotorAngles(double ain, double bin, double gin ){
	if(beta==0 ){
		alpha=ain;
		beta=bin;
		switch(RotationType)
		{
			case Spherical:
				spinnerAs[0]=OneType<complex>::one();
				spinnerAs[1]=OneType<complex>::one();
				spinnerAs[2]=OneType<complex>::one();
				spinnerAs[3]=OneType<complex>::one();
				spinnerAs[4]=OneType<complex>::one();
				break;
			case Cartesian:
				cartSpin.identity();
				break;
			default:
				spinnerAs[0]=OneType<complex>::one();
				spinnerAs[1]=OneType<complex>::one();
				spinnerAs[2]=OneType<complex>::one();
				spinnerAs[3]=OneType<complex>::one();
				spinnerAs[4]=OneType<complex>::one();
				cartSpin.identity();
		}

	}else if(alpha!=ain || beta!=bin){
		alpha=ain; beta=bin;
		switch(RotationType)
		{
			case Spherical:
				spinnerAs[0]=wigner2(alpha,beta,gin, -2,0);
				spinnerAs[1]=wigner2(alpha,beta,gin, -1,0);
				spinnerAs[2]=wigner2(alpha,beta,gin,  0,0);
				spinnerAs[3]=wigner2(alpha,beta,gin,  1,0);
				spinnerAs[4]=wigner2(alpha,beta,gin,  2,0);
				break;
			case Cartesian:
				cartSpin=rotationMatrix3D(alpha, beta);
				break;
			default:
				spinnerAs[0]=wigner2(alpha,beta,gin, -2,0);
				spinnerAs[1]=wigner2(alpha,beta,gin, -1,0);
				spinnerAs[2]=wigner2(alpha,beta,gin,  0,0);
				spinnerAs[3]=wigner2(alpha,beta,gin,  1,0);
				spinnerAs[4]=wigner2(alpha,beta,gin,  2,0);
				cartSpin=rotationMatrix3D(alpha, beta);
		}
		return;
	}
}

void Rotations::setSpinnerAsMat(double ain, double bin, double cin)
{	set_spinnerAsMat(ain, bin, cin);	}

void Rotations::set_spinnerAsMat(double ain, double bin, double cin)
{
	if(ain != alpha || bin !=beta || cin != chi)
	{
		alpha=ain;
		beta=bin;
		chi=cin;
		int m=-2, mp=-2,i=0,j=0;
		for(i=0,m=-2;m<=2;m++,i++){
			for(j=0,mp=-2;mp<=2;mp++,j++){
				//std::cout<<i<<" "<<j<<" "<<wigner2(phi, theta, gamma, m,mp)<<std::endl;
				spinnerAsMat(j,i)=wigner2(alpha, beta, chi, m,mp);
			}
		}
	}
}

Rotations &Rotations::operator=(const Rotations &rhs)
{
	if(this==&rhs) return *this;
	alpha=rhs.alpha;
	beta=rhs.beta;
	gamma=rhs.gamma;
	theta=rhs.theta;
	phi=rhs.phi;
	chi=rhs.chi;
	spinnerAs=rhs.spinnerAs;
	spinnerAsMat=rhs.spinnerAsMat;
	powderAs=rhs.powderAs;
	cartPowder=rhs.cartPowder;
	cartSpin=rhs.cartSpin;
	RotationType=rhs.RotationType;
	return *this;
}

END_BL_NAMESPACE



#endif
