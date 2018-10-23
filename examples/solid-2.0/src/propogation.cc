
/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-26-02
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
	Propogation.cc -->does the propogation on a list of
	Pulse data's...it does not care how they are generated...
	(that is the job for the SequenceParser and other classes)

	important things that need to be declared before using it

	A VEctor of Pulse Datas

	A Map of SolidSYstems sets

	A Map of Powder Angles sets

	A ParamSet class that acts as a bridge between
	the ParerGlobalVars and the system
	(things like wr, maxtstep, etc are defined in that class)

*/

#ifndef __Propogation_cc__
#define __Propogation_cc__ 1

#include "propogation.h"



Propogation::Propogation():
	Systems(NULL), Powders(NULL), Params(NULL),curPow(NULL), curSys(NULL),
	powChanged(false)
{
	//curPow=new powder;
	//curSys=new SolidSys;
}

Propogation::~Propogation()
{Systems=NULL; Powders=NULL; curPow=NULL; curSys=NULL; Params=NULL;}


void Propogation::update2Dtime(int pt)
{
	double curt=0;
	for(int i=0;i<pdata.size();++i){
		for(int j=0;j<pdata[i].amp.size();++j){
			pdata[i].t1=curt;
			//std::cout<<double(pt+1)*pdata[i].dtDelay<<std::endl;
			if(pdata[i].incDelay[j])
			{
				pdata[i].delay[j]=double(pt+1)*pdata[i].dtDelay[j];
			}
			//pdata[i].calcTime(curt);
			pdata[i].t2=pdata[i].t1+pdata[i].delay[j];
			curt=pdata[i].t2;
		}
	}
}

void Propogation::setPowder(std::string use)
{
	static std::string oldUse="";
	powChanged=false;
	if(use==oldUse) return;
	oldUse=use;
	curPow= Powders->getPowder(use);
	powChanged=true;
}

void Propogation::setSystem(std::string use)
{
	static std::string oldUse="";
	if(use==oldUse) return;
	oldUse=use;
	curSys= Systems->getSystem(use);
}

bool Propogation::is2D()
{
//	for(int i=0;i<pdata.size();++i){
//		if(pdata[i].incDelay && pdata[i].dtDelay!=0) return true;
//	}
	return false;
}

matrixs Propogation::propogate(matrixs &ro, Vector<nonUniqueMapEle<std::string,int> > &powpt, int pt)
{
	if(Params==NULL){
		BLEXCEPTION(" No ParamSet specified...")
	}

	if(pt!=-1) update2Dtime(pt);
	int i=0;
	double t1,t2, tau,wr=Params->getD("wr");
	static hmatrixs hh;
	static HamiltonianGen myGen;

	matrixs U=curSys->Fe();

	std::string lastPow="";
	int CurPowPt=-1;
	for(i=0;i<pdata.size();++i)
	{
		wr=Params->getD("wr");
		int tmpp=1,j;
		if(lastPow!=pdata[i].usepow) CurPowPt++;

		setPowder(pdata[i].usepow);
		setSystem(pdata[i].usesys);
		tau=(pdata[i].t2-pdata[i].t1);
		tmpp=(int)std::ceil(tau/Params->getD("maxtstep"));
		if(tau<=Params->getD("mintstep") || tmpp<1){
			tmpp=1;
		}
		if(wr==0) tmpp=1;
		if(tau >1e-10){
			U=curSys->Fe();
			curSys->setRotorAngles(0.0, pdata[i].rotorangle);
			curSys->setPowderAngles(curPow->theta(powpt[CurPowPt].value),
								   curPow->phi(powpt[CurPowPt].value),
								   curPow->gamma(powpt[CurPowPt].value));

			tau/=tmpp;
			t1=pdata[i].t1; t2=pdata[i].t1+tau;
			for(j=0;j<tmpp;j++){
				try{
					hh=(curSys->Hamiltonian(t1, t2,wr)+pdata[i].Pulse(*curSys));
					U*=(Mexp(hh, scomplexi*tau*PI2));
				}catch(BL_exception e){
					e.print(std::cerr);
					BLEXCEPTION("propgation failed");
				}
				t1+=tau;
				t2+=tau;
			}
			ro.prop(U);
			if(pdata[i].cycler !="Fe" && pdata[i].cycler !="Ie")
			{
				matrixs trac=(myGen.Hamiltonian(*curSys, pdata[i].cycler,
												 curPow->theta(powpt[CurPowPt].value),
												 curPow->phi(powpt[CurPowPt].value),
												 curPow->gamma(powpt[CurPowPt].value)));
				//scomplex norm=trace(trac);
				ro=trace(ro,trac)*trac;
			}

		}
	}
	return U;
}

matrixs Propogation::propogate(matrixs &ro, int powpt,  int pt)
{
	if(Params==NULL){
		BLEXCEPTION(" No ParamSet specified...")
	}

	if(pt!=-1) update2Dtime(pt);
	int i=0;
	double t1,t2, tau, wr;
	static hmatrixs hh;
	static HamiltonianGen myGen;

	matrixs U=curSys->Fe();

	for(i=0;i<pdata.size();++i)
	{
		wr=Params->getD("wr");
		int tmpp=1;
		setPowder(pdata[i].usepow);
		setSystem(pdata[i].usesys);
		tau=(pdata[i].t2-pdata[i].t1);

		tmpp=(int)std::ceil(tau/Params->getD("maxtstep"));

		if(tau<=Params->getD("mintstep") || tmpp<1){
			tmpp=1;
		}
		if(wr==0) tmpp=1;
		curSys->setRotorAngles(0.0, pdata[i].rotorangle);
		curSys->setPowderAngles(curPow->theta(powpt), curPow->phi(powpt), curPow->gamma(powpt));
		if(tmpp==1){
			try{
				hh=(curSys->Hamiltonian(pdata[i].t1, pdata[i].t2,wr)+pdata[i].Pulse(*curSys));
				U*=Mexp(hh, scomplexi*tau*PI2);
			}catch(BL_exception e){
				e.print(std::cerr);
				BLEXCEPTION("propgation failed");
			}
		}else{
			U=curSys->Fe();
			tau/=tmpp;
			t1=pdata[i].t1; t2=pdata[i].t1+tau;
			for(int j=0;j<tmpp;j++){
				try{
					hh=curSys->Hamiltonian(t1, t2,wr)+pdata[i].Pulse(*curSys);
					U*=Mexp(hh, scomplexi*tau*PI2);
				}catch(BL_exception e){
					e.print(std::cerr);
					BLEXCEPTION("propgation failed");
				}
				t1+=tau;
				t2+=tau;
			}
		}
		ro.prop(U);
		if(pdata[i].cycler !="Fe" && pdata[i].cycler !="Ie")
		{
			matrixs trac=(myGen.Hamiltonian(*curSys,
			                 pdata[i].cycler,curPow->theta(powpt),
			                 curPow->phi(powpt),
			                 curPow->gamma(powpt)));
			//scomplex norm=trace(trac);
			ro=trace(ro,trac)*trac;
		}
	}
	return U;
}

matrixs Propogation::propogator(int powpt,  int pt)
{
	if(Params==NULL){
		BLEXCEPTION(" No ParamSet specified...")
	}

	//if(pt!=-1) update2Dtime(pt);
	int i=0;
	double t1,t2, tau, wr;
	static hmatrixs hh;
	static dmatrixs dh;

	static matrixs U;
	U=curSys->Fe();

	for(i=0;i<pdata.size();++i)
	{
		wr=Params->getD("wr");
		t1=pdata[i].t1;
		t2=pdata[i].t2;
		tau=(t2-t1);
		int tmpp=1;
		setPowder(pdata[i].usepow);
		setSystem(pdata[i].usesys);
		tmpp=(int)std::ceil(tau/Params->getD("maxtstep"));
		if(tau<=Params->getD("mintstep") || tmpp<1 || wr==0)
		{		tmpp=1;		}

		curSys->setRotorAngles(0.0, pdata[i].rotorangle);
		curSys->setPowderAngles(curPow->theta(powpt), curPow->phi(powpt), curPow->gamma(powpt));
		if(tmpp==1){
			try{
				hh=curSys->Hamiltonian(pdata[i].t1, pdata[i].t2,wr)+pdata[i].Pulse(*curSys);
				U*=Mexp(hh, scomplexi*tau*PI2);
			}catch(BL_exception e){
				e.print(std::cerr);
				std::cerr<<"Current Pulse: "<<std::endl;
				pdata[i].print(std::cerr, 1);
				std::cerr<<std::endl<<"Pulse Hamiltonian::"<<pdata[i].Pulse(*curSys)<<std::endl;
				std::cerr<<std::endl<<"Current Hamiltonian:: "<<curSys->Hamiltonian(t1, t2,wr)<<std::endl;
				BLEXCEPTION("propgation failed");
			}
		}else{
			tau/=tmpp;
			t2=t1+tau;
			for(int j=0;j<tmpp;j++){
				try{
					hh=curSys->Hamiltonian(t1, t2,wr)+pdata[i].Pulse(*curSys);
					U*=Mexp(hh, scomplexi*tau*PI2);
				}catch(BL_exception e){
					e.print(std::cerr);
					std::cerr<<"Current Pulse: "<<std::endl;
					pdata[i].print(std::cerr, 1);
					std::cerr<<std::endl<<"Pulse Hamiltonian::"<<pdata[i].Pulse(*curSys)<<std::endl;
					std::cerr<<std::endl<<"Current Hamiltonian:: "<<curSys->Hamiltonian(t1, t2,wr)<<std::endl;
					BLEXCEPTION("propgation failed");
				}
				//U=hh*U;
				t1+=tau;
				t2+=tau;
			}
		}
	}
	return U;
}

matrixs Propogation::propogator(Vector<nonUniqueMapEle<std::string,int> > &powpt, int pt)
{
	if(Params==NULL){
		BLEXCEPTION(" No ParamSet specified...")
	}

	//if(pt!=-1) update2Dtime(pt);
	int i=0;
	double t1,t2, tau;
	static hmatrixs hh;

	matrixs U=curSys->Fe();

	std::string lastPow="";
	int CurPowPt=-1;
	for(i=0;i<pdata.size();++i)
	{
		int tmpp=1;
		if(lastPow!=pdata[i].usepow) CurPowPt++;

		setPowder(pdata[i].usepow);
		setSystem(pdata[i].usesys);
		while((pdata[i].t2-pdata[i].t1)/tmpp>=Params->getD("maxtstep")){
			tmpp++;
		}
		if((pdata[i].t2-pdata[i].t1)<=Params->getD("mintstep")){
			tmpp=1;
		}
		if(Params->getD("wr")==0) tmpp=1;
		if((pdata[i].t2-pdata[i].t1) >1e-10){
			curSys->setRotorAngles(0.0, pdata[i].rotorangle);
			curSys->setPowderAngles(curPow->theta(powpt[CurPowPt].value),
								   curPow->phi(powpt[CurPowPt].value),
								   curPow->gamma(powpt[CurPowPt].value));

			tau=(pdata[i].t2-pdata[i].t1)/tmpp;
			t1=pdata[i].t1; t2=pdata[i].t2+tau;
			for(int j=0;j<tmpp;j++){
				try{
					hh=curSys->Hamiltonian(t1, t2,Params->getD("wr"))+pdata[i].Pulse(*curSys);
					U*=Mexp(hh, scomplexi*tau*PI2);
				}catch(BL_exception e){
					e.print(std::cerr);
					BLEXCEPTION("Failed to propogate");
				}
				t1+=tau;
				t2+=tau;
			}
		}
	}
	return U;
}

std::ostream &operator<<(std::ostream &oo,const Propogation &out)
{
	if(out.Systems!=NULL)oo<<*(out.Systems)<<std::endl;
	if(out.Powders!=NULL) oo<<*(out.Powders)<<std::endl;
	if(out.Params!=NULL) oo<<*(out.Params)<<std::endl;
	oo<<out.pdata<<std::endl;
	return oo;
}


#endif

