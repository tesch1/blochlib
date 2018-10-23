
/* utils.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-8-01
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
 	utils.cc-->a collection of string and other utilities
 */



#ifndef _utils_cc_
#define _utils_cc_ 1

#include "utils/utils.h"
#include "utils/constants.h"
#include "container/complex.h"
#include "utils/blassert.h"
#include <string>
#include <time.h>

BEGIN_BL_NAMESPACE


//using namespace std;

void query_parameter(int argc,char* argv[],int par,const std::string& Q, char *V)
{
  if(argc>par) V = argv[par];		// If parameter input, set as std::string
  else { std::cout << Q; std::cin >> V; }		// else ask for it & read it
}


void query_parameter(int argc,char* argv[],int par,const std::string& Q,std::string& V)
{
	if(argc>par) V = argv[par];		// If parameter input, set as std::string
	else { std::cout << Q; std::cin >> V; }		// else ask for it & read it
}


void query_parameter(int argc,char* argv[],int par,const std::string& Q,double& V)
{
	if(argc>par) V = std::atof(argv[par]);	// If parameter input, set as float
	else { std::cout << Q; std::cin >> V; }		// else ask for it & read it
}


void query_parameter(int argc,char* argv[],int par,const std::string& Q,int& V)
{
	if(argc>par) V = std::atoi(argv[par]);	// If parameter input, then set as int
	else { std::cout << Q; std::cin >> V; }		// else ask for it & read it
}

Vector<std::string> parse_param(char * in){
	std::string moo=std::string(in);
	return parse_param(moo);
}

Vector<std::string> parse_param(char * in, char what){
	std::string moo=std::string(in);
	return parse_param(moo, what);
}

Vector<std::string> parse_param(std::string in){
	Vector<std::string> out;
	unsigned int len=in.length();
	//clear all white, tabs and returns in the std::string
	if(in[0]=='\0' || len<=0) {return out;	}
	std::string nn="";
	int gotsumum=0;

	unsigned int i=0,got=0,count=0;;
	//trim any begining white spaces


	for(i=0;i<len;i++){
		if(!isspace(in[i])){
			nn+=in[i];
			count++; got=0;
		}else if(count!=0){
			if(nn[count-1]!='!'){
				nn+='!'; got=1; gotsumum=1;
				count++;
			}
		}
	}
	if(gotsumum==0){out.push_back(nn); return out; }

	std::string hol="", nul="";
	nn+='!';
	unsigned int nc=0;
	for(i=0;i<=count;i++){
		if(nn[i]=='!'){
			out.push_back(std::string(hol));
			nc=0;
			hol=nul;
		}else if(nn[i]==' ' || nn[i]=='\t'){}
		else if(nn[i]!='!'){
			hol+=nn[i];
			nc++;
		}
	}
	if(out.size()>0){
		if(out[out.size()-1].empty()) out.resizeAndPreserve(out.size()-1);
	}
	return out;
}

Vector<std::string> parse_param(std::string in, char what){

	if(what==' ' || what=='\0'){	return parse_param(in);		}
	Vector<std::string> out;
	unsigned int len=in.length();
	//clear all white, tabs and returns in the std::string
	if(in[0]=='\0' || len<=0) {return out;	}
	std::string nn="";
	int gotsumum=0;

	unsigned int i=0,got=0,count=0;;
	for(i=0;i<len;i++){
		if(in[i]!=what){
			nn+=in[i];
			count++; got=0;
		}else if(count!=0){
			if(nn[count-1]!='!'){
				nn+='!'; got=1; gotsumum=1;
				count++;
			}
		}
	}

	if(gotsumum==0){out.push_back(nn); return out; }
	std::string hol="", nul="";
	nn+='!';
	unsigned int nc=0;
	for(i=0;i<=count;i++){
		if(nn[i]=='!'){
			out.push_back(std::string(hol));
			nc=0;
			hol=nul;
		}else if(nn[i]==' '){}
		else if(nn[i]!='!'){
			hol+=nn[i];
			nc++;
		}
	}
	return out;
}

std::string collapsVS(const Vector<std::string> &in){
	return collapsVS(in, 1, in.size());
}

std::string collapsVS(const Vector<std::string> &in, int lim2){
	return collapsVS(in, 1, lim2);
}

std::string collapsVS(const Vector<std::string> &in, int lim1, int lim2){
	if(lim1>int(in.size())) 	lim1=in.size();
	if(lim2>int(in.size())) 	lim2=in.size();
	std::string tt="";
	for(int i=lim1-1; i<lim2;i++)		tt+=in[i];

	return tt;
}

//the whitespace remover
std::string removeWhite(std::string in)
{
	std::string out;
	for(unsigned int i=0;i<in.size();++i){
		if(!isspace(in[i]) ) out+=in[i];
	}
	return out;
}

//gets the stuff inside a set of '(' and ')'
std::string getInside(std::string in)
{
	unsigned int i=0, len=in.size(), levct=0, rict=0,tmct=0;
	std::string tmp="";
	for(i=0;i<len;i++){
		if(in[i]=='('){ levct++;} //count all '('
	}
	if(levct==0) return in;		//no '('
	for(i=0;i<len;i++){
		if(in[i]==')') rict++;
		if(rict!=levct && tmct>0) tmp+=in[i];	//until we hit the same number of ')' as '(', we keep collecting
		else if(rict==levct) break;
		if(in[i]=='(' && tmct==0) tmct++;	//hit the first'('
	}

	if(rict!=levct)
	{
		BLEXCEPTION(std::string(" Parse Error at \"")+in+"\" a '(' ')' problem")
	}
	return tmp;
}

//parses /outermost set of commas
//i.e. 2, moo(4,5), 3
//will be parsed into::
// 2
// moo(4,5)
// 3
Vector<std::string> parseComma(const std::string &in)
{
	Vector<std::string> tmp;
	std::string tt="";
	if(in.empty()) return tmp;
	unsigned int i=0, len=in.size(), pact=0;
	for(i=0;i<len;i++){
		if(in[i]==',' && pact==0){
			tmp.push_back(tt);
			tt="";
		}else{
			tt+=in[i];
		}
		if(in[i]==')')	pact--;
		if(in[i]=='(')	pact++;
	}
	if(!tt.empty()) tmp.push_back(tt);
	return tmp;
}

//this strips out comments (lines starting with '#" and any substr starting with '#)
// and strips out any other text inbetween '{' and '}' sooo the input of
//
//  moo=34
//  loo=3 #moknkey
// #commet
// sec{
//	jj=0
// }

//will become

//  moo=34
//  loo=3
//
Vector<std::string> paramStrip(const Vector<std::string> &pset)
{
	Vector<std::string> out;
	std::string loopArg;
	for(int i=0;i<pset.size();++i)
	{
		loopArg=removeWhite(pset[i]);
		if(loopArg.size()>0 && loopArg[0]!='#' && loopArg[0]!='\n')
		{
			if(loopArg.find("{")<loopArg.size()){
				while( (i-1)<pset.size()){
					loopArg=removeWhite(pset[i]);
					if(loopArg.find("}")<loopArg.size()) break;
					++i;
				}
				if(i==pset.size())
				{
					BLEXCEPTION(" Bad bunch of '}' '{' either one too many or too few...")
				}

			}else{
				out.push_back(pset[i].substr(0,pset[i].find("#")));
			}
		}
	}
	return out;
}


void sttoch(std::string in, char tobe[]){

	int len=in.length();
	char *jj=new char[len];
	for(int i=0;i<len;i++){
		jj[i]=in[i];
	}
	tobe=jj;
}

std::string itost(int i){
	char buffer[80];
	sprintf(buffer, "%d", i);
  	return std::string(buffer);
}

std::string itost(int i, int len)
{
	std::string tm=itost(i);
	std::string loo="";
	if(int(tm.length())<len){
		for(unsigned int i=0;i<len-tm.length();i++){
			loo+=" ";
		}
		loo+=tm;
	}
	return loo;
}


std::string itost_form(const std::string& fmt, int i){
  char buffer[80];
  sprintf(buffer, fmt.c_str(), i);
  return std::string(buffer);
}

std::string dbtost_form(const std::string& fmt, double d){
  char buffer[80];
  sprintf(buffer, fmt.c_str(), d);
  return std::string(buffer);
}

std::string dbtost(double d,const std::string& fmt)
{	return dbtost_form(fmt, d);		}
//ONE D FFT
void FFT1D_(complex data[], int size, int isign){				//FFT OUTPUT WILL BE COMPEX!!!
	if(size==0) return ;
	complex w, wp, temp;
	double theta, sto2;
	int i,j,m,mmax,istep,ii;

	if(isign==-1){
		for(i=0,j=size/2; j<size; i++,j++)	swap_(data[i], data[j]);
	}
	for(j=0,i=0; i<size; i++){
		if(j>i )	swap_(data[i],data[j]);
		m = size/2;
		while((m>1) && (j>=m)){
			j -=m;
			m /= 2;
		}
		j+=m;
	}
	mmax=1;
	while (size>mmax){
		istep = 2*mmax;
		theta = isign*PI2/(mmax);
		sto2 = sin(theta*.5);
		wp = complex(-2.*sto2*sto2, sin(theta));
		w=1;
		for(ii=0; ii<mmax; ii++){
			for(i=ii; i<size; i+=istep){
				j = i+mmax;			// Danielson Lanczos
				if(j<size){
					temp = w*data[j];
					data[j] = data[i] - temp;
					data[i] += temp;
				}
			}
			w+=w*wp;
		}
		mmax=istep;
	}
	if(isign==1){
		for(i=0,j=size/2; j<size; i++,j++)	swap_(data[i], data[j]);
	}

}

//ONE IFFT (inverse FFT)
void IFFT1D_(complex data[], int size){
	FFT1D_(data, size, -1);
}



/*performs an FFT or IFFT on an N dimensional object (i.e. pass the referenace in data)
/ THe INPUT data is assumed real (although stored in a complex array as the output is
/ complex...
/	isign==1-->FFT
/	isign==-1-->IFFT

data souble be ordered like so...
	suppose a NxN matrix...the input should be a VECTOR ordered like

	i  j  data(n)
	-  -  -
	0  0  0
	0  1  1
	0  2  2
	  ...
	1  0  n
	1  1  n+1
	  ...
	i  j  i*n+j
	n  n   2n-1

*/


void FFTN_(complex data[], int nn[], int ndim,int isign)
{
	int idim;
	long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2, adv;
	long ibit,k1,k2,n,nprev,nrem,ntot;
	double theta,wtemp;
	complex w,wp, temp;

	for (ntot=1,idim=0;idim<ndim;idim++)	ntot *= nn[idim];
	nprev=1;
	for (idim=ndim-1;idim>=0;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		//cout<<"n: "<<n<<" nrem: "<<nrem<<" ip1: "<<ip1<<" ip2: "<<ip2<<" ip3: "<<ip3<<endl;
		for (i2=1;i2<=ip2;i2+=ip1) {
			//cout<<" i2: "<<i2<<endl;
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					//cout<<"i1: "<<i1<<endl;
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						//cout<<"i3: "<<i3<<" i3rev: "<<i3rev<<endl;

					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
				//cout<<"i2rev: "<<i2rev<<" ibit: "<<ibit<<endl;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		//cout<<"ifp1: "<<ifp1<<endl;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			//cout<<"ifp2: "<<ifp2<<endl;
			//theta=isign*6.28318530717959/(ifp2/ip1);
			//wtemp=sin(0.5*theta);
			//wpr = -2.0*wtemp*wtemp;
			///wpi=sin(theta);
			//wr=1.0;
			//wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				//cout<<"i3: "<<i3<<endl;
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					//cout<<"i1: "<<i1<<endl;
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						//cout<<"i2: "<<i2<<" k1: "<<k1<<" k2: "<<k2<<endl;
						//tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						//tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						//data[k2]=data[k1]-tempr;
						//data[k2+1]=data[k1+1]-tempi;
						//data[k1] += tempr;
						//data[k1+1] += tempi;
					}
				}
				//wr=(wtemp=wr)*wpr-wi*wpi+wr;
				//wi=wi*wpr+wtemp*wpi+wi;

			}
			ifp1=ifp2;
		}

		nprev *= n;

	}

	//cout<<"MININININI"<<endl;
	for (ntot=1,idim=0;idim<ndim;idim++)	ntot *= nn[idim];
	nprev=1;
	for (idim=ndim-1;idim>=0;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev;//<< 1;
		ip2=(nprev)*n;
		//if(ip1==0) ip2=n;
		ip3=(ip2)*nrem;
		i2rev=0;
		//cout<<"n: "<<n<<" nrem: "<<nrem<<" ip1: "<<ip1<<" ip2: "<<ip2<<" ip3: "<<ip3<<endl;
		adv=(ip1==0)?1:ip1;
		for (i2=0;i2<ip2;i2+=ip1) {
			//cout<<" i2: "<<i2<<endl;
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-1;i1++) {
					//cout<<"i1: "<<i1<<endl;
					for (i3=i1;i3<ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						swap_(data[i3],data[i3rev]);
						//cout<<"i3rev: "<<i3rev<<" i3: "<<i3<<endl;

					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit > ip1 && i2rev >= ibit) {
				i2rev -= ibit;
				ibit >>= 1;
				//cout<<"i2rev: "<<i2rev<<" ibit: "<<ibit<<endl;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
//cout<<"LLLLLL ifp1: "<<ifp1<<endl;

		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			//cout<<"ifp2: "<<ifp2<<endl;
			theta=isign*PI/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wp = complex(-2.0*wtemp*wtemp, sin(theta));
			w=1.0;
			for (i3=0;i3<ifp1;i3+=ip1) {
				//cout<<"i3: "<<i3<<endl;
					for (i1=i3;i1<=i3+ip1-1;i1++) {
					//cout<<"i1: "<<i1<<endl;
					for (i2=i1;i2<ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						//cout<<"i2: "<<i2<<" k1: "<<k1<<" k2: "<<k2<<endl;
						temp = w*data[k2];
						data[k2] = data[k1] - temp;
						data[k1] += temp;
					}
				}
				w+=w*wp;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}



bool BinaryWriteString(std::string out, std::fstream &oo){
	if(out.length()<=0){
		std::cerr<<std::endl<<"Error: BinaryWriteString(string)...bad string"<<std::endl;
	}
	for(unsigned int i=0;i<out.length();i++){
		oo.write(&out[i], sizeof(char));
	}
	return true;
}

std::string BinaryReadString(std::fstream &oo, int len){
	if(len<=0){
		std::cerr<<std::endl<<"Error: BinaryReadString(int)...bad legnth"<<std::endl;
	}
	char *tm; tm=new char[len];
	for(int i=0;i<len;i++){
		oo.read(&tm[i], sizeof(char));
	}
	return std::string(tm);
}


timer::timer() { resolution=HighRes;  reset();}
timer::timer(timer_type myRes) {resolution=myRes; reset();  }

void timer::reset() {
	if(resolution==HighRes){
		time_store=gettime();
	}else{
		low_time=std::time(NULL);
	}
}

double timer::operator()() const {
	if(resolution==HighRes){
		return (gettime()-time_store)/double(CLOCKS_PER_SEC);
	}else{
		std::time_t t2=std::time(NULL);
		return std::difftime(t2, low_time);
	}
}

clock_t timer::gettime() const
{
	return clock();
}


END_BL_NAMESPACE



#endif
