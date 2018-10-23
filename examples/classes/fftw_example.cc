/*****************************************************************************

Testing the FFTW fourier transform routines

 */

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

const int len=8;
int arrlen[len]={2, 8, 16, 24, 32, 35, 90,1110};


int main()
{

/*** uncomment to Test the Vector FFTs **/
    Vector<complex> ans, ff;
    int i=6;
 //   for(int i=0;i<len;++i){
		Vector<double> q(arrlen[i]);
		Vector<complex> qq;
		cout<<"%Original: length="<<arrlen[i]<<endl<<"org=[";
		for(int j=0;j<q.size();++j){
			q[j]=cos(400*2*Pi*j/q.size());//+cos(343*2*Pi*j/q.size());
			cout<<""<<q(j)<<" ";
		}
		cout<<"]"<<endl;

		ff=fftshift(fft(q));
		cout<<"%FFT: length="<<arrlen[i]<<endl<<"ans1=[";
		for(int j=0;j<ff.size();++j)		cout<<"complex"<<ff(j)<<" ";
		cout<<"]"<<endl;
		qq=ifft(ifftshift(ff));
		cout<<"% IFFT: length="<<arrlen[i]<<endl<<"ans2=[";
		for(int j=0;j<qq.size();++j)	cout<<"complex"<<qq(j)<<" ";
		cout<<"]"<<endl;
	//}


 /*   matrix ans, ff;
    int i=3;
	matrix q(arrlen[i],arrlen[i]);
	for(int j=0;j<q.rows();++j){
		for(int k=0;k<q.cols();++k){
			q(j,k)=cos(400*2*Pi*(j*q.rows()+k)/(q.rows()*q.cols()));//+cos(343*2*Pi*j/q.size());
		}
	}
	cout<<"%Original: length="<<arrlen[i]<<endl<<"org="<<q<<endl;

	ff=fftshift(fft(q));
	cout<<"%FFT: length="<<arrlen[i]<<endl<<"ans1="<<ff<<endl;

	q=ifft(ifftshift(ff*55));
	cout<<"% IFFT: length="<<arrlen[i]<<endl<<"ans2="<<q<<endl;*/


	return 0;
}

