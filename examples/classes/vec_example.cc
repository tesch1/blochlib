/*****************************************************************************

 */

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

int main()
{
    const int len = 20;

    Vector<int> a(len), c(len), d(len);
	Vector<double> b(len), y(len);
	Vector<complex> q(len);
    // Set up the vectors with some initial values
    int i;

    for (i=0; i < len; ++i)
    {
        a[i] = i;
        b[i] = i+1.5;
        c[i] = i+10;
        d[i] = i+1;
        q[i] = complex(6+i, 4);
    }

    // Evaluate a vector expression
    //y = (a+b)/(c-d)*9;
	 b=ceil(c+5-b-8+a+d*5/c*1.);
	cout<<"moo"<<endl;
	Vector<double> gg(b);
	Vector<complex> moogoo=FFT(gg);
	cout<<moogoo<<endl<<gg<<endl<<b<<endl;
	//cout<<"moo2"<<gg<<endl<<b+5<<endl<<FFT(gg)<<endl;
	Vector<Vector<complex> > kk(3);
	//kk[1]=b;
	kk[1]=b;kk[0]=gg; kk[2]=gg; kk[2]/=4;
	q=kk[1]+8*b;
	cout<<b<<endl<<kk<<endl<<q<<endl;
	//gg=b*8;
	//gg=b*1;//y=b+5;
	//b=a/8.*d/5-6+b;
	//cout <<" "<<b << endl<<" dot: "<< dot(gg/8, b*9)<<endl<<" Masx: "<<min(gg)<<endl<<"Norm: "<<norm(gg)<<endl;
	//cout<<(gg+5)/norm(gg+5)<<endl;
    return 1;
}

