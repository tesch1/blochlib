

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

double a=1.28;
double b=0.3;

void func(Vector<double> &X)
{
	static double tx;
	tx=X[0];
	X[0]=a-X[0]*X[0]+b*X[1];
	X[1]=tx;
}

rmatrix Jc(2,2,0);
void Jacobi(Vector<double> &X)
{
	Jc(0,0)=-2*X[0];
	Jc(0,1)=b;
	Jc(1,0)=1;
	Jc(1,1)=0;
}


int main()
{
	ofstream dat("data");
	ofstream inf("inf");
	ofstream ly("lyp");
	Vector<double> IC(2,0);
	IC[0]=0.78;
	IC[1]=0.94;

	rmatrix Jt(2,2,0);
	Jt.identity();

	Vector<double> lyp(2,0);

	int maxi=100, maxj=100;
	int i=0, j=0,kk=0;;
	double xdiv=2.5, ydiv=2.5;
	double xs=-xdiv, ys=-ydiv;
	double cap=100;
	while(j<maxj)
	{
		i=0;
		while(i<maxi)
		{
			kk=0;
			xs=-xdiv+xdiv*2.0*double(j)/double(maxj);
			IC[0]=xs;
			ys=-ydiv+ydiv*2.0*double(i)/double(maxi);
			IC[1]=ys;
			while(kk<10){
				func(IC);
				if(norm(IC)>cap){ inf<<xs<<" "<<ys<<" "<<kk<<endl; break; }
				kk++;
			}
			if(norm(IC)<cap)
			{
				kk=0;
				Jt.identity();
				lyp=0;
				while(kk<50)
				{
					if(norm(IC)>cap) break;
					func(IC);
					Jacobi(IC);
					Jt=Jt*Jc;
					Jt=GramSchmidt(Jt);
					lyp[0]+=log(norm(Jt.col(0)));
					Jt.putCol(0, Jt.col(0)/norm(Jt.col(0)));
					lyp[1]+=log(norm(Jt.col(1)));
					Jt.putCol(1, Jt.col(1)/norm(Jt.col(1)));
					if(kk==49) dat<<IC<<endl;
					kk++;
				}
				ly<<lyp/(50)<<endl;
			}
			i++;
		}

		j++;
	}
}

