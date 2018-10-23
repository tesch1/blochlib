#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

int main(){
	rmatrix fm(4,4, 2);
	rdmatrix dm(4,4,67);
	smatrix sm(4,4,10);
	fm=cos(dm)+sm;
	rmatrix::iterator j(fm.rows(), fm.cols());
	
	cout<<9*dm/fm*5;


//a special 'matrix*vector<coord<> >' multiply function
	Vector<coord<> > moo(3);
	Random<UniformRandom<coord<> > > myRC(-3,4);
	moo.apply(myRC);

	Random<UniformRandom<> > myR(-3,3);
	Vector<double > joo(9);
	int ct=0, i, k;
	for(k=0;k<moo.size();++k)
		for(i=0;i<3;++i)
			joo[ct++]=moo[k][i];

	rmatrix koo(4, 4), hhh(4,4), out(4,4);
	koo.apply(myR);
	hhh.apply(myR);
	cout<<"*KOO ..."<<koo<<endl;
	cout<<"*HHH ..."<<hhh<<endl;
	out = koo*hhh;
	cout<<"*OUT = KOO * HHH ..."<<out<<endl;
	for(int i=0;i<10;++i) out *= hhh;
	cout<<"#* out *= HHH 10 times..."<<out<<endl;
	for(int i=0;i<10;++i) out = hhh * out;
	cout<<"#* out = HHH *out 10 times..."<<out<<endl;
	
/*
	cout<<"vector<coord<> > before"<<endl;
	moo.print(cout, "\n");
	Vector<coord<> > loo(3, 56); loo=koo*moo;
	cout<<"matrix *vector<coord<> >"<<endl;
	loo.print(cout, "\n");

	cout<<"vector<double > before"<<endl;
	joo.print(cout, "\n");
	Vector<double > mooo(9, 56); mooo=koo*joo;
	cout<<"matrix *vector<double >"<<endl;
	mooo.print(cout, "\n");
*/
//LU decomposition
	rmatrix L,U;
	koo.resize(4,4);
	cout<<"*koo Before LU..."<<koo<<endl;
	koo.LU(L,U);//NON destructive
	cout<<"*L and U"<<L<<U;
	koo.LU(); //destrcutive (L&U are in koo)
	cout<<"*destructive LU"<<koo;

	Vector<Vector<double> > Vtest(4, Vector<double>(Spread<double>(1,5, 1)));
	Vtest.print(cout, "\n"); cout<<endl;
	Vector<double> Vsum; Vsum=sum(Vtest,2);
	Vsum.print(cout, "\n");
	Vsum=sum(Vtest,1);
	Vsum.print(cout, "\n");
	Vsum(Range(0,0)).print(cout, "\n"); cout<<endl;

}
