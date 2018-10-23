

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

int main(int argv, char **argc)
{


	ofstream oo("moo"), oo2("loo");


//a grid...
	coord<> mins(-0.01,-0.01,-0.01), maxs(0.01,0.01,0.01);
	coord<int> dims(20,20,20);

	typedef XYZrect TheShape1;
	typedef XYZplaneXZ TheShape2;
	typedef XYZrect TheShape3;
	typedef XYZshape<TheShape1> TheGrid;
	TheShape1 t1(mins, maxs);//(0,0.005, 0, PI2, -0.1, 0.1);
	TheShape2 t2(0,.005);
	TheShape3 t3(-0.01, 0, -0.01, 0, -0.01, 0);


	Grid<UniformGrid> gg(mins, maxs, dims);
	//t2.SetBelow();

	//t1.SetBelow();
	TheGrid jj(gg, t1);
	//jj.calculate((t1 || t3) &&t2);
	Edge<TheGrid> ed(jj, (t1 && t3) || t2);
	oo2<<ed<<endl;
	int nsp=jj.size();

	//Vector<double> vals(nsp, 300);
	//RandomOffsetGradFunc<GaussianRandom<> > mr(.1,2);
	int q=1;
	double slope, inter;
	query_parameter(argv, argc, q++, "slope: ", slope);
	query_parameter(argv, argc, q++, "intercept: ", inter);

	TemperatureGrad<LinearOffsetXGradFunc> mr(slope, inter);
	cout<<mr<<endl;
	int i=0;
	ListBlochParams<TheGrid> myp(nsp, "1H",jj);

//apply the gradient function...
	mr(myp, (t1 && t3) || t2);

	ListBlochParams<TheGrid>::iterator myit(myp);

//print sutff out for me to see
	while(myit)
	{
		oo<<myit.Point()<<" ";
		oo<<myit.temperature()<<endl;
		++i;
		++myit;
	}


}
