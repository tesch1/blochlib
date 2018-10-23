

#include "blochlib.h"


timer stopwatch;
void printTime(int nrounds=1){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}


int main(int argc, char **argv)
{

	MPIworld.start(argc, argv);
	typedef XYZrect TheShape;
	typedef XYZshape<TheShape> TheGrid;
	typedef double GradType ;
	typedef double FieldType;
//int nsGp=5;
	int nsp=100;
	coord<> mins(-1,-1,-1), maxs(1,1,1);
	coord<int> dims(5,5,5);

	Grid<UniformGrid> gg(mins, maxs, dims);

	TheShape tester(mins, maxs);
	TheGrid jj(gg,tester);

	if(MPIworld.master()) cout<<endl<<"Shape Size: "<<jj.size()<<endl;
	StencilPrep<TheGrid> sp(jj,tester);

	if(MPIworld.master()){
		ofstream oo("org"), oo2("div");
		printTime();
		stopwatch.reset();
		Vector<FieldType > moo(jj.size(),0), hold(moo.size());
		Vector<GradType > koo(moo.size(),0), koo2(koo), koo3(koo), koo4(koo);
		TheGrid::iterator It(jj);
		while(It){
			moo[It.curpos()]=It.x()+It.y()*It.y()+It.z()*It.z();
			++It;
		}

		hold=moo;
		int x=0, y=1, z=2;
		BoundaryCondition<ConstantBC<GradType > > myBC1(GradType(1),sp, BCparams::Zface | BCparams::PositiveFaceNorm | BCparams::WithEdges);
		BoundaryCondition<ConstantBC<GradType > > myBC2(GradType(0),sp, BCparams::Zface | BCparams::NegativeFaceNorm| BCparams::WithEdges);
		BoundaryCondition<ExtrapolateBC<1,GradType > > myBC3(GradType(0),sp);
		BoundaryCondition<EdgeBlendBC<GradType > > myBC4(GradType(0),sp);
		int iterct=1;
		double diff=1;
		int slice=0;
		for(int i=0;i<iterct;++i){
		//	myBC1.applyBC(sp,moo);
		//	myBC2.applyBC(sp,moo);
		//	string fn="div" + itost(i);
		//	ofstream oot(fn.c_str());
			printTime();
		//	for(int i=0;i<koo.size();++i)
		//	{	oot<<jj(i)<<" "<<moo[i]<<endl;	}
			Gradient(koo, sp, moo, x);
			Derivative_1_2n(koo2, sp, moo,x);
			Derivative_1_2n(koo3, sp, moo,y);
			Laplace2Dn(koo4, sp, moo,z);
			//myBC3.applyBC(sp, koo4);
			//myBC4.applyBC(sp, koo4);
			//moo=diff*koo;
		}

		//Edge<TheGrid> ed(jj, tester);	oo3<<ed;

		printTime();
		for(int i=0;i<koo.size();++i)
		{
			oo<<jj(i)<<" "<<hold[i]<<endl;
			oo2<<jj(i)<<" "<<koo4[i]<<endl;
		}
		BoundaryGeneral bg(sp, BCparams::AllFaces | BCparams::AllFaceNorms| BCparams::EdgesOnly);
		BoundaryGeneral::iterator it(bg);
		ofstream oo3("edge");
		while(it){
			oo3<<sp[it()].Point()<<endl;
			++it;
		}
		BoundaryGeneral Inter(sp, BCparams::WithInterior | BCparams::NoEdges); //grab the interior of the grid
		BoundaryGeneral::iterator itI(Inter);
		ofstream ooi("inside");
		while(itI){
			ooi<<sp[itI()].Point()<<" "<<koo[itI()]<<endl;
			++itI;
		}
	}
	MPIworld.end();
}

