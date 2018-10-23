

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

int main(){
	coord<> mins(-1,-1,-1), maxs(1,1,1);
	coord<int,3> dims(10,10,10);
	Grid<UniformGrid> mm(mins,maxs,dims);
	Grid<UniformGrid>::iterator it(mm);
	ofstream loo("fid");
	while(it){
		loo<<it.Point()<<endl;;
		++it;
	}

	maxs(1,PI, 2*PI);
	mins(0,0,0);
	dims(10,10,10);
	Grid<FullCubeGrid> mms(1,1,2);
	Grid<FullCubeGrid>::iterator its(mms);
	mms.Translate(10,0,0);
	ofstream loos("fids");
	while(its){
		loos<<its.Point()<<endl;;
		++its;
	}
}

