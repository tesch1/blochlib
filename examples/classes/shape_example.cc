

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

int main(int argv, char **argc){

	int which=1, q=1;
	query_parameter(argv, argc,q++, "Chose the grid [1,2,3]...", which);

	coord<> mins(-1,-1,-1), maxs(1,2,1);
	coord<int,3> dims(30,30,30);
	Grid<UniformGrid> mm(mins,maxs,dims);

	ofstream data("data");
	if(which==1){
		XYZshape<> mastg(mm, XYZfull());
		data<<mastg<<endl;
	}else if(which==2){
		XYZcylinder cyl(0.0,1.0, 0, 3.0*PI/2.0 , -2,2);
		XYZshape<XYZcylinder> mastg(mm, cyl);
		data<<mastg<<endl;
	}else if(which==3){
		XYZrect rect(coord<>(-1,-3,-1), coord<>(0,0,0)); //rectanlg shape
		XYZplaneXY plane(5, 1); //an XY plane
		plane.SetAbove();
		XYZcylinder cyl(0,.5, 0, PI/2, -2, 2); //a cylinderXYZshape<XYZcylinder> mastg(mm, cyl);
		XYZshape<> mastg(mm); //NOTE:: no 'XYZshape' placed in here...

		//THIS IS HOW TO CALCULATE THE GRIDS....
		//the means..."thing in both the rect and plane OR in the cyl"
		mastg.calculate( plane && rect || cyl);
		data<<mastg<<endl;
	}
}

