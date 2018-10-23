

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

timer stopwatch;
void printTime(int nrounds=1){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}


/* A little program that calculates diffusion in 1 Dimensions
via the 'driect' method (requires many times steps)
or
the 'Crank-Nichloson' method (requires fewer time steps but more memory)
*/

int main(int argc, char **argv)
{

	int rr=1;
	int len=500;

	double min=0, max=1, step=(max-min)/len;
	Vector<double> x0(len,0), sol(len,0), axis(Spread<double>(min,max,step));

	double dx=axis[1]-axis[0];
	double dxsq=dx*dx;


	//diffusion constant
	double diff=500;
	query_parameter(argc,argv,rr++, "Enter Diffusion constant: ", diff);
	//initial condition
	for(int i=0;i<len;i++)
	{
		//if(axis[i]>0.25 && axis[i]<0.75) x0[i]=axis[i]+2.*max;
		x0[i]=pow(sin(20*axis[i]*axis[i]), 2.0);
	}

	//boundaries
	//double BCm1=x0[0], BC1=x0[len-1];
	double BCm1=x0[len-1], BC1=x0[0];


	int dirch=0;
	query_parameter(argc,argv,rr++, "Boundar Condition [1=Dirchlet, 2=Neumann, 3=periodic]: ", dirch);


	int iter=100;


int inttype=1;
	query_parameter(argc,argv,rr++, "Integration method \n\t 1=2nd Ord Explicit, \n\t 2=Crank-Nicholson--::", inttype);

	if(inttype==1) //2nd Order Explicit method..typically quite slow...

	{
		double maxt=1, mint=0,t=mint;
		query_parameter(argc,argv,rr++, "Enter End Time: ", maxt);
		double dt=(maxt-mint)/double(iter);

		//this is the 'dt' criterion of numerical accuracy/stability...
		// 2D dt/dx^2 <1
		while(dt>dxsq/(2.0*abs(diff))){	dt/=2;	}

		iter=int(maxt/dt);
		cout<<endl<<"Simulation will require "<<iter<<" steps"<<endl;

		Vector<double> r(len,0), p(len,0),q(len,0);

		int grabmax=20;
		query_parameter(argc,argv,rr++, "capture 'n' points...n= ", grabmax);
		int maxiter=iter, grabmod=maxiter/grabmax, grabidx=1, ct=1;
		rmatrix data(len, grabmax);
		data.putCol(0,x0);

		double fact=dt*diff/dxsq;
		switch(dirch){
			case 1:	x0[0]=BCm1; x0[len-1]=BC1;	break; //dirichlet
			case 2:	x0[0]=x0[1]; x0[len-1]=x0[len-2]; break; //neumann
			case 3: x0[0]=x0[len-1]; break;	//periodic
		}


		while(ct<maxiter){


			x0=fact*Derivative_2_2n(x0)+x0;
			//x0=sol;
			t+=dt;
			ct++;
			switch(dirch){
				case 1:	x0[0]=BCm1; x0[len-1]=BC1;	break; //dirichlet
				case 2:	x0[0]=x0[1]; x0[len-1]=x0[len-2]; break; //neumann
				case 3: //periodic
					x0[0]=x0[0] + fact*(x0[1]-2.0*x0[0]+x0[len-1]);
					x0[len-1]=x0[len-1] + fact*(x0[0]-2.0*x0[len-1]+x0[len-2]);
					break;
			}

			if(ct%grabmod==0 && ct> grabmod)
			{
				cout<<ct<<" "<<t<<"\r"; cout.flush();
				data.putCol(grabidx,x0); grabidx++;
				if(grabidx > grabmax-1) break;
			}
		}
		matstream oo2("div.mat", ios::out | ios::binary);
		oo2.put("data", data);
		oo2.put("axis", axis);
	}else if(inttype==2){ //Crank-Nicholson scheme....
		int grabmax=20;
		query_parameter(argc,argv,rr++, "capture 'n' time slices...n= ", grabmax);
		int maxiter=grabmax,ct=1;

		rmatrix data(len, grabmax);
		data.putCol(0,x0);

		double maxt=1, mint=0;
		query_parameter(argc,argv,rr++, "Enter End Time: ", maxt);
		double dt=(maxt-mint)/double(maxiter), t=mint;

		double fact=dt*diff/dxsq;

		rtrimatrix A(len,len,0.0), B(len, len, 0.0);
		//cout<<A<<endl;

//set up the TriDiagonal Matrix to hold the
		Range Left(1, len), Right(0, len-1), Diag(0, len);
		A.putD(Diag, 1.0+fact);
		A.putR(Right, -fact/2.0);
		A.putL(Left, -fact/2.0);

		B.putD(Diag, 1.0-fact);
		B.putR(Right, fact/2.0);
		B.putL(Left, fact/2.0);


		//Decomp=A;

		Vector<int> perm(len,0);
		//Decomp.LUdecomp(perm);
		switch(dirch){
			case 1:	x0[0]=BCm1; x0[len-1]=BC1;	break; //dirichlet
			case 2:	x0[0]=x0[1]; x0[len-1]=x0[len-2]; break; //neumann
			case 3:
				x0[0]=x0[0] + B(0, 1)*x0[len-1] ;
				x0[len-1]=x0[len-1] + B(len-1, len-2)*x0[0];
				break;
		}
		while(ct<maxiter){
			sol=B*x0;
			//if(dirch){ sol[0]=BCm1; sol[len-1]=BC1;	}
			//else{ sol[0]=sol[1]; sol[len-1]=sol[len-2];	}
			//x0=sol;
			switch(dirch){
				case 1:	x0[0]=BCm1; x0[len-1]=BC1;	break; //dirichlet
				case 2:	x0[0]=x0[1]; x0[len-1]=x0[len-2]; break; //neumann
				case 3:
					x0[0]=x0[0] + B(1, 0)*x0[len-1] ;
					x0[len-1]=x0[len-1] + B(len-2, len-1)*x0[0];
					break;
			}
			x0=A.solve(sol); //A*x=b

			cout<<ct<<" "<<t<<"\r"; cout.flush();
			data.putCol(ct, x0);
			t+=dt;
			ct++;
		}
		matstream oo2("div.mat", ios::out | ios::binary);
		oo2.put("data", data);
		oo2.put("axis", axis);

	}
	printTime(1);

}
