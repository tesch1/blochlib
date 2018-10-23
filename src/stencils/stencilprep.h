/* stencilprep.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08-27-01
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
 	stencilprep.h-->contains the classes nessesary to
 	create the data strucutres to quickly impliment an arbitrary stencil
 	(a stencil being a metric that is applied to the 'array' input to result
 	in a single value return...things like gradients, derivative, etc)

 	becuase the XYZshape functions generate completely GENERAL shapes
 	with arbitrary edges and patterns, the stencil MUST be carfeully applied
 	(if the shape is completely rectangular in nature, the application of the stencil
 	is quite simple)...edges, bounds, interiors, neighbors for arbitrary shapes must
 	be known in advance for a stencil operations to run at maximum efficency.

 	this class sets up the nesseary data structure for the optimimum speed for
 	3D grids (i.e. cartesian space)...later these calsses may be applied to
 	arbitrary dimension grids....

 	the initialization cost for this class is VERY HIGH, due to the large amount
 	of searching it must do for every point in the grid...but hopefully this
 	class needs to be applied only ONCE, and the remaining taxing mathematical operations
 	(like solving large differential equations with stencils) should run MUCH faster


 	--The basic thrust of these classes is to find, for any given point in the grid, the nearest neighbors
 	 the next-nearest-neghbors in all 3 directions, and if it is an edge or not...

 	 it does NOT use any points explicitly, only indexes to the points within a grid data set



 	 //there are 3 data strucutre classes
 	 	1) the edge class (defined in 'edgedetect.h')
 	 	2) the nearest neighbor class
 	 	3) the next nearest neighbor class

 	 each sub class maintains a pointer to each nessesary data element
 	 in the grid (as well as its index)

 	 finallt 'StencilPrepData' merges all three into one large data structure


 */


#ifndef _stencil_prep_h_
#define _stencil_prep_h_ 1

#include "container/grids/coords.h"
#include "container/grids/edgedetect.h"
#include "container/grids/xyzshape.h"
#include "mpi/mpi_controller.h"


BEGIN_BL_NAMESPACE



//The nearest Neighbor data
template<class Data_t, int N=3>
class StencilNearest
{
	public:
		typedef Data_t numtype;

		// +dx (0)
		// -dx (1)
		// +dy (2)
		// -dy (3)
		// +dz (4)
		// -dz (5)
		coord<int,2*N> nearidx;	//nearest neighbors INDEX (2 possible the +dr and -dr)
		coord<Data_t, N> *neardata[2*N]; //the ptrs to the grid point data

		StencilNearest():
			nearidx(-1)
		{
			for(int i=0;i<2*N;++i) neardata[i]=NULL;
		}

		StencilNearest(coord<int,2*N> &in):
			nearidx(in)
		{
			for(int i=0;i<2*N;++i) neardata[i]=NULL;
		}

		StencilNearest(coord<Data_t, N> in[2*N]):
			nearidx(-1)
		{
			for(int i=0;i<2*N;++i) neardata[i]=&in[i];
		}


		StencilNearest(coord<int,2*N> &in, coord<Data_t, N> dat[2*N]):
			nearidx(in)
		{
			for(int i=0;i<2*N;++i) neardata[i]=&dat[i];
		}

		StencilNearest(const StencilNearest &cp):
			nearidx(cp.nearidx)
		{
			for(int i=0;i<2*N;++i) neardata[i]=cp.neardata[i];
		}


		~StencilNearest()
		{
			for(int i=0;i<2*N;++i) neardata[i]=NULL;
//			delete neardata;
		}


		StencilNearest &operator=(const StencilNearest &rhs)
		{
			if(&rhs==this)	return *this;
			for(int i=0;i<2*N;++i) neardata[i]=rhs.neardata[i];
			nearidx=rhs.nearidx;
			return *this;
		}

	//the send MPI commands for easy sends
		int put(MPIcontroller &sender, int to)
		{
			return sender.put(nearidx, to);
		}

	//the get MPI commands for easy gets
		int get(MPIcontroller &sender, int from)
		{
			return sender.get(nearidx, from);
		}

	//the send MPI commands for easy sends
		int send(MPIcontroller &sender, int to)
		{
			return sender.send(nearidx, to);
		}

	//the get MPI commands for easy gets
		int recv(MPIcontroller &sender, int from)
		{
			return sender.recv(nearidx, from);
		}

	//sets the Grid Pointers after an MPI call
	//where we could only pass the indexes
		template<class Shape_t>
		void setPointers(Shape_t &shape)
		{
			neardata[0]=(nearidx(0)!=-1)?(&shape.Point(nearidx(0))):NULL;
			neardata[1]=(nearidx(1)!=-1)?(&shape.Point(nearidx(1))):NULL;
			neardata[2]=(nearidx(2)!=-1)?(&shape.Point(nearidx(2))):NULL;
			neardata[3]=(nearidx(3)!=-1)?(&shape.Point(nearidx(3))):NULL;
			neardata[4]=(nearidx(4)!=-1)?(&shape.Point(nearidx(4))):NULL;
			neardata[5]=(nearidx(5)!=-1)?(&shape.Point(nearidx(5))):NULL;
		}

		bool HasNeighbor()
		{
			for(int i=0;i<2*N;++i){	if(nearidx[i]!=-1) return true;	}
			return false;
		}

		bool HasNeighbor(int i)
		{
			if(nearidx[i]==-1) return false;
			return true;
		}

		inline int neighborIndex(int idx, int dim)
		{
			if(idx==-1) return nearidx(dim*2+1);
			else if(idx==1)	return nearidx(dim*2);
			return -1;
		}

		inline coord<Data_t, N> &neighbor(int idx, int dim)
		{
			if(idx==-1) return neardata(dim*2+1);
			else if(idx==1)	return neardata(dim*2);
			return -1;
		}

		int operator()(int idx, int dim)
		{	return neighborIndex(idx,dim);	}


		template<class Shape_t>
		bool FindNeighbor(Shape_t *shape_, int ptidx)
		{
			if(!shape_){	return false;	}
			coord<numtype, N> pt=shape_->Point(ptidx);
			coord<numtype, N> xdp(pt),ydp(pt),zdp(pt), xdm(pt), ydm(pt), zdm(pt);

			typename Shape_t::iterator myit(shape_);

			while(myit)
			{
				xdp=pt; ydp=pt; zdp=pt;
				xdm=pt; ydm=pt; zdm=pt;
				xdp.x()+=myit.dx();
				ydp.y()+=myit.dy();
				zdp.z()+=myit.dz();
				xdm.x()-=myit.dx();
				ydm.y()-=myit.dy();
				zdm.z()-=myit.dz();
				if(myit.Point()==xdp){	nearidx(0)=myit.curpos();	neardata[0]=&myit.Point(); }
				if(myit.Point()==xdm){	nearidx(1)=myit.curpos();	neardata[1]=&myit.Point(); }
				if(myit.Point()==ydp){	nearidx(2)=myit.curpos();	neardata[2]=&myit.Point(); }
				if(myit.Point()==ydm){	nearidx(3)=myit.curpos();	neardata[3]=&myit.Point(); }
				if(myit.Point()==zdp){	nearidx(4)=myit.curpos();	neardata[4]=&myit.Point(); }
				if(myit.Point()==zdm){	nearidx(5)=myit.curpos();	neardata[5]=&myit.Point(); }
				++myit;
			}
			return true;
		}

	//function for determination of weather a given point is a next-nearest neighbor
	//performs this w/o the main loop....so thatone can use the function in another class
	//w/o the penatly of looping through the list 3-4 times...

	//NOTE:: YOU MUST make sure that the 'testpt' is coming from a VALID XYZshape OBJECT!!
	// for the ptr reference to make ANY sence....

		void FindNeighborSingle(
			coord<numtype, N> &pt, //curent point
			coord<numtype, N> &dr, //the delta r
			int pos,	//position of the TEST POINT!!!!
			coord<numtype, N> &testpt //the test point
		)
		{
			//static dec to avoid 'redeclaration costs'
			static coord<numtype, N> xdp(pt),ydp(pt),zdp(pt), xdm(pt), ydm(pt), zdm(pt);

			xdp=pt; ydp=pt; zdp=pt;
			xdm=pt; ydm=pt; zdm=pt;
			xdp.x()+=dr.x();
			ydp.y()+=dr.y();
			zdp.z()+=dr.z();
			xdm.x()-=dr.x();
			ydm.y()-=dr.y();
			zdm.z()-=dr.z();
			//cout<<(1.e6*(testpt-zdm))<<" "<<dr.z()<<endl;
			if(abs(testpt-xdp)<=1.0e-10){	nearidx(0)=pos;	neardata[0]=&testpt; return; }
			if(abs(testpt-xdm)<=1.0e-10){	nearidx(1)=pos;	neardata[1]=&testpt; return; }
			if(abs(testpt-ydp)<=1.0e-10){	nearidx(2)=pos;	neardata[2]=&testpt; return; }
			if(abs(testpt-ydm)<=1.0e-10){	nearidx(3)=pos;	neardata[3]=&testpt; return; }
			if(abs(testpt-zdp)<=1.0e-10){	nearidx(4)=pos;	neardata[4]=&testpt; return; }
			if(abs(testpt-zdm)<=1.0e-10){	nearidx(5)=pos;	neardata[5]=&testpt; return; }
		}

		// prints a little text based item depending on what it has
		void printInfo(std::ostream &oo)
		{

			std::string buf("--------------------------------------   +dr  --------------------------  -dr   -----");
			oo<<buf<<std::endl;
			int ct=0;
			for(int i=0;i<2*N;++i)
			{
				oo<<"Nearest Neighbors along direction "<<ct<<": ";
				if(neardata[i]!=NULL){
					oo<<"["<<*neardata[i]<<"] ";
				}else{
					oo<<setw(8)<<"      [ N/A ]        ";
				}
				++i;
				if(neardata[i]!=NULL){
					oo<<"["<<*neardata[i]<<"] ";
				}else{
					oo<<setw(8)<<"      [ N/A ]";
				}
				oo<<std::endl;
				++ct;
			}
		}
};

template<class data_t, int N>
std::ostream &operator<<(std::ostream &oo, StencilNearest<data_t, N> &out)
{
	out.printInfo(oo);
	return oo;
}



//The nearest Neighbor data
template<class Data_t, int N=3>
class StencilNextNearest
{
	public:

		typedef Data_t numtype;

		// +2 dx (0)
		// -2 dx (1)
		// +2 dy (2)
		// -2 dy (3)
		// +2 dz (4)
		// -2 dz (5)
		coord<int,2*N> nextnearidx;	//nearest neighbors INDEX (2 possible the +dr and -dr)
		coord<Data_t, N> *nextneardata[2*N]; //the ptrs to the grid point data

		StencilNextNearest():
			nextnearidx(-1)
		{
			for(int i=0;i<2*N;++i) nextneardata[i]=NULL;
		}

		StencilNextNearest(coord<int,2*N> &in):
			nextnearidx(in)
		{
			for(int i=0;i<2*N;++i) nextneardata[i]=NULL;
		}

		StencilNextNearest(coord<Data_t, N> in[2*N]):
			nextnearidx(-1)
		{
			for(int i=0;i<2*N;++i) nextneardata[i]=&in[i];
		}

		StencilNextNearest(coord<int,2*N> &in, coord<Data_t, N> dat[2*N]):
			nextnearidx(in)
		{
			for(int i=0;i<2*N;++i) nextneardata[i]=&dat[i];
		}


		StencilNextNearest(const StencilNextNearest &cp):
			nextnearidx(cp.nextnearidx)
		{
			for(int i=0;i<2*N;++i) nextneardata[i]=cp.nextneardata[i];
		}

		~StencilNextNearest()
		{
			for(int i=0;i<2*N;++i) nextneardata[i]=NULL;
//			delete nextneardata;
		}


		StencilNextNearest &operator=(const StencilNextNearest &rhs)
		{
			if(&rhs==this)	return *this;
			for(int i=0;i<2*N;++i) nextneardata[i]=rhs.nextneardata[i];
			//nextneardata=rhs.nextneardata;
			nextnearidx=rhs.nextnearidx;
			return *this;
		}

	//the send MPI commands for easy sends
		int put(MPIcontroller &sender, int to)
		{
			return sender.put(nextnearidx, to);
		}

	//the get MPI commands for easy gets
		int get(MPIcontroller &sender, int from)
		{
			return sender.get(nextnearidx, from);
		}

	//the send MPI commands for easy sends
		int send(MPIcontroller &sender, int to)
		{
			return sender.send(nextnearidx, to);
		}

	//the get MPI commands for easy gets
		int recv(MPIcontroller &sender, int from)
		{
			return sender.recv(nextnearidx, from);
		}

	//sets the Grid Pointers after an MPI call
	//where we could only pass the indexes
		template<class Shape_t>
		void setPointers(Shape_t &shape)
		{
			nextneardata[0]=(nextnearidx(0)!=-1)?&shape.Point(nextnearidx(0)):NULL;
			nextneardata[1]=(nextnearidx(1)!=-1)?&shape.Point(nextnearidx(1)):NULL;
			nextneardata[2]=(nextnearidx(2)!=-1)?&shape.Point(nextnearidx(2)):NULL;
			nextneardata[3]=(nextnearidx(3)!=-1)?&shape.Point(nextnearidx(3)):NULL;
			nextneardata[4]=(nextnearidx(4)!=-1)?&shape.Point(nextnearidx(4)):NULL;
			nextneardata[5]=(nextnearidx(5)!=-1)?&shape.Point(nextnearidx(5)):NULL;
		}

		bool HasNextNeighbor()
		{
			for(int i=0;i<2*N;++i){	if(nextnearidx[i]!=-1) return true;	}
			return false;
		}

		bool HasNextNeighbor(int i)
		{
			if(nextnearidx[i]==-1) return false;
			return true;
		}

		inline int nextneighborIndex(int idx, int dim)
		{
			if(idx==-2) return nextnearidx(dim*2+1);
			else if(idx==2)	return nextnearidx(dim*2);
			return -1;
		}

		inline coord<Data_t, N> &nextneighbor(int idx, int dim)
		{
			if(idx==-2) return nextneardata(dim*2+1);
			else if(idx==2)	return nextneardata(dim*2);
			return -1;
		}

		int operator()(int idx, int dim)
		{	return nextneighborIndex(idx,dim);	}

		template<class Shape_t>
		bool FindNextNeighbor(Shape_t *shape_, int ptidx)
		{
			if(!shape_){	return false;	}
			coord<typename Shape_t::numtype, N> pt=shape_->Point(ptidx);
			coord<typename Shape_t::numtype, N> xdp(pt),ydp(pt),zdp(pt), xdm(pt), ydm(pt), zdm(pt);

			typename Shape_t::iterator myit(shape_);

			while(myit)
			{
				xdp=pt; ydp=pt; zdp=pt;
				xdm=pt; ydm=pt; zdm=pt;
				xdp.x()+=2.0*myit.dx();
				ydp.y()+=2.0*myit.dy();
				zdp.z()+=2.0*myit.dz();
				xdm.x()-=2.0*myit.dx();
				ydm.y()-=2.0*myit.dy();
				zdm.z()-=2.0*myit.dz();
				if(abs(myit.Point()-xdp)<=1.0e-10){	nextnearidx(0)=pos;	nextneardata[0]=&myit.Point(); return; }
				if(abs(myit.Point()-xdm)<=1.0e-10){	nextnearidx(1)=pos;	nextneardata[1]=&myit.Point(); return; }
				if(abs(myit.Point()-ydp)<=1.0e-10){	nextnearidx(2)=pos;	nextneardata[2]=&myit.Point(); return; }
				if(abs(myit.Point()-ydm)<=1.0e-10){	nextnearidx(3)=pos;	nextneardata[3]=&myit.Point(); return; }
				if(abs(myit.Point()-zdp)<=1.0e-10){	nextnearidx(4)=pos;	nextneardata[4]=&myit.Point(); return; }
				if(abs(myit.Point()-zdm)<=1.0e-10){	nextnearidx(5)=pos;	nextneardata[5]=&myit.Point(); return; }

			}
			return true;
		}

	//function for determination of weather a given point is a next-nearest neighbor
	//performs this w/o the main loop....so thatone can use the function in another class
	//w/o the penatly of looping through the list 3-4 times...

	//NOTE:: YOU MUST make sure that the 'testpt' is coming from a VALID XYZshape OBJECT!!
	// for the ptr reference to make ANY sence....

		void FindNextNeighborSingle(
			coord<numtype, N> &pt, //curent point
			coord<numtype, N> &dr, //the delta r
			int pos,	//position of the TEST POINT!!!!
			coord<numtype, N> &testpt //the test point
		)
		{
			//static dec to avoid 'redeclaration costs'
			static coord<numtype, N> xdp(pt),ydp(pt),zdp(pt), xdm(pt), ydm(pt), zdm(pt);

			xdp=pt; ydp=pt; zdp=pt;
			xdm=pt; ydm=pt; zdm=pt;
			xdp.x()+=2.0*dr.x();
			ydp.y()+=2.0*dr.y();
			zdp.z()+=2.0*dr.z();
			xdm.x()-=2.0*dr.x();
			ydm.y()-=2.0*dr.y();
			zdm.z()-=2.0*dr.z();
			if(abs(testpt-xdp)<=1.0e-10){	nextnearidx(0)=pos;	nextneardata[0]=&testpt; return; }
			if(abs(testpt-xdm)<=1.0e-10){	nextnearidx(1)=pos;	nextneardata[1]=&testpt; return; }
			if(abs(testpt-ydp)<=1.0e-10){	nextnearidx(2)=pos;	nextneardata[2]=&testpt; return; }
			if(abs(testpt-ydm)<=1.0e-10){	nextnearidx(3)=pos;	nextneardata[3]=&testpt; return; }
			if(abs(testpt-zdp)<=1.0e-10){	nextnearidx(4)=pos;	nextneardata[4]=&testpt; return; }
			if(abs(testpt-zdm)<=1.0e-10){	nextnearidx(5)=pos;	nextneardata[5]=&testpt; return; }
		}

		// prints a little text based item depending on what it has
		void printInfo(std::ostream &oo)
		{

			std::string buf("--------------------------------------   +dr  --------------------------  -dr   -----");
			oo<<buf<<std::endl;
			int ct=0;
			for(int i=0;i<2*N;++i)
			{
				oo<<"Next Nearest Neighbors along direction "<<ct<<": ";
				if(nextneardata[i]!=NULL){
					oo<<"["<<*nextneardata[i]<<"] ";
				}else{
					oo<<setw(8)<<"      [ N/A ]        ";
				}
				++i;
				if(nextneardata[i]!=NULL){
					oo<<"["<<*nextneardata[i]<<"] ";
				}else{
					oo<<setw(8)<<"      [ N/A ]";
				}
				oo<<std::endl;
				++ct;
			}
		}
};


template<class Data_t, int N>
std::ostream &operator<<(std::ostream &oo, StencilNextNearest<Data_t, N> &out)
{
	out.printInfo(oo);
	return oo;
}


template<class Data_t, int N=3>
class StencilPrepData
{
	private:

	public:

		//order as follows...
		// +dx (0)
		// -dx (1)
		// +dy (2)
		// -dy (3)
		// +dz (4)
		// -dz (5)
		StencilNearest<Data_t,N> near;
		StencilNextNearest<Data_t,N> nextnear;
		int ptidx;		//the index to the 'current' point in the grid list
		coord<Data_t, N> *ptdata;	//the ptr to the grid point
		EdgeData<N> edge;	//is an edge or not...

		StencilPrepData():
			near(),
			nextnear(),
			ptidx(-1),
			ptdata(NULL),
			edge(false)
		{}

		StencilPrepData(int pt_):
			near(),
			nextnear(),
			ptidx(pt_),
			ptdata(NULL),
			edge(false)
		{}

		~StencilPrepData(){	ptdata=NULL;	}

		StencilPrepData &operator=(const StencilPrepData &rhs)
		{
			if(&rhs==this) return *this;
			near=rhs.near;
			nextnear=rhs.nextnear;
			ptidx=rhs.ptidx;
			ptdata=rhs.ptdata;
			edge=rhs.edge;
			return *this;
		}

//the send MPI commands for easy gets and sends
		int put(MPIcontroller &sender, int to)
		{
			int err=sender.put(ptidx, to);
			err=near.put(sender, to);
			err=nextnear.put(sender, to);
			err=edge.put(sender, to);
			return err;
		}

//the et MPI commands for easy gets and sends
		int get(MPIcontroller &sender, int from)
		{
			int err=sender.get(ptidx, from);
			err=near.get(sender, from);
			err=nextnear.get(sender, from);
			err=edge.get(sender, from);
			return err;
		}


//the send MPI commands for easy gets and sends
		int send(MPIcontroller &sender, int to)
		{
			int err=sender.send(ptidx, to);
			err=near.send(sender, to);
			err=nextnear.send(sender, to);
			err=edge.send(sender, to);
			return err;
		}

/*//the et MPI commands for easy gets and sends
		int recv(MPIcontroller &sender, int from)
		{
			int err=sender.recv(ptidx, from);
			err=near.recv(sender, from);
			err=nextnear.recv(sender, from);
			err=edge.recv(sender, from);
			return err;
		}
*/
		template<class Shape_t>
		void setPointers(Shape_t &shape)
		{
			ptdata=&shape.Point(ptidx);
			near.setPointers(shape);
			nextnear.setPointers(shape);
		}

		inline bool HasNeighbor(){		return  near.HasNeighbor();	}
		inline bool HasNeighbor(int i){	return  near.HasNeighbor(i);	}
		inline bool HasNextNeighbor(){		return  nextnear.HasNextNeighbor();	}
		inline bool HasNextNeighbor(int i){	return  nextnear.HasNextNeighbor(i);	}
		inline bool IsFace(){	return edge.IsFace();	}
		inline bool IsFace(int i){	return edge.IsFace(i);	}

		coord<> &Point()	{	return *ptdata;	}
		coord<> Point() const {	return *ptdata;	}

		template<class Shape_t>
		inline bool FindNextNearest(Shape_t *shape, int idx){	return nextnear.FindNextNearest(shape, idx);	}

		template<class Shape_t>
		inline bool FindNextNearest(Shape_t &shape, int idx){	return nextnear.FindNextNearest(&shape, idx);	}

		template<class Shape_t>
		inline bool FindNearest(Shape_t *shape, int idx){	return near.FindNearest(shape, idx);	}

		template<class Shape_t>
		inline bool FindNearest(Shape_t &shape, int idx){	return near.FindNearest(&shape, idx);	}


		void print(std::ostream &oo)
		{
			if(ptdata==NULL)
			{
				oo<<"No Grid data defined...";
				return;
			}

			oo<<"Point: "<<*ptdata<<std::endl;
			if(IsFace())	oo<<edge<<std::endl;
			if(HasNeighbor()) oo<<" "<<near<<std::endl;
			if(HasNextNeighbor()) oo<<" "<<nextnear<<std::endl;
		}
};

template<class Shape_t, int N>
std::ostream &operator<<(std::ostream &oo, StencilPrepData<Shape_t, N> &out)
{
	out.print(oo);
	return oo;
}

//predec for the iterator
template<class StencilPrep_t>
class StencilPrepIter;

//the vectorized version of the data above...to be used for an entire XYZshape
template<class Shape_t, int N=3>
class StencilPrep
{
	private:
		Vector<StencilPrepData<typename Shape_t::numtype, N> > data_;
		Shape_t *shape_;

		const void ShapeErr()
		{
			std::cerr<<"Error: StencilPrep<Shape_t, N>"<<std::endl;
			std::cerr<<" cannot prepare for stencil...(no shape/grid defined...)"<<std::endl;
		}

		const void AccssErr()
		{
			std::cerr<<"Error: StencilPrep<Shape_t, N>::calculate"<<std::endl;
			std::cerr<<" Access element is too large for the diesred grid..."<<std::endl;
		}

	public:
		typedef StencilPrep<Shape_t, N> StencilPrep_t;
		typedef Shape_t shape;
		typedef typename Shape_t::numtype numtype;
		typedef StencilPrepData<typename Shape_t::numtype, N> Data_t;
		typedef StencilPrepIter<StencilPrep_t> iterator;

		static const int datalength=N;

		friend class iterator;

		MPIcontroller StencilControler;

		StencilPrep():
			data_(0),
			shape_(0),
			StencilControler(MPIworld)
		{}

		StencilPrep(Shape_t &ins,const MPIcontroller &contol=MPIcontroller()):
			data_(ins.size(), StencilPrepData<numtype,N>()),
			shape_(&ins), StencilControler(contol)
		{}

		template<class ShapeExpr>
		StencilPrep(Shape_t &ins, ShapeExpr ine,const MPIcontroller &contol=MPIcontroller()):
			data_(ins.size(), StencilPrepData<numtype,N>()),
			shape_(&ins), StencilControler(contol)
		{
			CalcStencilPrep(ins, ine);
		}

		~StencilPrep(){	shape_=NULL; }

		bool FindNeighbor(int ptidx)
		{
			if(!shape_){	ShapeErr();	return false;	}

			if(data_.size() != shape_->size()) data_.resize(shape_->size(), StencilPrepData<numtype, N>());

			if(!data_(ptidx).near.FindNeighbor(shape_, ptidx)) return false;

			return true;
		}

		bool FindNextNeighbor(int ptidx)
		{
			if(!shape_){	ShapeErr();	return false;	}
			if(data_.size() != shape_->size()) data_.resize(shape_->size(), StencilPrepData<numtype, N>());

			if(!data_(ptidx).nextnear.FindNextNeighbor(shape_, ptidx)) return false;

			return true;
		}

		bool FindAllNeighbors(int ptidx)
		{
			if(!shape_){	ShapeErr();	return false;	}
			if(data_.size() != shape_->size()) data_.resize(shape_->size(), StencilPrepData<numtype, N>());

			if(!FindNeighbor(ptidx)) return false;
			if(!FindNextNeighbor(ptidx)) return false;
			return true;
		}

		template<class ShapeExpr>
		bool calculate(Shape_t &ins, ShapeExpr ine){	return CalcStencilPrep(ins, ine);	}

		template<class ShapeExpr>
		bool CalcStencilPrep(Shape_t &ins, ShapeExpr ine)
		{
			if(!ins){ ShapeErr(); return false;	}
			if(&ins!=shape_) shape_=&ins;
			if(data_.size() != ins.size()) data_.resize(ins.size(), StencilPrepData<numtype, N>());
			//test for Parrellel or Not
			typename Shape_t::iterator myit(ins);
			if(StencilControler.serial()){
				while(myit)
				{
					CalcStencilPrep(ins, ine, myit.curpos());
					++myit;
				}
				return true;
			}else{
				int ct=0, done=-1,onPt=1;
				if(StencilControler.master()){
				//send the initial points to calculate
					for(int i=1;i<StencilControler.size();++i){
						ct=myit.curpos();
						StencilControler.put(ct, i);
						++myit;
						if(!myit) break;
					}

				//keep sending to any proc that is willing to get one
					while(myit){
						int getP=StencilControler.getAny(onPt); //get the 'we calced' point
						data_[onPt].get(StencilControler, getP); //get the chunk of data
						data_[onPt].setPointers(ins);
						ct=myit.curpos();
						StencilControler.put(ct,getP); //send the proc a new point to do
						++myit;
					}
				//get the last returns
					for(int qq=1;qq<StencilControler.size();++qq){
						int getP=StencilControler.getAny(onPt); //get the 'we calced' point
						data_[onPt].get(StencilControler, getP); //get the chunk of data
						data_[onPt].setPointers(ins);
					}

				//put the final kills
					for(int qq=1;qq<StencilControler.size();++qq)
						StencilControler.put(done, qq);

				}else{ //the slaves
					while(1)
					{
						StencilControler.get(onPt,0);
						if(onPt==done) break;
						CalcStencilPrep(ins, ine, onPt);
						StencilControler.put(onPt,0);
						data_[onPt].put(StencilControler, 0);
					}
				}

				//now the master must send the data to all the other procs
				// this can be very time consuming...perhaps a better way?
				int curi,i,j, IntTag=20;
				if(StencilControler.master()){
					for(i=0;i<data_.size();++i){
						for(j=1;j<StencilControler.size()&&i<data_.size();++j){
							StencilControler.send(i,j, IntTag);
							data_[i].send(StencilControler,j);
						}
					}
					for(j=1;j<StencilControler.size();++j)
						StencilControler.send(done,j, IntTag);

				}else{
					while(1)
					{
						StencilControler.get(curi, 0, IntTag);
						if(curi==done) break;
						data_[curi].get(StencilControler, 0);
						data_[curi].setPointers(ins);
					}
				}
				StencilControler.barrier();
				return true;
			}
		}

	private: //To be used ONLY from inseide the class
	//the single point calculator for the MPI commands
		template<class ShapeExpr>
		bool CalcStencilPrep(Shape_t &ins, ShapeExpr ine, int which)
		{
			coord<numtype, N> curpt(0), dr(0);
			typename Shape_t::iterator myit(ins, which);
			typename Shape_t::iterator myit2(ins);

			dr=coord<numtype, 3>(myit.dx(), myit.dy(), myit.dz());
			data_(myit.curpos()).ptidx=myit.curpos();
			data_(myit.curpos()).ptdata=&myit.Point();
			data_(myit.curpos()).edge.CalcEdgesSingle(ins, ine, myit.Point(), dr);
			while(myit2)
			{
				data_(myit.curpos()).near.FindNeighborSingle(myit.Point(), dr, myit2.curpos(), myit2.Point());
				data_(myit.curpos()).nextnear.FindNextNeighborSingle(myit.Point(), dr, myit2.curpos(), myit2.Point());
				++myit2;
			}
			return true;
		}

	public:

		bool empty(){	return data_.size()<=0;	}

		inline int size() const {	return data_.size();	}

		Vector<StencilPrepData<numtype, N> > &data(){	return data_;	}
		Vector<StencilPrepData<numtype, N> > data()const{	return data_;	}

		StencilPrepData<numtype, N> &data(int i){	return data_(i);	}
		StencilPrepData<numtype, N> data(int i) const {	return data_(i);	}

		StencilPrepData<numtype, N> operator()(int i) const {	return data_(i);	}
		StencilPrepData<numtype, N> &operator()(int i)  {	return data_(i);	}
		StencilPrepData<numtype, N> operator[](int i) const {	return data_(i);	}
		StencilPrepData<numtype, N> &operator[](int i)  {	return data_(i);	}

		numtype d(int dir, int which){	return shape_->d(dir, which);	}

		numtype dx(int i){	return shape_->dx(i);	}
		numtype dy(int i){	return shape_->dy(i);	}
		numtype dz(int i){	return shape_->dz(i);	}

		StencilNextNearest<numtype,N> &nextnearest(int i)
		{
			return data(i).nextnear;
		}

		StencilNextNearest<numtype,N> &nextnear(int i)
		{
			return data(i).nextnear;
		}

		StencilNearest<numtype,N> &nearest(int i)
		{
			return data(i).near;
		}

		StencilNearest<numtype,N> &near(int i)
		{
			return data(i).near;
		}

		EdgeData<N> &edge(int i)
		{
			return data(i).edge;
		}



		void print(std::ostream &oo)
		{
			iterator myit(this);
			while(myit)
			{
				oo<<myit.Point()<<std::endl;
				++myit;
			}
		}

		void printNeighbors(std::ostream &oo)
		{
			iterator myit(this);
			while(myit)
			{
				oo<<myit.nearest()<<std::endl;
				++myit;
			}
		}

		void printEdge(std::ostream &oo)
		{
			iterator myit(this);
			while(myit)
			{
				oo<<myit.edge()<<std::endl;
				++myit;
			}
		}

		void printNextNeighbors(std::ostream &oo)
		{
			iterator myit(this);
			while(myit)
			{
				oo<<myit.nextnearest()<<std::endl;
				++myit;
			}
		}
};

template<class Shape_t, int N>
std::ostream &operator<<(std::ostream &oo, StencilPrep<Shape_t, N> &out)
{
	out.print(oo);
	return oo;
}

// The mighty iterator for the stencilPrep
template<class StencilPrep_t>
class StencilPrepIter
{
	private:
		StencilPrep_t *sp_;
		typename StencilPrep_t::Data_t *curpt_;

		bool notended_;
		int pos;

	public:

		typedef typename StencilPrep_t::numtype numtype;
		typedef typename StencilPrep_t::Data_t DataType;
		static const int datalength=StencilPrep_t::datalength;
		typedef StencilNextNearest<numtype,datalength> NextNearestPoint;
		typedef StencilNearest<numtype,datalength> NextPoint;
		typedef EdgeData<datalength> EdgePoint;


		StencilPrepIter():
			sp_(NULL),
			curpt_(NULL),
			notended_(false),
			pos(0)
		{}

		StencilPrepIter(StencilPrep_t &in):
			sp_(&in),
			notended_(false),
			pos(0)
		{
			if(!in.empty()){
				curpt_=&in.data(0);
				notended_=true;
			}else{
				notended_=false;
				curpt_=NULL;
			}
		}

		StencilPrepIter(StencilPrep_t *in):
			sp_(in),
			notended_(false),
			pos(0)
		{
			if(!in->empty()){
				curpt_=&in->data(0);
				notended_=true;
			}else{
				notended_=false;
				curpt_=NULL;
			}
		}

		~StencilPrepIter()
		{
			sp_=NULL;
			curpt_=NULL;
		}
		void operator++()
		{
			if(sp_ && pos<sp_->size()-1){
				pos++;
				curpt_=&sp_->data(pos);
			}else{
				notended_=false;
			}
		}

		void operator++(int){	operator++();	}

		inline operator bool(){	return notended_;	}

		typename StencilPrep_t::Data_t &Point(){	return *curpt_;	}
		typename StencilPrep_t::Data_t &operator()(){	return *curpt_;	}

		typename StencilPrep_t::Data_t Point() const {	return *curpt_;	}
		typename StencilPrep_t::Data_t operator()() const {	return *curpt_;	}

		typename StencilPrep_t::Data_t &Point(int i){	return sp_->data(i);	}
		typename StencilPrep_t::Data_t &operator()(int i){	return sp_->data(i);	}

		typename StencilPrep_t::Data_t Point(int i) const {	return sp_->data(i);	}
		typename StencilPrep_t::Data_t operator()(int i) const {	return sp_->data(i);	}

		void moveTo(int i)
		{
			//int &ref=const_cast<int &>( pos);
			//ref=i;
			//DataType &ref2_=const_cast<DataType &>( *curpt_);
			//ref2_=sp_->data(pos);
			pos=i;
			curpt_=&sp_->data(pos);
		}

		inline int curpos(){	return pos;	}

		StencilNextNearest<numtype,datalength> &nextnearest()
		{
			return curpt_->nextnear;
		}

		StencilNextNearest<numtype,datalength> &nextnear()
		{
			return curpt_->nextnear;
		}

		StencilNearest<numtype,datalength> &nearest()
		{
			return curpt_->near;
		}

		StencilNearest<numtype,datalength> &near()
		{
			return curpt_->near;
		}

		EdgeData<datalength> &edge()
		{
			return curpt_->edge;
		}
//Returns the edge strucutre for the nearest neighbor along direction 'dir'
		EdgeData<datalength> NearEdge(int idx,int dir)
		{
			if(!edge().IsFace(dir)){
				return sp_->data(curpt_->near(idx, dir)).edge;
			}else{
				return EdgeData<datalength>();
			}
		}

//Returns the edge strucuture for the next nearest neighbor along direction 'dir'
		EdgeData<datalength> NextNearEdge(int idx,int dir)
		{
			if(!edge().IsFace(dir)){
				return sp_->data(curpt_->nextnear(idx, dir)).edge;
			}else{
				return EdgeData<datalength>();
			}
		}

//Is the nearest neighbor along direction 'dir' an edge
		bool IsNearEdge(int dir)
		{
			if(!edge().IsFace(dir)){
				return NearEdge(-1,dir).IsFace()
					|| NearEdge(1,dir).IsFace();

			}else{
				return true;
			}
		}

//Is the next nearest neighbor along direction 'dir' an edge
		bool IsNextNearEdge(int dir)
		{
			if(!edge().IsFace(dir)){
				return NextNearEdge(-1,dir).IsFace()
					|| NextNearEdge(1,dir).IsFace();
			}else{
				return true;
			}
		}

//Is the nearest neighbor along Any direction an edge
		bool IsNearEdge()
		{

			for(int i=0;i<datalength;++i){
				if(IsNearEdge(i)){  return true;}
			}
			return false;
		}

//Is the NExt nearest neighbor along Any direction an edge
		bool IsNextNearEdge()
		{
			for(int i=0;i<datalength;++i){
				if(IsNextNearEdge(i)) return true;
			}
			return false;
		}

		numtype dx(){	return sp_->dx(pos);	}
		numtype dy(){	return sp_->dy(pos);	}
		numtype dz(){	return sp_->dz(pos);	}
		coord<numtype> dr(){	return coord<numtype>(dx(), dy(), dz());	}

		numtype dr(int which){	return dr()[which]; }

		coord<> &GridPoint()	{	return curpt_->Point();	}

		StencilNextNearest<numtype,datalength> &nextnearest(int i)
		{
			return sp_->data(i).nextnear;
		}

		StencilNextNearest<numtype,datalength> &nextnear(int i)
		{
			return sp_->data(i).nextnear;
		}

		StencilNearest<numtype,datalength> &nearest(int i)
		{
			return sp_->data(i).near;
		}

		StencilNearest<numtype,datalength> &near(int i)
		{
			return sp_->data(i).near;
		}

		EdgeData<datalength> &edge(int i)
		{
			return sp_->data(i).edge;
		}

		numtype dx(int i){	return sp_->dx(i);	}
		numtype dy(int i){	return sp_->dy(i);	}
		numtype dz(int i){	return sp_->dz(i);	}


		inline int size(){	return sp_->size();	}

		template<class OutType, int Idx>
		OutType &shift(int Idx);
};


/*template<class StencilPrep_t>
StencilPrepIter<StencilPrep_t>::NextPoint
	&StencilPrepIter<StencilPrep_t>::shift(int )
{
	return nearest();
}

template<class StencilPrep_t>
StencilPrepIter<StencilPrep_t>::NextPoint
	&StencilPrepIter<StencilPrep_t>::shift<typename  StencilPrepIter<StencilPrep_t>::NextPoint, -1>(int j=-1)
{
	return nearest();
}
*/



END_BL_NAMESPACE






#endif


