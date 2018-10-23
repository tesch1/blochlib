

/* boundary_condition.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-15-01
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

 	boundary_neumann.h-> implimentation that sets the first derivative at a boundary to 0
 */


#ifndef _boundary_neuamann_h_
#define _boundary_neumann_h_ 1

#include "stencils/stencilprep.h"
#include "stencils/boundary_condition.h"
#include "stencils/boundary_general.h"


#include "container/grids/edgedetect.h"
#include "container/grids/xyzshape.h"




BEGIN_BL_NAMESPACE



/************* NEUMANN BC *********************

These require the first derivate along direction 'r' to be 0
becuase we do not have a functional form for the first derivative
we approximate it by taking a 'Foward Difference' approach

starting at an edge...the approximate first derivative derivative
at the edge is....
dAo/dr=1/2[ 4 A1 - 3 Ao - A2]==g

so

Ao=2/3[A2-4A1]

This is how we calculate our boundry...of course we need a 'neighbor'
and a 'next-neighbor' to do it...so we check for those first...
*/

template<class Numt>
class NeumannBC;


class NeumannBC<double> :
	public BoundaryGeneral
{
	private:
		double const_;
		const void SizeErr()
		{
			BLEXCEPTION(" Data size and StencilPrep size are NOT the same...")
		}
	public:

		typedef double numtype;

		NeumannBC():
			BoundaryGeneral(),
			const_(1)
		{}

		NeumannBC(const numtype &in):
			BoundaryGeneral(),
			const_(in)
		{}

		template<class Shape_t, int Dim>
		NeumannBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph, inopts),
			const_(in)
		{}

		template<class Shape_t, int Dim, class Shapeexpr>
		NeumannBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph, const ShapeExpr<Shapeexpr> &se,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph,se,inopts),
			const_(in)
		{}


		~NeumannBC(){ 	}

		template<class Shape_t, int Dim>
		void applyBCOnce(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter, typename BoundaryGeneral::iterator &myit)
		{
			double ct=const_*2./3.;
			switch(faces)
			{
				case BCparams::AllFacesInt:
					for(int i=0;i<Dim;i++)
					{
						switch(facenormals)
						{
							case BCparams::AllFaceNormsInt:
								if(test(myit()).edge.IsFace(-1,i)){
									if(test(myit()).nextnear(2,i)!=-1 && test(myit()).near(1,i) !=-1){
										alter[myit()]=ct*(alter[test(myit()).nextnear(2,i)]-4.*alter[test(myit()).near(1,i)])/test.d(i, myit());
									}else if( test(myit()).near(1,i) !=-1){
										alter[myit()]=ct*(-4.*alter[test(myit()).near(1,i)])/test.d(i, myit());
									}
								//if there is an to the right of dim 'i' then the BC is calculated to the LEFT
								}else if(test(myit()).edge.IsFace(1,i)){
									if(test(myit()).nextnear(-2,i)!=-1 && test(myit()).near(-1,i) !=-1){
										alter[myit()]=ct*(alter[test(myit()).nextnear(-2,i)]-4.*alter[test(myit()).near(-1,i)])/test.d(i, myit());
									}else if( test(myit()).near(-1,i) !=-1){
										alter[myit()]=ct*(-4.*alter[test(myit()).near(-1,i)])/test.d(i, myit());
									}
								}
								break;
							case BCparams::PositiveFaceNormInt:
							case BCparams::NegativeFaceNormInt:
							default:
								//if(test(myit()).edge.IsFace(-facenormals,i))
								//{
									if(test(myit()).nextnear(-2*facenormals,i)!=-1 && test(myit()).near(-facenormals,i) !=-1){
										alter[myit()]=ct*(alter[test(myit()).nextnear(-2*facenormals,i)]-4.*alter[test(myit()).near(facenormals,i)])/test.d(i, myit());
									}else if( test(myit()).near(-facenormals,i) !=-1){
										alter[myit()]=ct*(-4.*alter[test(myit()).near(-facenormals,i)])/test.d(i, myit());
									}
								//}
								break;
						}
					}
					break;
				case BCparams::XfaceInt:
				case BCparams::YfaceInt:
				case BCparams::ZfaceInt:
				default:
					switch(facenormals)
					{
						case BCparams::AllFaceNormsInt:
							if(test(myit()).edge.IsFace(-1,faces)){
								if(test(myit()).nextnear(2,faces)!=-1 && test(myit()).near(1,faces) !=-1){
									alter[myit()]=ct*(alter[test(myit()).nextnear(2,faces)]-4.*alter[test(myit()).near(1,faces)])/test.d(faces, myit());
								}else if( test(myit()).near(1,faces) !=-1){
									alter[myit()]=ct*(-4.*alter[test(myit()).near(1,faces)])/test.d(faces, myit());
								}
							//if there is an to the right of dim 'i' then the BC is calculated to the LEFT
							}else if(test(myit()).edge.IsFace(1,faces)){
								if(test(myit()).nextnear(-2,faces)!=-1 && test(myit()).near(-1,faces) !=-1){
									alter[myit()]=ct*(alter[test(myit()).nextnear(-2,faces)]-4.*alter[test(myit()).near(-1,faces)])/test.d(faces, myit());
								}else if( test(myit()).near(-1,faces) !=-1){
									alter[myit()]=ct*(-4.*alter[test(myit()).near(-1,faces)])/test.d(faces, myit());
								}
							}
							break;
						case BCparams::PositiveFaceNormInt:
						case BCparams::NegativeFaceNormInt:
						default:
							//if(test(myit()).edge.IsFace(facenormals,faces))
							//{
								//goes the 'oposite' direction from the 'normal' edge (hence the negative signs)
								if(test(myit()).nextnear(-2*facenormals,faces)!=-1 && test(myit()).near(-facenormals,faces) !=-1){
									alter[myit()]=ct*(alter[test(myit()).nextnear(-2*facenormals,faces)]-4.*alter[test(myit()).near(-facenormals,faces)])/test.d(faces, myit());
								}else if( test(myit()).near(-facenormals,faces) !=-1){
									alter[myit()]=ct*(-4.*alter[test(myit()).near(-facenormals,faces)])/test.d(faces, myit());
								}
							//}
							break;
					}	//if there is an to the right of dim 'i' then the BC is calculated to the LEFT
				break;
			}
			cout<<"Post:" <<alter(myit())<<endl;
		}

		template<class Shape_t, int Dim>
		void applyBC(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{
			if(alter.size() != test.size()) SizeErr();

			//HELP ME
			typename BoundaryGeneral::iterator myit(*this);
			while(myit)
			{
				//Ao=2/3[A2-4A1]
				//need to find out which direction is the edge....
				applyBCOnce(test, alter, myit);
				++myit;
			}
		}


		template<class Shape_t, int Dim>
		inline void apply(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{	applyBC(test, alter);	}


		template<class Shape_t, int Dim>
		void operator()(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{	applyBC(test, alter);	}


};



END_BL_NAMESPACE
#endif


