


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

 	boundary_extrap.h-> the implimentation for the 'extrapolation' BC

 	this is really assuming 'infinite' dimensions along the face
 */


#ifndef _boundary_extrap_h_
#define _boundary_extrap_h_ 1

#include "stencils/stencilprep.h"
#include "stencils/boundary_general.h"

#include "container/grids/edgedetect.h"
#include "container/grids/xyzshape.h"


BEGIN_BL_NAMESPACE


/*********************

Extrapolation...fills in the BC as if it was 'infinite'
the stencils cannot go to the bounardy, so this 'fills' in the blanks..


*********************/

template<int Order, class Numt>
class ExtrapolateBC;

//first order only copies the neighbor point to the boundary
// a silly one, but nessesary for completeness
template<class Num_t>
class ExtrapolateBC<1, Num_t> :
	public BoundaryGeneral
{
	private:
		Num_t const_;
		const void SizeErr()
		{
			BLEXCEPTION(" Data size and StencilPrep size are NOT the same...")
		}
	public:

		typedef Num_t numtype;

		ExtrapolateBC():
			BoundaryGeneral(),
			const_(0)
		{}

		ExtrapolateBC(const numtype &in):
			BoundaryGeneral(),
			const_(in)
		{}

		template<class Shape_t, int Dim>
		ExtrapolateBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph,
			int opts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph, opts),
			const_(in)
		{}

		template<class Shape_t, int Dim, class Shapeexpr>
		ExtrapolateBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph, const ShapeExpr<Shapeexpr> &se,
			int opts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph,se,opts),
			const_(in)
		{}


		~ExtrapolateBC(){ 	}

		template<class Shape_t, int Dim>
		void applyBCOnce(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter, typename BoundaryGeneral::iterator &myit)
		{
			switch(faces)
			{
				case BCparams::AllFacesInt:
					for(int i=0;i<Dim;i++)
					{
						switch(facenormals)
						{
							case BCparams::AllFaceNormsInt:
								if(test(myit()).edge.IsFace(-1,i)){
									//if(test(myit()).near(1,i) !=-1){
										alter[myit()]=alter[test(myit()).near(1,i)]+const_;
									//}
								}else if(test(myit()).edge.IsFace(1,i)){
									//if(test(myit()).near(-1,i) !=-1){
										alter[myit()]=alter[test(myit()).near(-1,i)]+const_;
									//}
								}
								break;
							case BCparams::PositiveFaceNormInt:
							case BCparams::NegativeFaceNormInt:
							default:
								if(test(myit()).near(facenormals,i) !=-1){
									alter[myit()]=alter[test(myit()).near(facenormals,i)]+const_;
								}
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
							if(test(myit()).near(1,faces) !=-1){
								alter[myit()]=alter[test(myit()).near(1,faces)]+const_;
							}else if(test(myit()).near(-1,faces) !=-1){
								alter[myit()]=alter[test(myit()).near(-1,faces)]+const_;
							}

							break;
						case BCparams::PositiveFaceNormInt:
						case BCparams::NegativeFaceNormInt:
						default:
							if(test(myit()).near(facenormals,faces) !=-1){
								alter[myit()]=alter[test(myit()).near(facenormals,faces)]+const_;
							}
							break;
					}	//if there is an to the right of dim 'i' then the BC is calculated to the LEFT
				break;
			}

		}

		template<class Shape_t, int Dim>
		void applyBC(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{
			if(alter.size() != test.size()) SizeErr();

			typename BoundaryGeneral::iterator myit(*this);
			while(myit)
			{
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


//second order extrapolation--> (fit a line...needs near and next near neighbor)
template<class Num_t>
class ExtrapolateBC<2, Num_t> :
	public BoundaryGeneral
{
	private:
		Num_t const_;
		const void SizeErr()
		{
			BLEXCEPTION(" Data size and StencilPrep size are NOT the same...")
		}
	public:

		typedef Num_t numtype;

		ExtrapolateBC():
			BoundaryGeneral(),
			const_(1)
		{}

		ExtrapolateBC(const numtype &in):
			BoundaryGeneral(),
			const_(in)
		{}

		template<class Shape_t, int Dim>
		ExtrapolateBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph,
			int opts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph, opts),
			const_(in)
		{}

		template<class Shape_t, int Dim, class Shapeexpr>
		ExtrapolateBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph, const ShapeExpr<Shapeexpr> &se,
			int opts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph,se,opts),
			const_(in)
		{}


		~ExtrapolateBC(){ 	}

		template<class Shape_t, int Dim>
		void applyBCOnce(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter, typename BoundaryGeneral::iterator &myit)
		{
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
										alter[myit()]=2.0*alter[test(myit()).near(1,i)]-alter[test(myit()).nextnear(2,i)] + const_;
									}else if( test(myit()).near(1,i) !=-1){
										alter[myit()]=alter[test(myit()).near(1,i)] + const_;
									}
								//if there is an to the right of dim 'i' then the BC is calculated to the LEFT
								}else if(test(myit()).edge.IsFace(1,i)){
									if(test(myit()).nextnear(-2,i)!=-1 && test(myit()).near(-1,i) !=-1){
										alter[myit()]=2.0*alter[test(myit()).near(-1,i)]-alter[test(myit()).nextnear(-2,i)] + const_;
									}else if( test(myit()).near(-1,i) !=-1){
										alter[myit()]=alter[test(myit()).near(1,i)]+ const_;
									}
								}
								break;
							case BCparams::PositiveFaceNorm:
							case BCparams::NegativeFaceNorm:
							default:
								//if(test(myit()).edge.IsFace(-facenormals,i))
								//{
								if(test(myit()).nextnear(-2*facenormals,i)!=-1 && test(myit()).near(-facenormals,i) !=-1){
									alter[myit()]=2.0*alter[test(myit()).near(-facenormals,i)]-alter[test(myit()).nextnear(-2*facenormals,i)] + const_;
								}else if( test(myit()).near(-facenormals,i) !=-1){
									alter[myit()]=alter[test(myit()).near(-facenormals,i)]+ const_;
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
									alter[myit()]=2.0*alter[test(myit()).near(1,faces)]-alter[test(myit()).nextnear(2,faces)] + const_;
								}else if( test(myit()).near(1,faces) !=-1){
									alter[myit()]=alter[test(myit()).near(1,faces)]+ const_;
								}
							//if there is an to the right of dim 'i' then the BC is calculated to the LEFT
							}else if(test(myit()).edge.IsFace(1,faces)){
								if(test(myit()).nextnear(-2,faces)!=-1 && test(myit()).near(-1,faces) !=-1){
									alter[myit()]=2.0*alter[test(myit()).near(-1,faces)]-alter[test(myit()).nextnear(-2,faces)] + const_;
								}else if( test(myit()).near(-1,faces) !=-1){
									alter[myit()]=alter[test(myit()).near(-1,faces)] + const_;
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
									alter[myit()]=2.0*alter[test(myit()).near(-facenormals,faces)]-alter[test(myit()).nextnear(-2*facenormals,faces)] + const_;
								}else if( test(myit()).near(-facenormals,faces) !=-1){
									alter[myit()]=alter[test(myit()).near(-facenormals,faces)]+ const_;
								}
							//}
							break;
					}	//if there is an to the right of dim 'i' then the BC is calculated to the LEFT
				break;
			}
		}

		template<class Shape_t, int Dim>
		void applyBC(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{
			if(alter.size() != test.size()) SizeErr();

			typename BoundaryGeneral::iterator myit(*this);
			int i=0;
			while(myit)
			{
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


