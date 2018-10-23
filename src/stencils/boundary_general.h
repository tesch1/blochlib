


/* boundary_general.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-16-01
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

 	boundary_general.h-> given a shape/shape expression AND a previously calculated
 	'StencilPrep' object will create a list of ints that are the indexes of the
 	'edge' points from the stencil prep that are also in that shape..

 	this class should not be used directly...instead it is a base class for the
 	'specific' boundary conditions (i.e. the ones that apply the BC)
 */


#ifndef _boundary_general_h_
#define _boundary_general_h_ 1

#include "stencils/stencilprep.h"


#include "container/grids/edgedetect.h"
#include "container/grids/xyzshape.h"


BEGIN_BL_NAMESPACE

class BCparams{
	public:
		enum Faces
		{
			Xface=0x00001,
			Yface=0x00002,
			Zface=0x00004,
			AllFaces=0x00008,
			NoFaces=0x00010
		};

		enum FacesInt
		{
			XfaceInt=0,
			YfaceInt=1,
			ZfaceInt=2,
			AllFacesInt=-1,
			NoFacesInt=10
		};

		enum FaceNormals
		{
			PositiveFaceNorm=0x00100,
			NegativeFaceNorm=0x00200,
			AllFaceNorms=0x00400
		};

		enum FaceNormalsInt
		{
			PositiveFaceNormInt=1,
			NegativeFaceNormInt=-1,
			AllFaceNormsInt=0
		};

		enum Edges
		{
			NoEdges=0x01000,
			WithEdges=0x02000,
			EdgesOnly=0x04000,
			CornersOnly=0x08000
		};

		enum EdgesInt
		{
			NoEdgesInt=0,
			WithEdgesInt=1,
			EdgesOnlyInt=2,
			CornersOnlyInt=3
		};

		enum Interior
		{
			WithInterior=0x10000,
			NoInterior=0x20000
		};


};


class BoundaryGeneralIter;

class BoundaryGeneral {
	private:
		Vector<int> stpidx_;

	public:

		int facenormals;
		int faces;
		int edges;
		int interior;

		typedef BoundaryGeneralIter iterator;

		friend class BoundaryGeneralIter;

		BoundaryGeneral():
			stpidx_(0)
		{}

	//construcutor assuming ONLY the StPrep object (i.e. uses those edges...and
	// does not attmpt to calculate more
		template<class Shape_t, int Dim>
		BoundaryGeneral(StencilPrep<Shape_t, Dim> &stp,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges | BCparams::NoInterior)
		{
			setOptions(inopts);
			calculate(stp);
		}

	//construcutor assuming BOTH the StPrep object AND a shape/shapeExpr
		template<class Shape_t, int Dim, class Shapeexpr>
		BoundaryGeneral(StencilPrep<Shape_t, Dim> &stp, Shapeexpr test,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges | BCparams::NoInterior)
		{
			setOptions(inopts);
			calculate(stp, test);
		}


	//set the options...
		void setOptions(int inops)
		{
			if(inops & BCparams::Xface){ faces=BCparams::XfaceInt;	}
			else if(inops & BCparams::Yface){	 faces=BCparams::YfaceInt;	}
			else if(inops & BCparams::Zface){ faces=BCparams::ZfaceInt;	}
			else if(inops & BCparams::NoFaces){ faces=BCparams::NoFacesInt;	}
			else{	faces=BCparams::AllFacesInt;	}

			if(inops & BCparams::PositiveFaceNorm){ facenormals=BCparams::PositiveFaceNormInt;	}
			else if(inops & BCparams::NegativeFaceNorm){	 facenormals=BCparams::NegativeFaceNormInt;	}
			else{ facenormals=BCparams::AllFaceNormsInt;	}

			if(inops & BCparams::WithEdges){	 edges=BCparams::WithEdgesInt;	}
			else if(inops & BCparams::EdgesOnly) { edges=BCparams::EdgesOnlyInt;	}
			else if(inops & BCparams::CornersOnly) { edges=BCparams::CornersOnlyInt;	}
			else{ edges=BCparams::NoEdgesInt;	}

			if(inops & BCparams::WithInterior){	interior=BCparams::WithInterior;	}
			else{	interior=BCparams::NoInterior;	}
		}


	//calculate assumping we use ALL the edges in the stencil prep object
		template<class Shape_t, int Dim>
		void calculate(StencilPrep<Shape_t, Dim> &stp)
		{
			typename StencilPrep<Shape_t, Dim>::iterator myit(stp);
			while(myit)
			{
				//get the interios if need be
				switch(interior)
				{
					case BCparams::WithInterior:
						if(!myit.edge().IsFace()) stpidx_.push_back(myit.curpos());
						break;
					case BCparams::NoInterior:
					default:
						break;
				}

				//detect faces..so that these get hit first
				switch(faces)
				{
					case BCparams::AllFacesInt:
						switch(facenormals)
						{
							case BCparams::AllFaceNormsInt:
								switch(edges)
								{
									case BCparams::WithEdgesInt:
										if(myit.edge().IsFace()) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::EdgesOnlyInt:
										if(myit.edge().IsEdge()) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::CornersOnlyInt:
										if(myit.edge().IsCorner()) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::NoEdgesInt:
									default:
										if(myit.edge().IsFace() && !myit.edge().IsEdge()) stpidx_.push_back(myit.curpos());
										break;
								}
								break;
							case BCparams::PositiveFaceNormInt:
							case BCparams::NegativeFaceNormInt:
							default:
								switch(edges)
								{
									case BCparams::WithEdgesInt:
										for(int i=0;i<Dim;++i)
										{
											if(myit.edge().IsFace(facenormals,i)) stpidx_.push_back(myit.curpos());
										}
										break;
									case BCparams::EdgesOnlyInt:
										for(int i=0;i<Dim;++i)
										{
											if(myit.edge().IsEdge(facenormals,i)) stpidx_.push_back(myit.curpos());
										}
										break;
									case BCparams::CornersOnlyInt:
										if(myit.edge().IsCorner(facenormals)) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::NoEdgesInt:
									default:
										for(int i=0;i<Dim;++i)
										{
											if(myit.edge().IsFace(facenormals, i) && !myit.edge().IsEdge()) stpidx_.push_back(myit.curpos());
										}
										break;
								}
								break;
						}
						break;
					case BCparams::NoFacesInt: break;
					case BCparams::XfaceInt:
					case BCparams::YfaceInt:
					case BCparams::ZfaceInt:
					default:
						switch(facenormals)
						{
							case BCparams::AllFaceNormsInt:
								switch(edges)
								{
									case BCparams::WithEdgesInt:
										if(myit.edge().IsFace(facenormals,faces)) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::EdgesOnlyInt:
										if(myit.edge().IsEdge(facenormals,faces)) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::CornersOnlyInt:
										if(myit.edge().IsFace(facenormals,faces) && myit.edge().IsCorner()) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::NoEdgesInt:
									default:
										if(myit.edge().IsFace(facenormals, faces) && !myit.edge().IsEdge()) stpidx_.push_back(myit.curpos());
										break;
								}
								break;
							case BCparams::PositiveFaceNormInt:
							case BCparams::NegativeFaceNormInt:
							default:
								switch(edges)
								{
									case BCparams::WithEdgesInt:
										if(myit.edge().IsFace(facenormals,faces)) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::EdgesOnlyInt:
										if(myit.edge().IsEdge(facenormals,faces)) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::CornersOnlyInt:
										if(myit.edge().IsFace(facenormals,faces) && myit.edge().IsCorner()) stpidx_.push_back(myit.curpos());
										break;
									case BCparams::NoEdgesInt:
									default:
										if(myit.edge().IsFace(facenormals, faces) && !myit.edge().IsEdge()) stpidx_.push_back(myit.curpos());
										break;
								}
								break;
						}
				}
				++myit;
			}
		}


	//calculate assumping we use ALL the edges in the stencil prep object
		template<class Shape_t, int Dim, class Shapeexpr>
		void calculate(StencilPrep<Shape_t, Dim> &stp,const ShapeExpr<Shapeexpr> &test)
		{
			typename StencilPrep<Shape_t, Dim>::iterator myit(stp);
			while(myit)
			{
				switch(faces)
				{
					case BCparams::AllFaces:
						switch(facenormals)
						{
							case BCparams::AllFaceNorms:
								if(myit.edge().IsFace()  && test.ShapeFunc(myit.GridPoint())) stpidx_.push_back(myit.curpos());
								break;
							case PositiveFaceNorm:
							case NegativeFaceNorm:
							default:
								for(int i=0;i<Dim;++i)
								{
									if(myit.edge().IsFace(facenormals, i)  && test.ShapeFunc(myit.GridPoint())) stpidx_.push_back(myit.curpos());
								}
								break;
						}
						break;
					case BCparams::Xface:
					case BCparams::Yface:
					case BCparams::Zface:
					default:
						switch(facenormals)
						{
							case BCparams::AllFaceNorms:
								if(myit.edge().IsFace(faces)  && test.ShapeFunc(myit.GridPoint())) stpidx_.push_back(myit.curpos());
								break;
							case BCparams::PositiveFaceNorm:
							case BCparams::NegativeFaceNorm:
							default:
								if(myit.edge().IsFace(facenormals,faces)  && test.ShapeFunc(myit.GridPoint())) stpidx_.push_back(myit.curpos());
								break;
						}
				}
				++myit;
			}
		}

		inline int size() const {	return stpidx_.size();	}

		inline int operator()(int i) const	{	return stpidx_(i);	}

		inline Vector<int> &data(){	return stpidx_;	}
		inline Vector<int> data()	const	{	return stpidx_;	}

		inline int &data(int i){	return stpidx_(i);	}
		inline int data(int i)	const	{	return stpidx_(i);	}
		inline void SetFace(int in) {	faces=in; }
		inline void SetFaceNormals(int in){	facenormals=in;	}

};



class BoundaryGeneralIter{
	private:
		BoundaryGeneral *bg_;
		bool notended_;
		int curpos;
		int end;
	public:
		BoundaryGeneralIter():
			bg_(NULL),
			notended_(false), curpos(0), end(0)
		{}

		BoundaryGeneralIter( BoundaryGeneral &in):
			bg_(&in),
			notended_((in.size())?true:false), curpos(0), end(in.size())
		{}

		BoundaryGeneralIter( BoundaryGeneral *in):
			bg_(in),
			notended_((in->size())?true:false), curpos(0), end(in->size())
		{}

		~BoundaryGeneralIter(){	bg_=NULL;	}


		inline operator bool(){	return notended_;	}

		void operator++()
		{	if(curpos<end-1){	++curpos;}else{ notended_=false;}	}

		void operator++(int)
		{	operator++();	}

		int operator()(){	return bg_->data(curpos);	}

		void reset(){	curpos=0; notended_=true;	}

		inline int size() const {	return bg_->size();	}
};



END_BL_NAMESPACE

#endif
