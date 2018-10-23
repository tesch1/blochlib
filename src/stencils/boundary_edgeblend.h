


/* boundary_condition.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-10-01
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

 	boundary_edgeblend.h-> takes an edge and corner of an object
 	and 'blends' the vales of its neighbors...

 	since mose BCs at edges and corners of shapes are ill defined
 	this class makes the assumption that at an edge you want a 'equal weighting'
 	along both faces (corners get equal weighting along the 3 directions)

 	of course you can specify which edges and such

 */


#ifndef _boundary_edgeblend_h_
#define _boundary_edgeblend_h_ 1

#include "stencils/stencilprep.h"
#include "stencils/boundary_condition.h"
#include "stencils/boundary_general.h"
#include "stencils/boundary_cornerblend.h"

#include "container/grids/edgedetect.h"
#include "container/grids/xyzshape.h"


BEGIN_BL_NAMESPACE



/*********** Constant Numbers.....

there are 2 defined here...one for 'coord<>' and one for 'double'

***************/
//template<class  Numt >
//class EdgeBlendBC;

template<class Num_t>
class EdgeBlendBC :
	public BoundaryGeneral
{
	private:
		CornerBlendBC<Num_t > cb;
		Num_t const_;
		const void SizeErr()
		{
			BLEXCEPTION(" Data size and StencilPrep size are NOT the same...")
		}


	public:

		typedef Num_t numtype;

		EdgeBlendBC():
			BoundaryGeneral(),
			const_(0)
		{}

		EdgeBlendBC(const numtype &in):
			BoundaryGeneral(),
			const_(in)
		{}

		template<class Shape_t, int Dim>
		EdgeBlendBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms ):
			BoundaryGeneral(sph,(!(inopts&BCparams::EdgesOnly)?(inopts | BCparams::EdgesOnly):inopts)),
			cb(in, sph, inopts),
			const_(in)
		{}

		template<class Shape_t, int Dim, class Shapeexpr>
		EdgeBlendBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph, const ShapeExpr<Shapeexpr> &se,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms):
			BoundaryGeneral(sph,se,(!(inopts&BCparams::EdgesOnly)?(inopts | BCparams::EdgesOnly):inopts) ),
			cb(in, sph, se,inopts),
			const_(in)
		{}


		~EdgeBlendBC(){ 	}

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
			cb.applyBC(test, alter);
		}

/*		template<class Shape_t, int Dim>
		void applyBC(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{
			if(alter.size() != test.size()) SizeErr();

			typename BoundaryGeneral::iterator myit(*this);
			while(myit)
			{
				for(int i=0;i<Dim;++i){
					if(test.edge(myit()).IsFace(i)){	alter[myit()][i]=const_[i]; 	}
				//	cout<<myit.edge()<<"|| "<<myit.edge().IsFace(i) <<" || " <<alter(myit.curpos())(i)<<endl;
				}
				++myit;
			}
			cb.applyBC(test, alter);
		}
*/
		template<class Shape_t, int Dim>
		void applyBCOnce(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter, typename BoundaryGeneral::iterator &myit)
		{
			//the 'edge selection' should have already been made in the 'BoundaryGeneral' class



			int i=0;
			//if(!test(myit()).edge.IsCorner()){
				EdgePos<Dim-1> mye=test(myit()).edge.EdgeIdx();

				if(mye.idx != -1)
				{
					double we=1./2.;
					alter[myit()]=0.0;
					for(i=0;i<Dim-1;++i)
					{ alter[myit()]+=we*alter[test(myit()).near(mye.dir[i], mye.idx[i])];	}
					return;
				}

				//first detect if it is a corner
			//}
		}

		template<class Shape_t, int Dim>
		inline void apply(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{	applyBC(test, alter);	}


		template<class Shape_t, int Dim>
		void operator()(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{	applyBC(test, alter);	}

//apply the BC from the



};



END_BL_NAMESPACE



/***************** END EdgeBlend *********************/

#endif

