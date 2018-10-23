
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

 	boundary_dirichlet.h-> implimentation for the 'contstant' BC (called the
 	dirichelt BC) holds the face fixed at a constant value.

 */


#ifndef _boundary_dirchlet_h_
#define _boundary_dirchlet_h_ 1

#include "stencils/stencilprep.h"
#include "stencils/boundary_condition.h"
#include "stencils/boundary_general.h"

#include "container/grids/edgedetect.h"
#include "container/grids/xyzshape.h"



BEGIN_BL_NAMESPACE




/*********** Constant Numbers.....

there are 2 defined here...one for 'coord<>' and one for 'double'

***************/
template<class  Numt >
class ConstantBC;

template<>
class ConstantBC<coord<> > :
	public BoundaryGeneral
{
	private:
		coord<> const_;
		const void SizeErr()
		{
			BLEXCEPTION(" Data size and StencilPrep size are NOT the same...")
		}
	public:

		typedef coord<> numtype;

		ConstantBC():
			BoundaryGeneral(),
			const_(0)
		{}

		ConstantBC(const numtype &in):
			BoundaryGeneral(),
			const_(in)
		{}

		template<class Shape_t, int Dim>
		ConstantBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph,inopts),
			const_(in)
		{}

		template<class Shape_t, int Dim, class Shapeexpr>
		ConstantBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph, const ShapeExpr<Shapeexpr> &se,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph,se,inopts),
			const_(in)
		{}


		~ConstantBC(){ 	}

		template<class Shape_t, int Dim>
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
		}

		template<class Shape_t, int Dim>
		inline void apply(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{	applyBC(test, alter);	}


		template<class Shape_t, int Dim>
		void operator()(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{	applyBC(test, alter);	}

//apply the BC from the



};

class ConstantBC<double > :
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

		ConstantBC():
			BoundaryGeneral(),
			const_(0)
		{}

		ConstantBC(const numtype &in):
			BoundaryGeneral(),
			const_(in)
		{}

		template<class Shape_t, int Dim>
		ConstantBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph,inopts),
			const_(in)
		{}

		template<class Shape_t, int Dim, class Shapeexpr>
		ConstantBC(const numtype &in,StencilPrep<Shape_t, Dim> &sph, const ShapeExpr<Shapeexpr> &se,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BoundaryGeneral(sph,se,inopts),
			const_(in)
		{}


		~ConstantBC(){ 	}

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
		void applyBCOnce(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter, typename BoundaryGeneral::iterator &myit)
		{
				alter[myit()]=const_;
		}

		template<class Shape_t, int Dim>
		inline void apply(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{	applyBC(test, alter);	}


		template<class Shape_t, int Dim>
		void operator()(StencilPrep<Shape_t, Dim> &test,  Vector<numtype> &alter)
		{	applyBC(test, alter);	}


};


/***************** END CONSTANT *********************/


END_BL_NAMESPACE

#endif

