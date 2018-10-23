

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

 	boundary_condition.h-> a series of classes that apply boundry conditions
 	given a 'stencilprep' and/or a XYZshpae expression to determin
 	the edges
 */


#ifndef _boundary_condition_h_
#define _boundary_condition_h_ 1

#include "stencils/stencilprep.h"
#include "stencils/boundary_general.h"

#include "container/grids/edgedetect.h"
#include "container/grids/xyzshape.h"


BEGIN_BL_NAMESPACE


//the master container for a give boundary condition
template<class BCeng_t>
class BoundaryCondition:
	public BCeng_t
{
	public:

		typedef typename BCeng_t::numtype numtype;

		BoundaryCondition():
			BCeng_t()
		{}

		BoundaryCondition(const numtype &in):
			BCeng_t(in)
		{}

		template<class Shape_t, int Dim>
		BoundaryCondition(const numtype &in,StencilPrep<Shape_t, Dim> &sph,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BCeng_t(in,sph,inopts)
		{}

		template<class Shape_t, int Dim, class Shapeexpr>
		BoundaryCondition(const numtype &in,StencilPrep<Shape_t, Dim> &sph, const ShapeExpr<Shapeexpr> &se,
			int inopts=BCparams::AllFaces | BCparams::AllFaceNorms | BCparams::NoEdges):
			BCeng_t(in,sph,se,inopts)
		{}
};




END_BL_NAMESPACE

#endif


