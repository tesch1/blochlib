
/* Master Header file for the bloch lib */

#ifndef _blochlib_h_
#define _blochlib_h_ 1


//this is defined in the CodeWarrior include
#ifdef ON_WINDOWS
 #include "blochconfigCW.h"
#else
 #include "blochconfig.h"
#endif

/*#ifdef USE_MINUIT
#ifdef __cplusplus
extern "C"{
#endif
	#include "minuit/cfortran.h"
	#include "minuit/minuitfcn.h"
	#include "minuit/minuit.h"
#ifdef __cplusplus
}
#endif
#endif
*/

#ifdef HAVE_MPI
	#include "mpi.h"
#endif

#include "mpi/mpi_config.h"
#include "mpi/mpi_packer.h"
#include "mpi/mpi_tools.h"
#include "mpi/mpi_controller.h"


#ifdef HAVE_FFTW
	#include "fftw.h"
#endif

#include "version.h"
#include "utils/blassert.h"
#include "utils/plotter.h"
#include "utils/random.h"
#include "utils/utils.h"
#include "utils/constants.h"
#include "utils/params.h"
#include "utils/matlab5.h"
#include "utils/vnmrstream.h"
#include "utils/xwinnmr.h"
#include "utils/wavestream.h"
#include "utils/spinsight.h"
#include "utils/endians.h"
#include "utils/parser.h"
#include "utils/scriptparse.h"



#include "container/range.h"
#include "container/rankType.h"
#include "container/MemChunk.h"
#include "container/operations.h"
#include "container/complex.h"
#include "container/grids/coords.h"
#include "container/Vector/Vector.h"

#include "container/matrix/matrix.h"

#include "container/grids/grids.h"
#include "container/grids/fullcubegrid.h"
#include "container/grids/fullspheregrid.h"
#include "container/grids/surfcubegrid.h"
#include "container/grids/surfspheregrid.h"
#include "container/grids/unigrid.h"
#include "container/grids/gengrid.h"

#include "container/grids/xyzshape.h"
#include "container/grids/edgedetect.h"

#include "driver/bs.h"
#include "driver/stiffbs.h"
#include "driver/ckrk.h"
#include "driver/gear.h"
#include "driver/integrate.h"

#include "QMspins/spinsys.h"
#include "QMspins/spinsyspars.h"
#include "QMspins/singlespinop.h"
#include "QMspins/multispinop.h"
#include "QMspins/spin_ten.h"
#include "QMspins/space_ten.h"
#include "QMspins/j.h"
#include "QMspins/csa.h"
#include "QMspins/dip.h"
#include "QMspins/qua.h"
#include "QMspins/sys.h"

#include "QMspins/driver/powder.h"
#include "QMspins/driver/ham_gen.h"
#include "QMspins/driver/compute.h"
#include "QMspins/driver/oneFID.h"

#include "timetrain/timetrain.h"
#include "timetrain/unitrain.h"
#include "timetrain/consttrain.h"
#include "timetrain/filetrain.h"


#include "stencils/stencilexpr.h"
#include "stencils/stencilprep.h"
#include "stencils/stencilextent.h"
#include "stencils/stencilvector.h"
#include "stencils/stencilprep_func.h"
#include "stencils/boundary_condition.h"
#include "stencils/boundary_dirichlet.h"
#include "stencils/boundary_neumann.h"
#include "stencils/boundary_extrap.h"
#include "stencils/boundary_edgeblend.h"
#include "stencils/boundary_cornerblend.h"

#include "bloch/Isotope.h"
#include "bloch/blochParamsBasic.h"
#include "bloch/blochParams.h"
#include "bloch/pulse.h"
#include "bloch/listblochpars.h"
#include "bloch/gradientgrid.h"
#include "bloch/rotating_grid.h"
#include "bloch/bloch_basic.h"
#include "bloch/bloch.h"
#include "bloch/blochsolver.h"
#include "bloch/bloch_interac.h"
#include "bloch/bloch_int_dipdip.h"
#include "bloch/bloch_int_raddamp.h"
#include "bloch/bloch_int_bulks.h"
#include "bloch/bloch_int_relax.h"

#include "bloch/bloch_int_offset.h"
#include "bloch/biot.h"
#include "bloch/mag_fields.h"
#include "bloch/bloch_int_bfield.h"
#include "bloch/bloch_int_demag.h"
#include "bloch/bloch_int_dimlessdemag.h"
#include "bloch/bloch_int_demagmodulate.h"
//#include "bloch/bloch_int_diffu.h"
#include "bloch/bloch_int_dipdip.h"
//#include "bloch/bloch_int_qmdip.h"
#include "bloch/bloch_int_dimlessdip.h"
#include "bloch/scalefunc.h"
#include "bloch/gradfunc.h"
#include "bloch/gradimps.h"


#endif

