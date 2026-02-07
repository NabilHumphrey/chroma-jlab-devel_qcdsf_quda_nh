/*! \file
 *  \brief QUDA solver for UNPRECONDITIONED_SLRC fermion action.
 *
 *  Registration and solve implementation.
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_slrc_quda_w.h"
#include "io/aniso_io.h"

#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <quda.h>
#include "actions/ferm/invert/quda_solvers/quda_mg_utils.h"

namespace Chroma
{
  namespace LinOpSysSolverQUDASLRCEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_SLRC_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }

    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverQUDASLRC(A, state, SysSolverQUDAMULTIGRIDCloverParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

  SystemSolverResults_t
  LinOpSysSolverQUDASLRC::qudaInvert(const CloverTermT<T, U>& clover,
				      const CloverTermT<T, U>& invclov,
				      const T& chi_s,
				      T& psi_s) const
  {

    SystemSolverResults_t ret;

    // For the unpreconditioned SLRC, we pass full-field spinors
    // QUDA handles internal even-odd preconditioning via DIRECT_PC_SOLVE + MAT_SOLUTION

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    void* spinorIn  = (void *)&(chi_s.elem(all.start()).elem(0).elem(0).real());
    void* spinorOut = (void *)&(psi_s.elem(all.start()).elem(0).elem(0).real());
#else
    void* spinorIn;
    void* spinorOut;
    std::vector<QDPCache::ArgKey> ids;
    ids.push_back(chi_s.getId());
    ids.push_back(psi_s.getId());
    auto dev_ptr = GetMemoryPtr(ids);
    spinorIn  = dev_ptr[0];
    spinorOut = dev_ptr[1];
#endif

    // Do the solve here
    StopWatch swatch1;
    swatch1.reset();
    swatch1.start();
    invertQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
    swatch1.stop();

    QDPIO::cout << solver_string << "Total Time (incl. load gauge)=" << swatch1.getTimeInSeconds() << " s" << std::endl;

    ret.n_count = quda_inv_param.iter;
    ret.resid = quda_inv_param.true_res[0];
    return ret;

  }

}
