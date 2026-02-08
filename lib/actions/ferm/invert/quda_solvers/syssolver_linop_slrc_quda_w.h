// -*- C++ -*-
/*! \file
 *  \brief QUDA solver for SLRC fermion action.
 *
 *  Handles the SLRC decomposition by:
 *  1. Extracting thin (unsmeared) links from SLICFermState for the clover term
 *  2. Extracting fat (stout-smeared) links for the hopping term
 *  3. Computing the clover term on the Chroma side using thin links
 *  4. Passing the fat links as the gauge field to QUDA
 *  5. Passing the pre-computed clover term to QUDA via loadCloverQuda()
 *
 *  Uses MATPC_SOLUTION (odd-checkerboard) for QDPJIT compatibility.
 *  Supports multigrid preconditioning with subspace caching.
 */

#ifndef __syssolver_linop_quda_slrc_h__
#define __syssolver_linop_quda_slrc_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA
#include <quda.h>

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include "quda_mg_utils.h"
#include <string>
#include <sstream>
#include "util/gauge/reunit.h"
#ifdef QDP_IS_QDPJIT
#include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#endif

namespace Chroma
{

	//! SLRC QUDA solver namespace
	namespace LinOpSysSolverQUDASLRCEnv
	{
		//! Register the syssolver
		bool registerAll();
	}

	//! Solve an SLRC Fermion System using the QUDA inverter
	/*! \ingroup invert
	 *** This solver is for UNPRECONDITIONED_SLRC fermions with SLIC_FERM_STATE ***
	 *
	 *  The SLRC action uses different gauge fields for the hopping term
	 *  (stout-smeared / fat links) and the clover term (unsmeared / thin links).
	 *  This solver handles the decomposition by pre-computing the clover term
	 *  from thin links on the Chroma side and loading it separately to QUDA,
	 *  while the fat links are loaded as the gauge field for the hopping term.
	 */

	class LinOpSysSolverQUDASLRC : public LinOpSystemSolver<LatticeFermion>
	{
	public:
		typedef LatticeFermion T;
		typedef LatticeColorMatrix U;
		typedef multi1d<LatticeColorMatrix> Q;

		typedef LatticeFermionF TF;
		typedef LatticeColorMatrixF UF;
		typedef multi1d<LatticeColorMatrixF> QF;

		typedef LatticeFermionF TD;
		typedef LatticeColorMatrixF UD;
		typedef multi1d<LatticeColorMatrixF> QD;

		typedef WordType<T>::Type_t REALT;

		//! Constructor
		/*!
		 * \param A_        Linear operator (UnprecSLRCLinOp)
		 * \param state_    The SLIC ferm state (contains both thin and fat links)
		 * \param invParam_ Solver parameters (reuses SysSolverQUDAMULTIGRIDCloverParams)
		 */
		LinOpSysSolverQUDASLRC(Handle< LinearOperator<T> > A_,
				Handle< FermState<T,Q,Q> > state_,
				const SysSolverQUDAMULTIGRIDCloverParams& invParam_) :
		A(A_), invParam(invParam_), clov(new CloverTermT<T, U>() ), invclov(new CloverTermT<T, U>())
		{
			StopWatch init_swatch;
			init_swatch.reset(); init_swatch.start();

			// Set the solver string
			{
				std::ostringstream solver_string_stream;
				solver_string_stream << "QUDA_SLRC_SOLVER( "
				    << invParam.SaveSubspaceID << " ): ";
				solver_string = solver_string_stream.str();
			}
			QDPIO::cout << solver_string << "Initializing SLRC+QUDA solver" << std::endl;

			// 1) Work out cpu_prec, cuda_prec, cuda_prec_sloppy
			int s = sizeof( WordType<T>::Type_t );
			if (s == 4) {
				cpu_prec = QUDA_SINGLE_PRECISION;
			}
			else {
				cpu_prec = QUDA_DOUBLE_PRECISION;
			}

			// Work out GPU precision
			switch( invParam.cudaPrecision ) {
				case HALF:
				gpu_prec = QUDA_HALF_PRECISION;
				break;
				case SINGLE:
				gpu_prec = QUDA_SINGLE_PRECISION;
				break;
				case DOUBLE:
				gpu_prec = QUDA_DOUBLE_PRECISION;
				break;
				default:
				gpu_prec = cpu_prec;
				break;
			}

			// Work out GPU Sloppy precision
			switch( invParam.cudaSloppyPrecision ) {
				case HALF:
				gpu_half_prec = QUDA_HALF_PRECISION;
				break;
				case SINGLE:
				gpu_half_prec = QUDA_SINGLE_PRECISION;
				break;
				case DOUBLE:
				gpu_half_prec = QUDA_DOUBLE_PRECISION;
				break;
				default:
				gpu_half_prec = gpu_prec;
				break;
			}

			// 2) Create GAUGE and Invert params
			q_gauge_param = newQudaGaugeParam();
			quda_inv_param = newQudaInvertParam();

			// 3) Set lattice size
			const multi1d<int>& latdims = Layout::subgridLattSize();
			q_gauge_param.X[0] = latdims[0];
			q_gauge_param.X[1] = latdims[1];
			q_gauge_param.X[2] = latdims[2];
			q_gauge_param.X[3] = latdims[3];

			// 5) Set QUDA_WILSON_LINKS, QUDA_GAUGE_ORDER
			q_gauge_param.type = QUDA_WILSON_LINKS;
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
			q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;
#else
			q_gauge_param.location = QUDA_CUDA_FIELD_LOCATION;
			q_gauge_param.gauge_order = QUDA_QDPJIT_GAUGE_ORDER;
#endif

			// 6) Set t_boundary
			if( invParam.AntiPeriodicT ) {
				q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
			}
			else {
				q_gauge_param.t_boundary = QUDA_PERIODIC_T;
			}

			// Set cpu_prec, cuda_prec, reconstruct and sloppy versions
			q_gauge_param.cpu_prec = cpu_prec;
			q_gauge_param.cuda_prec = gpu_prec;

			switch( invParam.cudaReconstruct ) {
				case RECONS_NONE:
				q_gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
				break;
				case RECONS_8:
				q_gauge_param.reconstruct = QUDA_RECONSTRUCT_8;
				break;
				case RECONS_12:
				q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
				break;
				default:
				q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
				break;
			};

			q_gauge_param.cuda_prec_sloppy = gpu_half_prec;
			q_gauge_param.cuda_prec_precondition = gpu_half_prec;

			switch( invParam.cudaSloppyReconstruct ) {
				case RECONS_NONE:
				q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
				break;
				case RECONS_8:
				q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_8;
				break;
				case RECONS_12:
				q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
				break;
				default:
				q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
				break;
			};
			q_gauge_param.reconstruct_precondition = q_gauge_param.reconstruct_sloppy;

			// ============================================================
			// SLRC-specific: Extract thin and fat links from SLIC state
			// ============================================================
			QDPIO::cout << solver_string
			            << "Extracting thin and fat links from SLICFermState"
			            << std::endl;

			// Fat links: for hopping term / QUDA gauge field
			Q fat_links(Nd);
			for(int mu=0; mu < Nd; mu++) {
				fat_links[mu] = (state_->getLinks())[mu];
			}

			// Thin links: for clover term computation
			Q thin_links;
			try {
				const StoutFermState<T,Q,Q>& slic_state =
					dynamic_cast<const StoutFermState<T,Q,Q>&>(*state_);
				thin_links = slic_state.getThinLinks();
			}
			catch (std::bad_cast& e) {
				QDPIO::cerr << solver_string
				            << "ERROR: FermState is not a StoutFermState/SLICFermState! "
				            << "This solver requires SLIC_FERM_STATE."
				            << std::endl;
				QDP_abort(1);
			}

			// Diagnostic: verify thin and fat links differ
			{
				Double thin_plaq = zero;
				Double fat_plaq = zero;
				for (int mu = 1; mu < Nd; mu++) {
					for (int nu = 0; nu < mu; nu++) {
						thin_plaq += sum(real(trace(
							thin_links[mu] * shift(thin_links[nu], FORWARD, mu)
							* adj(shift(thin_links[mu], FORWARD, nu)) * adj(thin_links[nu])
						)));
						fat_plaq += sum(real(trace(
							fat_links[mu] * shift(fat_links[nu], FORWARD, mu)
							* adj(shift(fat_links[mu], FORWARD, nu)) * adj(fat_links[nu])
						)));
					}
				}
				thin_plaq /= Double(Layout::vol() * Nc * Nd * (Nd-1) / 2);
				fat_plaq  /= Double(Layout::vol() * Nc * Nd * (Nd-1) / 2);

				QDPIO::cout << solver_string << "Thin link plaquette = " << thin_plaq << std::endl;
				QDPIO::cout << solver_string << "Fat  link plaquette = " << fat_plaq << std::endl;

				if ( toBool( fabs(thin_plaq - fat_plaq) < 1.0e-12 ) ) {
					QDPIO::cerr << solver_string
					            << "WARNING: Thin and fat link plaquettes are identical! "
					            << "Check that SLIC_FERM_STATE is configured correctly."
					            << std::endl;
				}
			}

			// GaugeFix (applied to fat links for QUDA)
			if( invParam.axialGaugeP ) {
				temporalGauge(fat_links, GFixMat, Nd-1);
				for(int mu=0; mu < Nd; mu++) {
					fat_links[mu] = GFixMat*(state_->getLinks())[mu]*adj(shift(GFixMat, FORWARD, mu));
				}
				q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_YES;
			}
			else {
				q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
			}

			// Gauge Anisotropy
			const AnisoParam_t& aniso = invParam.CloverParams.anisoParam;
			if( aniso.anisoP ) {
				Real gamma_f = aniso.xi_0 / aniso.nu;
				q_gauge_param.anisotropy = toDouble(gamma_f);
			}
			else {
				q_gauge_param.anisotropy = 1.0;
			}

			// ============================================================
			// SLRC-specific: Compute clover term from THIN links
			// ============================================================
			// Create FermState from thin links BEFORE anisotropy rescaling
			Handle<FermState<T,Q,Q> > thin_fstate(
				new PeriodicFermState<T,Q,Q>(thin_links));

			// Apply anisotropy rescaling to fat links (for QUDA gauge field)
			if( aniso.anisoP ) {
				multi1d<Real> cf = makeFermCoeffs(aniso);
				for(int mu=0; mu < Nd; mu++) {
					fat_links[mu] *= cf[mu];
				}
			}

			// Now onto the inv param:
			quda_inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;

			// Solver type: GCR for MG-preconditioned, configurable otherwise
			if( invParam.MULTIGRIDParamsP ) {
				quda_inv_param.inv_type = QUDA_GCR_INVERTER;
				solver_string_type = "GCR";
			}
			else {
				switch( invParam.solverType ) {
				case CG:
					quda_inv_param.inv_type = QUDA_CG_INVERTER;
					solver_string_type = "CG";
					break;
				case BICGSTAB:
					quda_inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
					solver_string_type = "BICGSTAB";
					break;
				case GCR:
					quda_inv_param.inv_type = QUDA_GCR_INVERTER;
					solver_string_type = "GCR";
					break;
				default:
					quda_inv_param.inv_type = QUDA_GCR_INVERTER;
					solver_string_type = "GCR";
					break;
				}
			}

			// Kappa=0.5 with KAPPA_NORMALIZATION gives A - (1/2)D
			quda_inv_param.kappa = 0.5;
			quda_inv_param.clover_coeff = 1.0; // Dummy, we provide pre-computed clover
			quda_inv_param.Ls = 1;

			quda_inv_param.tol = toDouble(invParam.RsdTarget);
			quda_inv_param.maxiter = invParam.MaxIter;
			quda_inv_param.reliable_delta = toDouble(invParam.Delta);
			quda_inv_param.pipeline = invParam.Pipeline;

			// Use MATPC_SOLUTION: QUDA solves on odd checkerboard
			// (required for QDPJIT which doesn't support full-site spinor fields)
			quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
			quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
			quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;

			quda_inv_param.dagger = QUDA_DAG_NO;
			quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

			quda_inv_param.cpu_prec = cpu_prec;
			quda_inv_param.cuda_prec = gpu_prec;
			quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
			quda_inv_param.cuda_prec_precondition = gpu_half_prec;

			quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
			quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
			quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
			quda_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
			quda_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
#else
			quda_inv_param.dirac_order = QUDA_QDPJIT_DIRAC_ORDER;
			quda_inv_param.input_location = QUDA_CUDA_FIELD_LOCATION;
			quda_inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
#endif

			// Setup padding
			multi1d<int> face_size(4);
			face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
			face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
			face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
			face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;

			int max_face = face_size[0];
			for(int i=1; i <=3; i++) {
				if ( face_size[i] > max_face ) {
					max_face = face_size[i];
				}
			}
			q_gauge_param.ga_pad = max_face;
			quda_inv_param.sp_pad = 0;

			// Clover precision and order
			quda_inv_param.clover_cpu_prec = cpu_prec;
			quda_inv_param.clover_cuda_prec = gpu_prec;
			quda_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;
			quda_inv_param.clover_cuda_prec_precondition = gpu_half_prec;
			quda_inv_param.cl_pad = 0;

			// MG preconditioner precision
			if( invParam.MULTIGRIDParamsP ) {
				const MULTIGRIDSolverParams& ip = *(invParam.MULTIGRIDParams);

				switch( ip.prec ) {
					case HALF:
					quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
					quda_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
					q_gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
					break;
					case SINGLE:
					quda_inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
					quda_inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
					q_gauge_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
					break;
					case DOUBLE:
					quda_inv_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
					quda_inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
					q_gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
					break;
					default:
					quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
					quda_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
					q_gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
					break;
				}

				switch( ip.reconstruct ) {
					case RECONS_NONE:
					q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_NO;
					break;
					case RECONS_8:
					q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_8;
					break;
					case RECONS_12:
					q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;
					break;
					default:
					q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;
					break;
				};
			}

			// ============================================================
			// Load FAT links as gauge field to QUDA
			// ============================================================
			void* gauge[4];
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
			for(int mu=0; mu < Nd; mu++) {
				gauge[mu] = (void *)&(fat_links[mu].elem(all.start()).elem().elem(0,0).real());
			}
#else
			GetMemoryPtrGauge(gauge, fat_links);
#endif
			loadGaugeQuda((void *)gauge, &q_gauge_param);

			// ============================================================
			// Compute and load clover term from THIN links
			// ============================================================
			if( invParam.MULTIGRIDParamsP ) {
				const MULTIGRIDSolverParams& ip = *(invParam.MULTIGRIDParams);
				quda_inv_param.tol_precondition = toDouble(ip.tol[0]);
				quda_inv_param.maxiter_precondition = ip.maxIterations[0];
				quda_inv_param.gcrNkrylov = ip.outer_gcr_nkrylov;
				quda_inv_param.residual_type = static_cast<QudaResidualType>(QUDA_L2_RELATIVE_RESIDUAL);

				switch( ip.schwarzType ) {
					case ADDITIVE_SCHWARZ :
					quda_inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
					break;
					case MULTIPLICATIVE_SCHWARZ :
					quda_inv_param.schwarz_type = QUDA_MULTIPLICATIVE_SCHWARZ;
					break;
					default:
					quda_inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
					break;
				}
				quda_inv_param.precondition_cycle = 1;
			}

			if( invParam.verboseP ) {
				quda_inv_param.verbosity = QUDA_VERBOSE;
			}
			else {
				quda_inv_param.verbosity = QUDA_SUMMARIZE;
			}

			if( invParam.MULTIGRIDParamsP ) {
				quda_inv_param.verbosity_precondition = QUDA_SILENT;
				quda_inv_param.inv_type_precondition = QUDA_MG_INVERTER;
			}
			else {
				quda_inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
				quda_inv_param.verbosity_precondition = QUDA_SILENT;
				quda_inv_param.tol_precondition = 1.0e-1;
				quda_inv_param.maxiter_precondition = 1000;
				quda_inv_param.gcrNkrylov = 10;
				q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_NO;
			}

			QDPIO::cout << solver_string
			            << "Creating CloverTerm from thin links" << std::endl;
			clov->create(thin_fstate, invParam_.CloverParams);
			invclov->create(thin_fstate, invParam_.CloverParams);

			QDPIO::cout << solver_string << "Inverting CloverTerm" << std::endl;
			invclov->choles(0);
			invclov->choles(1);

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
			quda_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;

			multi1d<QUDAPackedClovSite<REALT> > packed_clov;
			packed_clov.resize(all.siteTable().size());
			clov->packForQUDA(packed_clov, 0);
			clov->packForQUDA(packed_clov, 1);

			multi1d<QUDAPackedClovSite<REALT> > packed_invclov(all.siteTable().size());
			invclov->packForQUDA(packed_invclov, 0);
			invclov->packForQUDA(packed_invclov, 1);

			loadCloverQuda(&(packed_clov[0]), &(packed_invclov[0]), &quda_inv_param);
#else
			quda_inv_param.clover_location = QUDA_CUDA_FIELD_LOCATION;
			quda_inv_param.clover_order = QUDA_QDPJIT_CLOVER_ORDER;

			void *clover[2];
			void *cloverInv[2];

			GetMemoryPtrClover(clov->getOffId(),clov->getDiaId(),invclov->getOffId(),invclov->getDiaId());

			loadCloverQuda( (void*)(clover), (void *)(cloverInv), &quda_inv_param);
#endif

			if( invParam.MULTIGRIDParamsP ) {
				quda_inv_param.omega = toDouble((*(invParam.MULTIGRIDParams)).relaxationOmegaOuter);
			}

			// ============================================================
			// MG Subspace: create or recover (only if MG is enabled)
			// ============================================================
			if( invParam.MULTIGRIDParamsP ) {
				if(TheNamedObjMap::Instance().check(invParam.SaveSubspaceID))
				{
					StopWatch update_swatch;
					update_swatch.reset(); update_swatch.start();
					QDPIO::cout << solver_string
					            << "Recovering subspace..." << std::endl;
					subspace_pointers = TheNamedObjMap::Instance().getData< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID);
					MULTIGRIDSolverParams ip = *(invParam.MULTIGRIDParams);
					for(int j=0; j < ip.mg_levels-1; ++j) {
						(subspace_pointers->mg_param).setup_maxiter_refresh[j] = 0;
					}
					updateMultigridQuda(subspace_pointers->preconditioner, &(subspace_pointers->mg_param));
					update_swatch.stop();
					QDPIO::cout << solver_string << "subspace_update_time = "
					            << update_swatch.getTimeInSeconds() << " sec." << std::endl;
				}
				else
				{
					StopWatch create_swatch;
					create_swatch.reset(); create_swatch.start();
					QDPIO::cout << solver_string
					            << "Creating Subspace" << std::endl;
					subspace_pointers = QUDAMGUtils::create_subspace<T>(invParam);

					XMLBufferWriter file_xml;
					push(file_xml, "FileXML");
					pop(file_xml);

					int foo = 5;
					XMLBufferWriter record_xml;
					push(record_xml, "RecordXML");
					write(record_xml, "foo", foo);
					pop(record_xml);

					TheNamedObjMap::Instance().create< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID);
					TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setFileXML(file_xml);
					TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setRecordXML(record_xml);
					TheNamedObjMap::Instance().getData< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID) = subspace_pointers;

					create_swatch.stop();
					QDPIO::cout << solver_string << "subspace_create_time = "
					            << create_swatch.getTimeInSeconds() << " sec." << std::endl;
				}
				quda_inv_param.preconditioner = subspace_pointers->preconditioner;
			}
			else {
				QDPIO::cout << solver_string
				            << "No multigrid: using direct " << solver_string_type
				            << " solver" << std::endl;
				subspace_pointers = nullptr;
			}

			init_swatch.stop();
			QDPIO::cout << solver_string << "init_time = "
			            << init_swatch.getTimeInSeconds() << " sec." << std::endl;
		}

		//! Destructor
		~LinOpSysSolverQUDASLRC()
		{
			QDPIO::cout << solver_string << "Destructing" << std::endl;
			quda_inv_param.preconditioner = nullptr;
			subspace_pointers = nullptr;
			freeGaugeQuda();
			freeCloverQuda();
		}

		//! Return the subset on which the operator acts
		const Subset& subset() const {return A->subset();}

		//! Solve the linear system
		/*!
		 * \param psi      solution ( Modify )
		 * \param chi      source ( Read )
		 * \return syssolver results
		 */
		SystemSolverResults_t operator() (T& psi, const T& chi) const
		{
			SystemSolverResults_t res;

			START_CODE();
			StopWatch swatch;
			swatch.start();

			psi = zero; // Zero initial guess

			// ========================================
			// Prepare Schur complement source
			// ========================================
			// QUDA with MATPC_SOLUTION expects the preconditioned source:
			//   b_hat_o = chi_o - M_oe M_ee^{-1} chi_e
			// where M_ee = A_ee (clover from thin links), M_oe uses fat links.
			// We compute M_oe M_ee^{-1} chi_e by:
			//   1) Apply A_ee^{-1} to even sites of chi
			//   2) Apply full operator A to get M_oe on odd sites
			T mod_chi;
			{
				T tmp1;
				tmp1 = zero;
				invclov->apply(tmp1, chi, PLUS, 0); // tmp1[rb[0]] = A_ee^{-1} chi_e
				// tmp1[rb[1]] = 0, so (*A)(tmp2, tmp1) gives tmp2[rb[1]] = M_oe * tmp1[rb[0]]
				T tmp2;
				(*A)(tmp2, tmp1, PLUS);
				mod_chi[rb[1]] = chi - tmp2; // b_hat_o = chi_o - M_oe M_ee^{-1} chi_e
			}

			if ( invParam.axialGaugeP ) {
				T g_chi, g_psi;

				// Gauge Fix source and initial guess
				g_chi[ rb[1] ] = GFixMat * mod_chi;
				g_psi[ rb[1] ] = GFixMat * psi;
				res = qudaInvert(*clov,
						*invclov,
						g_chi,
						g_psi);
				psi[ rb[1] ] = adj(GFixMat)*g_psi;
			}
			else {
				res = qudaInvert(*clov,
						*invclov,
						mod_chi,
						psi);
			}

			swatch.stop();

			// ========================================
			// Reconstruct even-site solution
			// ========================================
			// QUDA with MATPC_SOLUTION only returns the odd-checkerboard solution.
			// For the unpreconditioned SLRC operator, reconstruct even sites:
			//   psi_e = A_ee^{-1} (chi_e - M_eo * psi_o)
			// where A_ee is clover (thin links), M_eo = -(1/2) D_eo (fat links).
			// Since psi_e = 0, applying the full operator gives tmp_e = M_eo * psi_o.
			{
				T tmp;
				(*A)(tmp, psi, PLUS);
				T even_rhs;
				even_rhs[rb[0]] = chi - tmp;
				invclov->apply(psi, even_rhs, PLUS, 0);
			}

			// ========================================
			// Chroma-side residual verification
			// ========================================
			Double rel_resid;

			if( invParam.SolutionCheckP ) {
				T r;
				r[A->subset()] = chi;
				T tmp;
				(*A)(tmp, psi, PLUS);
				r[A->subset()] -= tmp;
				res.resid = sqrt(norm2(r, A->subset()));

				rel_resid = res.resid / sqrt(norm2(chi, A->subset()));

				QDPIO::cout << solver_string
				            << "QUDA  true residual = " << quda_inv_param.true_res
				            << std::endl;
				QDPIO::cout << solver_string
				            << "Chroma true residual = " << rel_resid
				            << std::endl;
				QDPIO::cout << solver_string << res.n_count
				            << " iterations. Rsd = " << res.resid
				            << " Relative Rsd = " << rel_resid << std::endl;

				// Warn if Chroma residual is much worse than QUDA residual
				if ( toBool( rel_resid > 1.0e-6 ) ) {
					QDPIO::cerr << solver_string
					            << "WARNING: Chroma-side residual " << rel_resid
					            << " is large. Check SLRC gauge/clover field decomposition!"
					            << std::endl;
				}
			}
			else {
				rel_resid = res.resid;
			}

			// Convergence Check/Blow Up
			if ( ! invParam.SilentFailP ) {
				if ( toBool( rel_resid > invParam.RsdToleranceFactor*invParam.RsdTarget) ) {
					QDPIO::cerr << solver_string
					            << "ERROR: Solver residuum is outside tolerance: QUDA resid="
					            << rel_resid << " Desired =" << invParam.RsdTarget
					            << " Max Tolerated = "
					            << invParam.RsdToleranceFactor*invParam.RsdTarget
					            << std::endl;
					QDP_abort(1);
				}
			}

			END_CODE();
			return res;
		}

	private:
		// Hide default constructor
		LinOpSysSolverQUDASLRC() {}

#if 1
		Q links_orig;
#endif

		U GFixMat;
		QudaPrecision_s cpu_prec;
		QudaPrecision_s gpu_prec;
		QudaPrecision_s gpu_half_prec;

		Handle< LinearOperator<T> > A;
		const SysSolverQUDAMULTIGRIDCloverParams invParam;
		QudaGaugeParam q_gauge_param;
		mutable QudaInvertParam quda_inv_param;
		mutable QUDAMGUtils::MGSubspacePointers* subspace_pointers;

		Handle< CloverTermT<T, U> > clov;
		Handle< CloverTermT<T, U> > invclov;

		SystemSolverResults_t qudaInvert(const CloverTermT<T, U>& clover,
				const CloverTermT<T, U>& inv_clov,
				const T& chi_s,
				T& psi_s
		) const;

		std::string solver_string;
		std::string solver_string_type;
	};

} // End namespace

#endif // BUILD_QUDA
#endif
