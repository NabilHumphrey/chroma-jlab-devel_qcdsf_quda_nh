# Integrating QUDA with UNPRECONDITIONED_SLRC in Chroma (QDP-JIT)

## Overview

This guide covers two things:

1. **XML workflow**: A 5-mass multigrid QUDA configuration for your unperturbed propagators, integrated into your existing Feynman-Hellmann pipeline.
2. **C++ code changes**: The adapter layer needed to make QUDA's clover solver work with your `UNPRECONDITIONED_SLRC` fermion action and `SLIC_FERM_STATE`.

The core challenge is that SLRC uses **different gauge fields** for the hopping term (relevant, stout-smeared) and the clover term (irrelevant, unsmeared or differently smeared), while QUDA's clover solver expects to build everything from a single gauge field. The solution is to **pre-compute the clover term on the Chroma side** and pass both it and the appropriate hopping-term gauge field to QUDA separately.

---

## Part 1: Code Changes to Chroma

### 1.1 Understanding the SLRC Decomposition

In your current setup, `SLIC_FERM_STATE` constructs two gauge fields internally:

- **Thin (original) links**: Used for the clover term (the "irrelevant" operator)
- **Stout-smeared links** (`n_smear=1`, `rho=0.1`): Used for the hopping/Dslash term (the "relevant" operator)

The standard Chroma-QUDA clover interface (`syssolver_linop_clover_quda_w.h`) does something like:

```cpp
// Pseudo-code of what the standard QUDA clover interface does:
QudaInvertParam inv_param;
QudaGaugeParam gauge_param;
// ...
loadGaugeQuda(gauge_field_ptr, &gauge_param);   // loads ONE gauge field
loadCloverQuda(NULL, NULL, &inv_param);          // tells QUDA to compute clover from that gauge field
invertQuda(solution_ptr, source_ptr, &inv_param);
```

This won't work for SLRC because QUDA would compute the clover term from the smeared links rather than the thin links (or vice versa).

### 1.2 Strategy: Pre-computed Clover + Separate Gauge Field

QUDA supports receiving a **pre-computed clover field** via `loadCloverQuda()`. The key API call is:

```cpp
// QUDA API signature:
void loadCloverQuda(void *h_clover, void *h_clovInv, QudaInvertParam *inv_param);
```

When you pass non-NULL pointers for `h_clover` and/or `h_clovInv`, QUDA uses your pre-computed clover field instead of computing its own. Combined with loading the smeared gauge field for the hopping term, this gives us exactly the SLRC decomposition.

### 1.3 New Files to Create

You need a new QUDA system solver registered for the SLRC action. The cleanest approach is to create a dedicated syssolver rather than modifying the existing clover one.

#### File: `lib/actions/ferm/invert/quda_solvers/syssolver_linop_slrc_quda_w.h`

This is the main header. The key differences from the standard clover QUDA solver are:

```cpp
#ifndef __syssolver_linop_slrc_quda_w_h__
#define __syssolver_linop_slrc_quda_w_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <quda.h>

namespace Chroma
{
  //! System solver for SLRC + QUDA
  /*!
   *  This solver handles the SLRC decomposition by:
   *  1. Extracting the thin links from the SLIC FermState for the clover term
   *  2. Extracting the fat (stout-smeared) links for the hopping term
   *  3. Computing the clover term on the Chroma side using thin links
   *  4. Passing the fat links as the gauge field to QUDA
   *  5. Passing the pre-computed clover term to QUDA via loadCloverQuda()
   */

  class LinOpSysSolverQUDASLRC : public LinOpSystemSolver<LatticeFermion>
  {
  public:
    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> Q;
    typedef multi1d<LatticeColorMatrix> P;

    //! Constructor
    /*!
     * \param A_        The linear operator (UnprecSLRCLinOp)
     * \param state_    The SLIC ferm state (contains both thin and fat links)
     * \param invParam_ Solver parameters (reuse SysSolverQUDACloverParams)
     */
    LinOpSysSolverQUDASLRC(
      Handle< LinearOperator<T> > A_,
      Handle< FermState<T,P,Q> > state_,
      const SysSolverQUDACloverParams& invParam_);

    //! Destructor
    ~LinOpSysSolverQUDASLRC();

    //! Return the subset on which the operator acts
    const Subset& subset() const { return all; }

    //! Solve
    SystemSolverResults_t operator()(T& psi, const T& chi) const;

  private:
    Handle< LinearOperator<T> > A;
    Handle< FermState<T,P,Q> > state;
    SysSolverQUDACloverParams invParam;

    QudaGaugeParam q_gauge_param;
    QudaInvertParam quda_inv_param;

    // Storage for the gauge and clover fields in QUDA-compatible format
    Handle< CloverTermT<T, P, Q> > clov;      // Chroma clover term
    Handle< CloverTermT<T, P, Q> > invclov;   // Inverse clover term

    // The two gauge fields from the SLIC state
    Q thin_links;   // For clover term computation
    Q fat_links;    // For hopping term (passed to QUDA as gauge field)

    // Multigrid setup handle (for subspace caching)
    void* mg_preconditioner;
    bool mg_setup_done;

    // Private methods
    void extractSLICGaugeFields(Handle< FermState<T,P,Q> > state_);
    void computeCloverTerm();
    void setupQUDAParams();
    void loadFieldsToQUDA();
  };

} // namespace Chroma

#endif // BUILD_QUDA
#endif
```

#### File: `lib/actions/ferm/invert/quda_solvers/syssolver_linop_slrc_quda_w.cc`

The implementation. The critical method is the constructor and solver:

```cpp
#include "chroma_config.h"

#ifdef BUILD_QUDA

#include "actions/ferm/invert/quda_solvers/syssolver_linop_slrc_quda_w.h"
#include "actions/ferm/fermstates/slic_fermstate_w.h"

namespace Chroma
{
  // ------------------------------------------------------------------
  // Extract thin and fat links from the SLIC FermState
  // ------------------------------------------------------------------
  void LinOpSysSolverQUDASLRC::extractSLICGaugeFields(
    Handle< FermState<T,P,Q> > state_)
  {
    // The SLIC FermState stores both the original (thin) and
    // stout-smeared (fat) links. We need to access both.
    //
    // The key insight: state_->getLinks() returns the FAT links
    // (these are what the hopping term sees). For the thin links,
    // we need access to the underlying unsmeared gauge field.
    //
    // How to get the thin links depends on your SLIC FermState
    // implementation. Common approaches:
    //
    // Option A: If your SLICFermState exposes getThinLinks():
    //   thin_links = dynamic_cast<SLICFermState<T,P,Q>&>(*state_).getThinLinks();
    //
    // Option B: If it only exposes getLinks() (fat) and you need to
    //   reconstruct thin links, you need to modify SLICFermState to
    //   store and expose them. See Section 1.4 below.

    // Get fat links (for hopping term / QUDA gauge field)
    fat_links = state_->getLinks();

    // Get thin links (for clover term computation)
    // This requires the SLICFermState to expose them - see Section 1.4
    try {
      const SLICFermState<T,P,Q>& slic_state =
        dynamic_cast<const SLICFermState<T,P,Q>&>(*state_);
      thin_links = slic_state.getThinLinks();
    }
    catch (std::bad_cast& e) {
      QDPIO::cerr << "SLRC QUDA solver: FermState is not a SLICFermState!"
                  << std::endl;
      QDP_abort(1);
    }
  }

  // ------------------------------------------------------------------
  // Compute the clover term using thin links
  // ------------------------------------------------------------------
  void LinOpSysSolverQUDASLRC::computeCloverTerm()
  {
    // Build a temporary simple FermState from the thin links
    // so the clover term constructor sees unsmeared links
    Handle< FermBC<T,P,Q> > fbc = state->getFermBC();
    Handle< FermState<T,P,Q> > thin_state(
      new PeriodicFermState<T,P,Q>(thin_links));

    // Apply boundary conditions to the thin state
    // (the clover term needs BC-aware links)
    fbc->modify(thin_links);

    // Create the clover term from thin links
    // clovCoeff comes from your action parameters
    CloverFermActParams clparam;
    clparam.clovCoeffR = invParam.CloverParams.clovCoeffR;
    clparam.clovCoeffT = invParam.CloverParams.clovCoeffT;

    clov = new CloverTermT<T,P,Q>();
    clov->create(thin_state, clparam);

    // Also compute the inverse clover term (QUDA can use this)
    invclov = new CloverTermT<T,P,Q>();
    invclov->create(thin_state, clparam, clov->getTriBuffer());
    // NOTE: The exact API for getting/passing the triangular buffer
    // depends on your Chroma version. You may need to just create
    // invclov separately and call choles() on it.
  }

  // ------------------------------------------------------------------
  // Constructor
  // ------------------------------------------------------------------
  LinOpSysSolverQUDASLRC::LinOpSysSolverQUDASLRC(
    Handle< LinearOperator<T> > A_,
    Handle< FermState<T,P,Q> > state_,
    const SysSolverQUDACloverParams& invParam_)
    : A(A_), state(state_), invParam(invParam_),
      mg_preconditioner(nullptr), mg_setup_done(false)
  {
    QDPIO::cout << "LinOpSysSolverQUDASLRC: Initializing SLRC+QUDA solver"
                << std::endl;

    // Step 1: Extract the two gauge fields from SLIC state
    extractSLICGaugeFields(state_);

    // Step 2: Compute clover term from thin links
    computeCloverTerm();

    // Step 3: Set up QUDA parameter structures
    setupQUDAParams();

    // Step 4: Load gauge field (fat links) and clover field into QUDA
    loadFieldsToQUDA();

    QDPIO::cout << "LinOpSysSolverQUDASLRC: Initialization complete"
                << std::endl;
  }

  // ------------------------------------------------------------------
  // Load fields into QUDA
  // ------------------------------------------------------------------
  void LinOpSysSolverQUDASLRC::loadFieldsToQUDA()
  {
    // --- Load the FAT links as the gauge field ---
    // QUDA will use these for the hopping term (Dslash)
    //
    // Convert fat_links to the format QUDA expects.
    // This is the same conversion the standard clover solver does,
    // but we're explicitly passing the fat (stout-smeared) links.

    // The gauge field packing code - reuse from existing QUDA interface
    // (typically in quda_interface_utils.h or similar in your branch)
    void* gauge[4];
    for(int mu=0; mu < Nd; mu++) {
      gauge[mu] = (void*)&(fat_links[mu].elem(all.start()).elem().elem(0,0).real());
    }
    loadGaugeQuda(gauge, &q_gauge_param);

    // --- Load the pre-computed clover term ---
    // This is the critical SLRC-specific step.
    //
    // QUDA expects the clover field packed in a specific format.
    // When we pass non-NULL h_clover, QUDA uses our field instead
    // of computing its own from the loaded gauge field.
    //
    // Set the flag to tell QUDA we're providing our own clover field:
    quda_inv_param.clover_coeff = 0.0;  // We already included it
    quda_inv_param.compute_clover = QUDA_BOOLEAN_FALSE;

    // Pack the Chroma clover term into QUDA's expected layout
    // QUDA wants the clover field as a flat array in a specific order.
    // The packing depends on your QUDA version - see packCloverField() below.
    multi1d<QUDAPackedClovSite<REAL>> packed_clover(Layout::sitesOnNode());
    multi1d<QUDAPackedClovSite<REAL>> packed_clover_inv(Layout::sitesOnNode());

    clov->packForQUDA(packed_clover, 0);   // cb=0
    clov->packForQUDA(packed_clover, 1);   // cb=1
    invclov->packForQUDA(packed_clover_inv, 0);
    invclov->packForQUDA(packed_clover_inv, 1);

    // NOTE: The packForQUDA method may not exist in your Chroma branch.
    // In that case, you need to write a packing routine - see Section 1.6.

    loadCloverQuda(
      (void*)packed_clover.slice(),
      (void*)packed_clover_inv.slice(),
      &quda_inv_param
    );
  }

  // ------------------------------------------------------------------
  // Solve
  // ------------------------------------------------------------------
  SystemSolverResults_t LinOpSysSolverQUDASLRC::operator()(
    T& psi, const T& chi) const
  {
    SystemSolverResults_t res;

    // Pack source
    // (same as standard QUDA clover solver)
    void* spinorIn = (void*)&(chi.elem(all.start()).elem(0).elem(0).real());
    void* spinorOut = (void*)&(psi.elem(all.start()).elem(0).elem(0).real());

    // Solve
    invertQuda(spinorOut, spinorIn, &quda_inv_param);

    // Extract residual info
    res.n_count = quda_inv_param.iter;
    res.resid = quda_inv_param.true_res;

    QDPIO::cout << "QUDA_SLRC_SOLVER: iters = " << res.n_count
                << "  residual = " << res.resid << std::endl;

    return res;
  }

} // namespace Chroma

#endif // BUILD_QUDA
```

### 1.4 Modifying SLICFermState to Expose Thin Links

Your `SLICFermState` (typically in `lib/actions/ferm/fermstates/slic_fermstate_w.h`) stores the smeared links but may not expose the original thin links. You need to add a getter:

```cpp
// In slic_fermstate_w.h, inside the SLICFermState class:

// Add a member to store the thin (unsmeared) links if not already present:
private:
  Q thin_links_;   // Original unsmeared gauge links

// In the constructor, BEFORE smearing, save a copy:
//   thin_links_ = u;  // where u is the original gauge field

public:
  //! Get the original unsmeared (thin) links
  const Q& getThinLinks() const { return thin_links_; }
```

**Check your existing code first** — some versions of `SLICFermState` already store the original links internally (often called `original_links` or `u_thin`). If so, you just need to add a public accessor.

### 1.5 Registering the New Solver

You need to register the new solver type so Chroma's XML parser recognizes it. In the solver registration file (typically `lib/actions/ferm/invert/quda_solvers/quda_solver_registration.cc` or wherever your existing QUDA solvers are registered):

```cpp
#include "actions/ferm/invert/quda_solvers/syssolver_linop_slrc_quda_w.h"

// In the registration function:
namespace {
  //! Callback function
  LinOpSystemSolver<LatticeFermion>*
  createSLRCQUDA(XMLReader& xml_in,
                  const std::string& path,
                  Handle< FermState<LatticeFermion,
                    multi1d<LatticeColorMatrix>,
                    multi1d<LatticeColorMatrix>> > state,
                  Handle< LinearOperator<LatticeFermion> > A)
  {
    return new LinOpSysSolverQUDASLRC(A, state,
      SysSolverQUDACloverParams(xml_in, path));
  }

  //! Register
  const std::string name("QUDA_SLRC_CLOVER_INVERTER");
  bool registered = TheLinOpFermSystemSolverFactory::Instance()
    .registerObject(name, createSLRCQUDA);
}
```

You also need to make sure the `UNPRECONDITIONED_SLRC` FermAct knows that it can use this solver. In your SLRC fermion action class (likely `unprec_slrc_fermact_w.h/.cc`), ensure the `invLinOp` method can dispatch to the QUDA solver. This typically means ensuring the FermAct's `qprop` or solver-creation machinery respects the `invType` string from the XML and checks the system solver factory.

### 1.6 Clover Field Packing for QUDA

If your Chroma branch doesn't already have a `packForQUDA` method on the clover term, you need a packing routine. QUDA expects the clover field as 72 real numbers per site (two 6×6 Hermitian matrices for the two chiral blocks, each stored as 36 reals):

```cpp
// Utility: Pack Chroma clover term into QUDA layout
// This goes in a utility header, e.g., quda_slrc_utils.h

template<typename T, typename P, typename Q>
void packChromaCloverForQUDA(
  const CloverTermT<T,P,Q>& clover,
  double* quda_clover,      // output: 72 doubles per site
  const Subset& sub)
{
  // QUDA clover layout per site:
  //   First 36 reals: upper-left chiral block (6x6 hermitian)
  //   Next 36 reals:  lower-right chiral block (6x6 hermitian)
  //
  // Hermitian 6x6 stored as: diagonal (6 reals) then
  // off-diagonal upper triangle (15 complex = 30 reals)
  //
  // The exact packing order must match QUDA's expectations.
  // See QUDA's clover_field.cpp for the reference layout.
  //
  // In practice, look at how your existing Chroma QUDA interface
  // packs the clover field in syssolver_linop_clover_quda_w.cc.
  // The code is typically in a function called something like
  // qudaPackClover() or similar. You can reuse that code directly,
  // just feeding it the clover term you computed from thin links
  // rather than letting QUDA compute its own.
}
```

> **Practical shortcut**: The existing `syssolver_linop_clover_quda_w.cc` in your branch almost certainly already has clover packing code (it may compute the clover in Chroma for verification or for certain preconditioners). Find that code and reuse it. The only change is the input: your clover term computed from thin links.

---

## Part 2: XML Workflow — 5-Mass Multigrid

### 2.1 Concept

For 5 κ values (call them κ₁ through κ₅, lightest to heaviest):

1. **Build the MG subspace** on the first (lightest) mass. This is the expensive setup.
2. **Solve** the lightest mass propagator using MG-preconditioned GCR.
3. **Reuse the MG subspace** for the remaining 4 masses — each solve is just the cheap MG-preconditioned GCR, no setup cost.
4. After all 5 unperturbed propagators are done, proceed with FH inversions using the existing BiCGStab solver (with `initial_guess_id` pointing at the relevant unperturbed prop).

The MG subspace is stored as a named object in Chroma's object map, so it persists across propagator tasks within the same Chroma run.

### 2.2 Full XML Template

Below is the complete XML for a 5-mass workflow. For clarity, only one FH block (one operator, one lambda) is shown per mass — your XML generator would replicate these across all operator/lambda combinations as in your current setup.

```xml
<?xml version="1.0"?>
<!--
  Feynman-Hellmann Spectrum with QUDA Multigrid
  5 masses with shared MG subspace for unperturbed propagators
  FH propagators use QDP-JIT BiCGStab with initial guess
  Lattice: 32^3 x 64
-->

<chroma>
    <Param>
        <InlineMeasurements>

            <!-- ============================================================ -->
            <!-- SOURCE (shared across all masses)                            -->
            <!-- ============================================================ -->

            <elem>
                <n>MAKE_SOURCE</n>
                <Frequency>1</Frequency>
                <Param>
                    <version>6</version>
                    <Source>
                        <version>3</version>
                        <SourceType>SHELL_SOURCE</SourceType>
                        <j_decay>3</j_decay>
                        <t_srce>0 0 0 48</t_srce>
                        <quark_smear_lastP>false</quark_smear_lastP>
                        <SmearingParam>
                            <wvf_kind>GAUGE_INV_JACOBI</wvf_kind>
                            <wvf_param>0.21</wvf_param>
                            <wvfIntPar>75</wvfIntPar>
                            <no_smear_dir>3</no_smear_dir>
                        </SmearingParam>
                        <Displacement>
                            <version>1</version>
                            <DisplacementType>NONE</DisplacementType>
                        </Displacement>
                        <LinkSmearing>
                            <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
                            <link_smear_fact>0.1</link_smear_fact>
                            <link_smear_num>2</link_smear_num>
                            <no_smear_dir>3</no_smear_dir>
                        </LinkSmearing>
                    </Source>
                </Param>
                <NamedObject>
                    <gauge_id>default_gauge_field</gauge_id>
                    <source_id>sh_source_1</source_id>
                </NamedObject>
            </elem>

            <!-- ============================================================ -->
            <!-- MASS 1 (lightest): BUILD MG subspace + solve                 -->
            <!-- Kappa: 0.120900 (your lightest mass)                         -->
            <!-- ============================================================ -->

            <elem>
                <annotation>Unperturbed prop, mass 1 (lightest) - builds MG subspace</annotation>
                <n>PROPAGATOR</n>
                <Frequency>1</Frequency>
                <Param>
                    <version>10</version>
                    <quarkSpinType>FULL</quarkSpinType>
                    <obsvP>false</obsvP>
                    <numRetries>1</numRetries>
                    <FermionAction>
                        <FermAct>UNPRECONDITIONED_SLRC</FermAct>
                        <Kappa>0.120900</Kappa>
                        <clovCoeff>2.65</clovCoeff>
                        <FermState>
                            <n>SLIC_FERM_STATE</n>
                            <n_smear>1</n_smear>
                            <rho>0.1</rho>
                            <orthog_dir>5</orthog_dir>
                            <FermionBC>
                                <FermBC>SIMPLE_FERMBC</FermBC>
                                <boundary>1 1 1 -1</boundary>
                            </FermionBC>
                        </FermState>
                    </FermionAction>
                    <InvertParam>
                        <invType>QUDA_SLRC_CLOVER_INVERTER</invType>
                        <SolverType>GCR</SolverType>
                        <RsdTarget>1.0e-12</RsdTarget>
                        <MaxIter>10000</MaxIter>
                        <AntiPeriodicT>true</AntiPeriodicT>
                        <Delta>1.0e-1</Delta>
                        <Pipeline>8</Pipeline>
                        <Verbosity>VERBOSE</Verbosity>

                        <!-- Precision control (mixed precision is key) -->
                        <CudaPrecision>DOUBLE</CudaPrecision>
                        <CudaSloppyPrecision>SINGLE</CudaSloppyPrecision>
                        <CudaRefinementPrecision>SINGLE</CudaRefinementPrecision>

                        <!-- ======================================== -->
                        <!-- MULTIGRID CONFIGURATION                  -->
                        <!-- ======================================== -->
                        <MULTIGRIDParams>
                            <!-- Create and store the subspace -->
                            <SubspaceCreate>true</SubspaceCreate>
                            <SubspaceId>mg_subspace_slrc</SubspaceId>

                            <!-- 2 or 3 levels; start with 2 for a 32^3x64 -->
                            <Levels>2</Levels>

                            <!-- Level 0 -> Level 1 blocking -->
                            <!-- 32^3x64 -> 8^3x16 with 4^4 blocks -->
                            <Blocking>
                                <elem>4 4 4 4</elem>
                            </Blocking>

                            <!-- Null vectors per level -->
                            <NullVectors>
                                <elem>24</elem>
                            </NullVectors>

                            <!-- Setup solver on each level -->
                            <SetupSolveType>
                                <elem>BICGSTAB</elem>
                            </SetupSolveType>
                            <SetupTol>
                                <elem>5.0e-6</elem>
                            </SetupTol>
                            <SetupMaxIter>
                                <elem>500</elem>
                            </SetupMaxIter>
                            <NumSetupIter>1</NumSetupIter>

                            <!-- Smoother on each level -->
                            <SmootherType>
                                <elem>CA_GCR</elem>
                                <elem>CA_GCR</elem>
                            </SmootherType>
                            <SmootherTol>
                                <elem>0.25</elem>
                                <elem>0.25</elem>
                            </SmootherTol>
                            <SmootherSchwarzCycle>
                                <elem>1</elem>
                                <elem>1</elem>
                            </SmootherSchwarzCycle>

                            <!-- Coarse solver -->
                            <CoarseSolverType>
                                <elem>GCR</elem>
                            </CoarseSolverType>
                            <CoarseSolverTol>
                                <elem>0.1</elem>
                            </CoarseSolverTol>
                            <CoarseSolverMaxIter>
                                <elem>12</elem>
                            </CoarseSolverMaxIter>

                            <!-- Pre/post smoothing -->
                            <PreSmooth>
                                <elem>0</elem>
                                <elem>0</elem>
                            </PreSmooth>
                            <PostSmooth>
                                <elem>8</elem>
                                <elem>8</elem>
                            </PostSmooth>

                            <!-- MG cycle type: 1 = V-cycle -->
                            <CycleType>1</CycleType>

                            <!-- Sloppy precision on coarse grids -->
                            <CoarsePrecision>SINGLE</CoarsePrecision>

                            <!-- Subspace generation precision -->
                            <SubspacePrecision>SINGLE</SubspacePrecision>
                        </MULTIGRIDParams>
                    </InvertParam>
                </Param>
                <NamedObject>
                    <gauge_id>default_gauge_field</gauge_id>
                    <source_id>sh_source_1</source_id>
                    <prop_id>prop_m1</prop_id>
                </NamedObject>
            </elem>

            <!-- Sink smear mass 1 -->
            <elem>
                <n>SINK_SMEAR</n>
                <Frequency>1</Frequency>
                <Param>
                    <version>5</version>
                    <Sink>
                        <version>2</version>
                        <SinkType>SHELL_SINK</SinkType>
                        <j_decay>3</j_decay>
                        <SmearingParam>
                            <wvf_kind>GAUGE_INV_JACOBI</wvf_kind>
                            <wvf_param>0.21</wvf_param>
                            <wvfIntPar>75</wvfIntPar>
                            <no_smear_dir>3</no_smear_dir>
                        </SmearingParam>
                        <Displacement>
                            <version>1</version>
                            <DisplacementType>NONE</DisplacementType>
                        </Displacement>
                        <LinkSmearing>
                            <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
                            <link_smear_fact>0.1</link_smear_fact>
                            <link_smear_num>2</link_smear_num>
                            <no_smear_dir>3</no_smear_dir>
                        </LinkSmearing>
                    </Sink>
                </Param>
                <NamedObject>
                    <gauge_id>default_gauge_field</gauge_id>
                    <prop_id>prop_m1</prop_id>
                    <smeared_prop_id>prop_m1_smeared</smeared_prop_id>
                </NamedObject>
            </elem>

            <!-- FH propagators for mass 1 (BiCGStab, QDP-JIT) -->
            <!-- Repeat your existing FH blocks here, pointing at prop_m1 -->
            <elem>
                <annotation>FH propagator (mass 1, mu=1, lambda=0.025)</annotation>
                <n>PROPAGATOR</n>
                <Frequency>1</Frequency>
                <Param>
                    <version>10</version>
                    <quarkSpinType>FULL</quarkSpinType>
                    <obsvP>false</obsvP>
                    <numRetries>1</numRetries>
                    <FermionAction>
                        <FermAct>UNPREC_SLRC_FEYNHELL</FermAct>
                        <Kappa>0.120900</Kappa>
                        <FeynHellParam>
                            <elem>
                                <LambdaReal>0.025</LambdaReal>
                                <LambdaImag>0.0</LambdaImag>
                                <Operator>1</Operator>
                                <Momentum>4 1 0</Momentum>
                                <Source>0 0 0 48</Source>
                                <NoiseReal>1.0</NoiseReal>
                                <NoiseImag>0.0</NoiseImag>
                            </elem>
                        </FeynHellParam>
                        <clovCoeff>2.65</clovCoeff>
                        <FermState>
                            <n>SLIC_FERM_STATE</n>
                            <n_smear>1</n_smear>
                            <rho>0.1</rho>
                            <orthog_dir>5</orthog_dir>
                            <FermionBC>
                                <FermBC>SIMPLE_FERMBC</FermBC>
                                <boundary>1 1 1 -1</boundary>
                            </FermionBC>
                        </FermState>
                    </FermionAction>
                    <InvertParam>
                        <!-- FH props stay on QDP-JIT BiCGStab -->
                        <invType>BICGSTAB_INVERTER</invType>
                        <RsdBiCGStab>1.0e-12</RsdBiCGStab>
                        <MaxBiCGStab>20000</MaxBiCGStab>
                        <MROver>1.1</MROver>
                    </InvertParam>
                </Param>
                <NamedObject>
                    <gauge_id>default_gauge_field</gauge_id>
                    <source_id>sh_source_1</source_id>
                    <prop_id>fh_prop_m1_mu1_lam1</prop_id>
                    <initial_guess_id>prop_m1</initial_guess_id>
                </NamedObject>
            </elem>

            <!-- ... more FH blocks for mass 1 (all mu, lambda combos) ... -->
            <!-- ... sink smear, contractions, cleanup for mass 1 FH props ... -->

            <!-- ============================================================ -->
            <!-- MASS 2: REUSE MG subspace                                    -->
            <!-- Kappa: 0.120800 (example)                                    -->
            <!-- ============================================================ -->

            <elem>
                <annotation>Unperturbed prop, mass 2 - reuses MG subspace</annotation>
                <n>PROPAGATOR</n>
                <Frequency>1</Frequency>
                <Param>
                    <version>10</version>
                    <quarkSpinType>FULL</quarkSpinType>
                    <obsvP>false</obsvP>
                    <numRetries>1</numRetries>
                    <FermionAction>
                        <FermAct>UNPRECONDITIONED_SLRC</FermAct>
                        <Kappa>0.120800</Kappa>
                        <clovCoeff>2.65</clovCoeff>
                        <FermState>
                            <n>SLIC_FERM_STATE</n>
                            <n_smear>1</n_smear>
                            <rho>0.1</rho>
                            <orthog_dir>5</orthog_dir>
                            <FermionBC>
                                <FermBC>SIMPLE_FERMBC</FermBC>
                                <boundary>1 1 1 -1</boundary>
                            </FermionBC>
                        </FermState>
                    </FermionAction>
                    <InvertParam>
                        <invType>QUDA_SLRC_CLOVER_INVERTER</invType>
                        <SolverType>GCR</SolverType>
                        <RsdTarget>1.0e-12</RsdTarget>
                        <MaxIter>10000</MaxIter>
                        <AntiPeriodicT>true</AntiPeriodicT>
                        <Delta>1.0e-1</Delta>
                        <Pipeline>8</Pipeline>
                        <Verbosity>SUMMARIZE</Verbosity>

                        <CudaPrecision>DOUBLE</CudaPrecision>
                        <CudaSloppyPrecision>SINGLE</CudaSloppyPrecision>
                        <CudaRefinementPrecision>SINGLE</CudaRefinementPrecision>

                        <MULTIGRIDParams>
                            <!-- KEY DIFFERENCE: load existing subspace -->
                            <SubspaceCreate>false</SubspaceCreate>
                            <SubspaceId>mg_subspace_slrc</SubspaceId>

                            <!-- All other MG params same as mass 1 -->
                            <Levels>2</Levels>
                            <Blocking>
                                <elem>4 4 4 4</elem>
                            </Blocking>
                            <NullVectors>
                                <elem>24</elem>
                            </NullVectors>

                            <SmootherType>
                                <elem>CA_GCR</elem>
                                <elem>CA_GCR</elem>
                            </SmootherType>
                            <SmootherTol>
                                <elem>0.25</elem>
                                <elem>0.25</elem>
                            </SmootherTol>
                            <SmootherSchwarzCycle>
                                <elem>1</elem>
                                <elem>1</elem>
                            </SmootherSchwarzCycle>

                            <CoarseSolverType>
                                <elem>GCR</elem>
                            </CoarseSolverType>
                            <CoarseSolverTol>
                                <elem>0.1</elem>
                            </CoarseSolverTol>
                            <CoarseSolverMaxIter>
                                <elem>12</elem>
                            </CoarseSolverMaxIter>

                            <PreSmooth>
                                <elem>0</elem>
                                <elem>0</elem>
                            </PreSmooth>
                            <PostSmooth>
                                <elem>8</elem>
                                <elem>8</elem>
                            </PostSmooth>

                            <CycleType>1</CycleType>
                            <CoarsePrecision>SINGLE</CoarsePrecision>
                            <SubspacePrecision>SINGLE</SubspacePrecision>
                        </MULTIGRIDParams>
                    </InvertParam>
                </Param>
                <NamedObject>
                    <gauge_id>default_gauge_field</gauge_id>
                    <source_id>sh_source_1</source_id>
                    <prop_id>prop_m2</prop_id>
                </NamedObject>
            </elem>

            <!-- Sink smear, FH props, contractions for mass 2 ... -->
            <!-- (same pattern as mass 1, with Kappa=0.120800 and prop_m2) -->

            <!-- ============================================================ -->
            <!-- MASS 3: REUSE MG subspace (Kappa: 0.120700)                  -->
            <!-- ============================================================ -->

            <!-- Same pattern as mass 2: SubspaceCreate=false,
                 SubspaceId=mg_subspace_slrc, Kappa=0.120700 -->

            <!-- ============================================================ -->
            <!-- MASS 4: REUSE MG subspace (Kappa: 0.120600)                  -->
            <!-- ============================================================ -->

            <!-- Same pattern -->

            <!-- ============================================================ -->
            <!-- MASS 5: REUSE MG subspace (Kappa: 0.120500)                  -->
            <!-- ============================================================ -->

            <!-- Same pattern. This is the heaviest mass and will converge
                 fastest with the MG preconditioner. -->

            <!-- ============================================================ -->
            <!-- CLEANUP                                                      -->
            <!-- ============================================================ -->

            <!-- Delete MG subspace when all masses are done -->
            <elem>
                <annotation>Delete MG subspace</annotation>
                <n>ERASE_NAMED_OBJECT</n>
                <Frequency>1</Frequency>
                <NamedObject>
                    <object_id>mg_subspace_slrc</object_id>
                    <object_type>MGSubspace</object_type>
                </NamedObject>
            </elem>

            <!-- Delete source and remaining props -->
            <elem>
                <annotation>Delete source</annotation>
                <n>ERASE_NAMED_OBJECT</n>
                <Frequency>1</Frequency>
                <NamedObject>
                    <object_id>sh_source_1</object_id>
                    <object_type>LatticePropagator</object_type>
                </NamedObject>
            </elem>

            <!-- ... delete all remaining props ... -->

        </InlineMeasurements>
        <nrow>32 32 32 64</nrow>
    </Param>

    <RNG>
        <Seed>
            <elem>11</elem>
            <elem>11</elem>
            <elem>11</elem>
            <elem>0</elem>
        </Seed>
    </RNG>

    <Cfg>
        <cfg_type>SCIDAC</cfg_type>
        <cfg_file>/path/to/your/config.lime</cfg_file>
    </Cfg>

</chroma>
```

### 2.3 Memory Management Notes

For 5 masses on a 32³×64 lattice, memory management is critical. Each propagator is 32³ × 64 × 12 × 2 × 8 bytes ≈ 1.2 GB (double precision, 12 spin-colour, complex). You're going to have several of these in memory simultaneously.

**Recommended workflow order per mass** (minimizes peak memory):

```
For each mass:
  1. Solve unperturbed prop (QUDA MG)        → prop_mN in memory
  2. Sink smear                               → prop_mN_smeared in memory
  3. For each (mu, lambda):
     a. Solve FH prop (BiCGStab + initial guess from prop_mN)
     b. Sink smear FH prop
     c. Contract (meson, baryon, multi-hadron)
     d. ERASE_NAMED_OBJECT the smeared FH prop
     e. ERASE_NAMED_OBJECT the unsmeared FH prop
  4. ERASE_NAMED_OBJECT prop_mN (unsmeared unperturbed)
     (keep prop_mN_smeared only if needed for contractions with next mass)
  5. ERASE_NAMED_OBJECT prop_mN_smeared when done
```

This is essentially what your current XML already does for the single-mass case. The key addition is that the MG subspace (`mg_subspace_slrc`) persists across all 5 masses and is only deleted at the very end.

**GPU memory budget** (rough estimates for a 32³×64 lattice):

| Item | GPU Memory |
|------|-----------|
| QUDA gauge field | ~150 MB |
| QUDA clover field | ~300 MB |
| QUDA solver workspace | ~500 MB–1 GB |
| MG subspace (24 null vectors) | ~1–2 GB |
| QDP-JIT pool (gauge + props) | ~2–4 GB |
| **Total** | **~4–8 GB** |

An A100 (40 or 80 GB) will be comfortable. A V100 (32 GB) should work but will be tighter. Set `QUDA_RESOURCE_PATH` to cache autotuning results.

---

## Part 3: Environment Variables and Runtime

### 3.1 Essential Environment Variables

```bash
# QUDA tuning cache (avoids re-autotuning every run)
export QUDA_RESOURCE_PATH=/path/to/quda_tunecache

# QDP-JIT GPU memory pool limit (leave room for QUDA)
# Adjust based on your GPU memory — start conservative
export QDP_JIT_GPU_POOL_SIZE=4000   # in MB

# QUDA verbosity (useful during debugging, reduce later)
export QUDA_ENABLE_TUNING=1
export QUDA_TUNE_VERSION_CHECK=0

# Ensure QUDA and QDP-JIT use the same GPU
export CUDA_VISIBLE_DEVICES=0       # or appropriate GPU ID
```

### 3.2 MPI Launch

Ensure the grid decomposition is consistent:

```bash
# Example: 4 MPI ranks, 2x2x1x1 decomposition
mpirun -np 4 chroma -i input.xml -geom 2 2 1 1
```

QUDA picks up the grid topology from Chroma's initialization. If you get topology mismatch errors, check that the `-geom` values multiply to your MPI rank count.

---

## Part 4: Tuning the Multigrid Parameters

The MG parameters in the XML above are reasonable starting points for a 32³×64 lattice. Here is guidance for tuning:

### 4.1 Blocking

For a 32³×64 lattice with 2 levels: `4 4 4 4` gives a coarse grid of 8³×16 = 8192 sites. This is a good size. If you go to 3 levels, use `4 4 4 4` then `2 2 2 2` for a coarsest grid of 4³×8 = 512 sites.

### 4.2 Null Vectors

Start with 24 null vectors. If convergence is poor (outer GCR takes >50 iterations), increase to 32. Beyond 32 rarely helps and costs setup time and memory.

### 4.3 Smoother

`CA_GCR` with 8 post-smoothing steps and tolerance 0.25 is a good starting point. If the smoother dominates runtime, try reducing to 4 steps or loosening to 0.5.

### 4.4 When to Add a Third Level

For heavier-than-physical pion masses on 32³×64, 2 levels is usually sufficient. You'd want 3 levels if you go to larger volumes (48³ or 64³) or significantly lighter quarks. The cost of adding a third level is more null vectors to compute and more memory, but the payoff is fewer outer iterations.

### 4.5 Validation

Run with `<Verbosity>VERBOSE</Verbosity>` initially. You want to see:

- MG setup completing in O(10–30) seconds
- Outer GCR converging in O(10–50) iterations for the lightest mass
- Heavier masses converging in fewer iterations
- True residual matching the requested `1.0e-12`

If the outer solver stalls or diverges, the usual fixes are: more null vectors, tighter smoother tolerance, or switching to 3 levels.

---

## Part 5: Validating the QUDA Propagators

Introducing a new solver backend and a non-trivial adapter layer (the SLRC decomposition) means you need to carefully verify that QUDA is inverting the correct operator. There are several levels of validation, ordered from cheapest to most thorough.

### 5.1 Level 0: Solver Residual Check (Automatic)

QUDA reports two residuals: the **iterated residual** (accumulated during the solve) and the **true residual** (recomputed from scratch at the end: ||b − Dx||/||b||). These are printed in the Chroma log when `<Verbosity>VERBOSE</Verbosity>` is set.

What to look for:

- The true residual should be at or below your `<RsdTarget>1.0e-12</RsdTarget>`.
- The iterated and true residuals should agree to within an order of magnitude. If the true residual is much larger than the iterated residual, this indicates a problem — either the operator QUDA is inverting doesn't match the one used for the residual check, or there's a numerical instability.

**However**, this check has a critical blind spot for your case: QUDA computes the true residual using its own internal operator, which is built from whatever gauge and clover fields you loaded. If your SLRC adapter has a bug (e.g., it accidentally passes the fat links for both the clover and hopping terms), QUDA will happily report a tiny residual — it's self-consistent, just wrong. So this check is necessary but not sufficient.

### 5.2 Level 1: Chroma-Side Residual Check (Essential)

This is the most important test. After QUDA returns the solution ψ, apply the **Chroma** SLRC Dirac operator to ψ and check the residual on the Chroma side. This verifies that QUDA's solution actually satisfies the correct SLRC equation, not just QUDA's internal operator.

Add this to your `LinOpSysSolverQUDASLRC::operator()` method:

```cpp
SystemSolverResults_t LinOpSysSolverQUDASLRC::operator()(
  T& psi, const T& chi) const
{
  SystemSolverResults_t res;

  // ... (QUDA solve as before) ...

  // ========================================
  // Chroma-side residual verification
  // ========================================
  {
    T r;
    (*A)(r, psi, PLUS);          // r = D_SLRC * psi (using Chroma's operator)
    r -= chi;                     // r = D_SLRC * psi - b
    Double r_norm = sqrt(norm2(r));
    Double b_norm = sqrt(norm2(chi));
    Double chroma_resid = r_norm / b_norm;

    QDPIO::cout << "QUDA_SLRC_SOLVER: QUDA  true residual = "
                << res.resid << std::endl;
    QDPIO::cout << "QUDA_SLRC_SOLVER: Chroma true residual = "
                << chroma_resid << std::endl;

    // These should agree to within mixed-precision rounding
    // (typically within a factor of 2-5 of each other)
    if (toDouble(chroma_resid) > 1.0e-10) {
      QDPIO::cerr << "WARNING: Chroma-side residual " << chroma_resid
                  << " significantly exceeds target. "
                  << "Check SLRC gauge/clover field decomposition!"
                  << std::endl;
    }
  }

  return res;
}
```

**Interpreting the results:**

| QUDA residual | Chroma residual | Diagnosis |
|:--|:--|:--|
| ~1e-12 | ~1e-12 | Correct. Both operators agree. |
| ~1e-12 | ~1e-1 to ~1e-5 | Bug in the SLRC adapter. QUDA is inverting the wrong operator (e.g., wrong gauge field for clover or hopping term). |
| ~1e-12 | ~1e-10 to ~1e-11 | Acceptable. Small discrepancy from mixed precision. |
| ~1e-8 | ~1e-8 | Solver didn't converge. Check MG parameters. |

The second row is the failure mode specific to your SLRC adapter. If you see it, the most common causes are: thin and fat links being swapped, boundary conditions not applied to one of the gauge fields, or the clover coefficient being applied twice (once by Chroma and once by QUDA).

### 5.3 Level 2: Propagator Comparison Against Reference (Recommended for Initial Validation)

Run the same inversion with both solvers and compare the propagator site-by-site. This is the gold standard for verifying your adapter.

Create a validation XML that computes the same propagator twice:

```xml
<!-- Solve with existing BiCGStab (known-good reference) -->
<elem>
    <annotation>Reference propagator (BiCGStab)</annotation>
    <n>PROPAGATOR</n>
    <Frequency>1</Frequency>
    <Param>
        <version>10</version>
        <quarkSpinType>FULL</quarkSpinType>
        <obsvP>false</obsvP>
        <numRetries>1</numRetries>
        <FermionAction>
            <FermAct>UNPRECONDITIONED_SLRC</FermAct>
            <Kappa>0.120900</Kappa>
            <clovCoeff>2.65</clovCoeff>
            <FermState>
                <n>SLIC_FERM_STATE</n>
                <n_smear>1</n_smear>
                <rho>0.1</rho>
                <orthog_dir>5</orthog_dir>
                <FermionBC>
                    <FermBC>SIMPLE_FERMBC</FermBC>
                    <boundary>1 1 1 -1</boundary>
                </FermionBC>
            </FermState>
        </FermionAction>
        <InvertParam>
            <invType>BICGSTAB_INVERTER</invType>
            <RsdBiCGStab>1.0e-14</RsdBiCGStab>
            <MaxBiCGStab>20000</MaxBiCGStab>
            <MROver>1.1</MROver>
        </InvertParam>
    </Param>
    <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>sh_source_1</source_id>
        <prop_id>prop_reference</prop_id>
    </NamedObject>
</elem>

<!-- Solve with QUDA SLRC solver (under test) -->
<elem>
    <annotation>Test propagator (QUDA SLRC)</annotation>
    <n>PROPAGATOR</n>
    <Frequency>1</Frequency>
    <Param>
        <version>10</version>
        <quarkSpinType>FULL</quarkSpinType>
        <obsvP>false</obsvP>
        <numRetries>1</numRetries>
        <FermionAction>
            <FermAct>UNPRECONDITIONED_SLRC</FermAct>
            <Kappa>0.120900</Kappa>
            <clovCoeff>2.65</clovCoeff>
            <FermState>
                <n>SLIC_FERM_STATE</n>
                <n_smear>1</n_smear>
                <rho>0.1</rho>
                <orthog_dir>5</orthog_dir>
                <FermionBC>
                    <FermBC>SIMPLE_FERMBC</FermBC>
                    <boundary>1 1 1 -1</boundary>
                </FermionBC>
            </FermState>
        </FermionAction>
        <InvertParam>
            <invType>QUDA_SLRC_CLOVER_INVERTER</invType>
            <SolverType>GCR</SolverType>
            <RsdTarget>1.0e-14</RsdTarget>
            <MaxIter>10000</MaxIter>
            <AntiPeriodicT>true</AntiPeriodicT>
            <!-- Use a direct solver (no MG) for cleanest comparison -->
            <CudaPrecision>DOUBLE</CudaPrecision>
            <CudaSloppyPrecision>DOUBLE</CudaSloppyPrecision>
            <CudaRefinementPrecision>DOUBLE</CudaRefinementPrecision>
        </InvertParam>
    </Param>
    <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>sh_source_1</source_id>
        <prop_id>prop_quda</prop_id>
    </NamedObject>
</elem>
```

Note: both solvers are set to a tighter-than-production tolerance (1e-14) and the QUDA solver uses full double precision (no mixed precision) to eliminate rounding as a source of discrepancy.

Then compare the two propagators. If your Chroma branch has an inline measurement for propagator comparison, use it. Otherwise, you can add a simple inline comparison or write a small standalone tool. The comparison should compute:

```
||prop_quda - prop_reference|| / ||prop_reference||
```

over the full lattice. This should be O(1e-12) or smaller — essentially at the level of the solver tolerance. If it's O(1e-6) or worse, there's a systematic error in the operator.

If you don't have a propagator comparison measurement, a practical alternative is to compute correlators from both propagators and compare them. This is what Level 3 does.

### 5.4 Level 3: Correlator Comparison (Production Validation)

This is the physicist's validation: compute hadron correlators using both solvers and compare the numbers. This catches everything from operator bugs to subtle boundary condition issues, and has the advantage that you can use your existing measurement infrastructure.

Run your full measurement pipeline (meson, baryon, multi-hadron FH correlators) twice: once with the BiCGStab propagator and once with the QUDA propagator, on the same configuration. Compare the output correlator files.

**What to compare:**

For each correlator C(t):

```
relative_diff(t) = |C_QUDA(t) - C_BiCGStab(t)| / |C_BiCGStab(t)|
```

**Expected results:**

- For the unperturbed propagator (which both solvers compute): the correlators should agree to ~1e-10 or better at all timeslices when using double precision, or ~1e-6 to ~1e-8 with mixed precision. The agreement may be slightly worse at large time separations where the correlator signal is exponentially suppressed and floating-point cancellations matter more.
- For FH correlators built from the same FH propagator: these should be identical (the FH inversions still use BiCGStab in both cases). This serves as a sanity check that your measurement pipeline is deterministic.
- For FH correlators where the unperturbed propagator differs: the FH solver's `initial_guess_id` points at the unperturbed prop, so a different unperturbed prop leads to a different iteration trajectory for BiCGStab. However, the final FH propagator should agree to within the solver tolerance, and the correlators should agree to ~1e-10 or better.

**A quick script to compare correlator files:**

```python
#!/usr/bin/env python3
"""Compare two Chroma correlator output files."""
import sys
import numpy as np

def read_correlator(filename):
    """Read correlator data - adapt parsing to your output format."""
    data = {}
    with open(filename) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                t = int(parts[0])
                val = complex(float(parts[1]), float(parts[2]))
                data[t] = val
    return data

ref = read_correlator(sys.argv[1])
test = read_correlator(sys.argv[2])

print(f"{'t':>4} {'|ref|':>14} {'|test|':>14} {'rel_diff':>14}")
print("-" * 50)
max_diff = 0
for t in sorted(ref.keys()):
    if t in test:
        r, q = ref[t], test[t]
        diff = abs(q - r) / abs(r) if abs(r) > 0 else abs(q - r)
        max_diff = max(max_diff, diff)
        flag = " <<<" if diff > 1e-8 else ""
        print(f"{t:4d} {abs(r):14.8e} {abs(q):14.8e} {diff:14.8e}{flag}")

print(f"\nMax relative difference: {max_diff:.3e}")
if max_diff < 1e-10:
    print("PASS: Excellent agreement (double precision level)")
elif max_diff < 1e-6:
    print("PASS: Good agreement (mixed precision level)")
else:
    print("FAIL: Significant disagreement - investigate")
```

### 5.5 Level 4: Gauge Field Verification (Debugging Aid)

If Levels 1–3 reveal a problem, you need to verify that the correct gauge fields are reaching QUDA. Add diagnostic output to your `extractSLICGaugeFields()` and `loadFieldsToQUDA()` methods:

```cpp
void LinOpSysSolverQUDASLRC::extractSLICGaugeFields(
  Handle< FermState<T,P,Q> > state_)
{
  // ... extraction code ...

  // Diagnostic: print plaquettes of both gauge fields
  // These should differ (fat links have larger plaquette)
  for (int mu = 0; mu < Nd; mu++) {
    for (int nu = mu+1; nu < Nd; nu++) {
      Double thin_plaq = sum(real(trace(
        thin_links[mu] * shift(thin_links[nu], FORWARD, mu)
        * adj(shift(thin_links[mu], FORWARD, nu)) * adj(thin_links[nu])
      ))) / Double(Layout::vol() * Nc);

      Double fat_plaq = sum(real(trace(
        fat_links[mu] * shift(fat_links[nu], FORWARD, mu)
        * adj(shift(fat_links[mu], FORWARD, nu)) * adj(fat_links[nu])
      ))) / Double(Layout::vol() * Nc);

      QDPIO::cout << "SLRC_DIAG: plaq[" << mu << "," << nu << "]"
                  << "  thin=" << thin_plaq
                  << "  fat=" << fat_plaq << std::endl;
    }
  }

  // If thin_plaq == fat_plaq for all (mu,nu), the SLIC state
  // is returning the same links for both -- this is a bug.
}
```

Also verify the clover term by computing tr(σ_μν F_μν) on a few sites and comparing against a known-good calculation.

### 5.6 Recommended Validation Sequence

When bringing up the QUDA SLRC solver for the first time:

1. **Start with a small lattice** (e.g., 4⁴ or 8⁴) for fast turnaround. Port your full action to this small lattice.

2. **Run Level 4** (gauge field diagnostics) to verify thin/fat link extraction is correct.

3. **Run Level 2** (propagator comparison) in full double precision with no multigrid. This isolates the SLRC adapter from MG complications. The propagators should agree to ~1e-12 or better.

4. **Run Level 1** (Chroma-side residual) on every solve during development. Keep this check in production code as a cheap safety net — it adds negligible overhead compared to the solve.

5. **Run Level 3** (correlator comparison) on a single production-size configuration (32³×64). This is the end-to-end validation that everything works at scale.

6. **Enable mixed precision** (single-precision sloppy) and repeat Level 3. Verify that correlators still agree to ~1e-8 or better.

7. **Enable multigrid** and repeat Level 3. The MG solver should give the same correlators as the direct solver to within solver tolerance.

8. Once validated, remove the `<Verbosity>VERBOSE</Verbosity>` and the Level 4 diagnostics (they have some overhead), but **keep the Level 1 Chroma-side residual check** permanently. It's cheap and catches regressions.

### 5.7 Common Failure Modes and Diagnosis

| Symptom | Likely cause | Fix |
|:--|:--|:--|
| QUDA residual good, Chroma residual ~O(1) | Gauge fields swapped (fat links used for clover or vice versa) | Check `extractSLICGaugeFields()` and verify with plaquette diagnostics |
| QUDA residual good, Chroma residual ~O(0.01–0.1) | Boundary conditions not applied to one gauge field | Ensure `fbc->modify()` is called on thin links before clover construction |
| QUDA residual good, Chroma residual ~O(1e-6) | Clover coefficient applied twice (Chroma + QUDA) | Set `quda_inv_param.clover_coeff = 0` and `compute_clover = false` |
| Both residuals ~O(1e-6) | Mixed precision rounding | Normal for single-precision sloppy; use double sloppy to verify |
| MG solver diverges | Null space poorly captured | Increase null vectors, tighten setup tolerance, or try 3 levels |
| Correlators agree at small t but diverge at large t | Solver tolerance too loose | Normal if relative diff is below solver tolerance; tighten `RsdTarget` if needed |
| Everything matches for one mass but not another | Mass-dependent bug in clover term | Check that κ is correctly propagated to the clover construction |

---

## Part 6: Summary of Changes

### Files to create:
1. `syssolver_linop_slrc_quda_w.h` — Header for the SLRC QUDA solver
2. `syssolver_linop_slrc_quda_w.cc` — Implementation

### Files to modify:
3. `slic_fermstate_w.h` — Add `getThinLinks()` accessor
4. Solver registration file — Register `QUDA_SLRC_CLOVER_INVERTER`
5. `Makefile.am` or `CMakeLists.txt` — Add new source files to build

### XML changes:
6. Unperturbed `<InvertParam>` blocks: `BICGSTAB_INVERTER` → `QUDA_SLRC_CLOVER_INVERTER` with MG params
7. FH `<InvertParam>` blocks: **unchanged** (stay on `BICGSTAB_INVERTER`)

### Build changes:
8. Rebuild QUDA with `QUDA_INTERFACE_QDP=ON`, `QUDA_MULTIGRID=ON`
9. Rebuild Chroma with `--with-quda=/path/to/quda/install`
