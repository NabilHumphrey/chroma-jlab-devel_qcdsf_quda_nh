// -*- C++ -*-
/**
 * @file nucleon_block_w.h
 * @brief Non-degenerate nucleon contraction functions
 *
 * Provides correct nucleon 2-point contraction with separate quark propagators
 * for each flavor. Matches the standard nucl2pt() from baryon_w.cc but supports
 * non-degenerate up/down quarks as needed for Feynman-Hellmann calculations.
 *
 * Also provides the "nucleon block" -- the partially-contracted baryon with
 * sink spin indices left open as a LatticeSpinMatrix. This is the building
 * block for multi-hadron contractions (dibaryon, etc.).
 */

#ifndef __nucleon_block_w_h__
#define __nucleon_block_w_h__

#include "chromabase.h"

namespace Chroma
{
  namespace NucleonBlocks
  {
    //! Nucleon block with sink spin indices open
    /*!
     * Computes the color-contracted nucleon with sink spin indices left
     * as a LatticeSpinMatrix, ready for spin projection or multi-hadron
     * contraction. Includes BOTH Wick contraction terms.
     *
     * Interpolator: epsilon_{abc} (q1^a sp q2^b) q3^c
     *
     * Result: B_{alpha,beta}(x) where alpha = sink spin, beta = source spin
     *
     * B = traceColor(q3 * traceSpin(diquark)) + traceColor(q3 * diquark)
     *
     * where diquark = quarkContract13(q1 * sp, sp * q2)
     *
     * \param q1   First diquark quark propagator (Read)
     * \param q2   Second diquark quark propagator (Read)
     * \param q3   Spectator quark propagator (Read)
     * \param sp   Diquark spin matrix, e.g. Cg5 (Read)
     *
     * \return LatticeSpinMatrix with sink/source spin indices open
     */
    LatticeSpinMatrix nucleonBlock(
        const LatticePropagator& q1,
        const LatticePropagator& q2,
        const LatticePropagator& q3,
        const SpinMatrix& sp);


    //! Non-degenerate nucleon 2-point correlator (fully traced)
    /*!
     * Correct non-degenerate version of nucl2pt() from baryon_w.cc.
     * Computes the nucleon block and traces with the spin projector T.
     *
     * Interpolator: epsilon_{abc} (q1^a sp q2^b) q3^c
     *
     * For a proton: q1 = u, q2 = d, q3 = u, sp = BaryonSpinMats::Cg5()
     * For a neutron: q1 = d, q2 = u, q3 = d, sp = BaryonSpinMats::Cg5()
     *
     * \param q1   First diquark quark propagator (Read)
     * \param q2   Second diquark quark propagator (Read)
     * \param q3   Spectator quark propagator (Read)
     * \param T    Spin projector matrix, e.g. Tunpol, Tmixed (Read)
     * \param sp   Diquark spin matrix, e.g. Cg5 (Read)
     *
     * \return LatticeComplex correlator field ready for momentum projection
     */
    LatticeComplex nucl2pt(
        const LatticePropagator& q1,
        const LatticePropagator& q2,
        const LatticePropagator& q3,
        const SpinMatrix& T,
        const SpinMatrix& sp);

  } // namespace NucleonBlocks

} // namespace Chroma

#endif
