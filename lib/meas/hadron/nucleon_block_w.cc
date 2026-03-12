/**
 * @file nucleon_block_w.cc
 * @brief Non-degenerate nucleon contraction functions
 *
 * Implements the correct nucleon 2-point contraction with separate quark
 * propagators, matching the structure of nucl2pt() in baryon_w.cc.
 */

#include "meas/hadron/nucleon_block_w.h"

namespace Chroma
{
  namespace NucleonBlocks
  {
    //! Nucleon block with sink spin indices open
    LatticeSpinMatrix nucleonBlock(
        const LatticePropagator& q1,
        const LatticePropagator& q2,
        const LatticePropagator& q3,
        const SpinMatrix& sp)
    {
#if QDP_NC == 3
      // Form diquark: epsilon contraction of q1 and q2 with spin structure sp
      // This matches baryon_w.cc nucl2pt: quarkContract13(prop * sp, sp * prop)
      // but with separate propagators for non-degenerate quarks
      LatticePropagator di_quark = quarkContract13(q1 * sp, sp * q2);

      // Both Wick contraction terms, with spin indices left open:
      //   Term 1 (direct):   traceColor(q3 * traceSpin(di_quark))
      //   Term 2 (exchange): traceColor(q3 * di_quark)
      return LatticeSpinMatrix(traceColor(q3 * traceSpin(di_quark))
                             + traceColor(q3 * di_quark));
#else
      LatticeSpinMatrix a;
      a = zero;
      return a;
#endif
    }


    //! Non-degenerate nucleon 2-point correlator (fully traced)
    LatticeComplex nucl2pt(
        const LatticePropagator& q1,
        const LatticePropagator& q2,
        const LatticePropagator& q3,
        const SpinMatrix& T,
        const SpinMatrix& sp)
    {
#if QDP_NC == 3
      // Form diquark with correct spin structure on both sides
      LatticePropagator di_quark = quarkContract13(q1 * sp, sp * q2);

      // Both Wick contraction terms, traced with spin projector T
      // This matches baryon_w.cc nucl2pt exactly:
      //   trace(T * traceColor(prop * traceSpin(di_quark)))
      //   + trace(T * traceColor(prop * di_quark))
      return LatticeComplex(trace(T * traceColor(q3 * traceSpin(di_quark)))
                          + trace(T * traceColor(q3 * di_quark)));
#else
      LatticeComplex a;
      a = zero;
      return a;
#endif
    }

  } // namespace NucleonBlocks

} // namespace Chroma
