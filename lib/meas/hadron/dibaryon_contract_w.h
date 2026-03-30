// -*- C++ -*-
/**
 * @file dibaryon_contract_w.h
 * @brief Full 6-quark dibaryon contraction with exchange diagrams
 *
 * Implements the 16 exchange Wick contraction topologies for the I=0
 * deuteron (proton-neutron) system with product source (independent
 * T_unpol-projected baryons). The direct term (topology 0) has net
 * coefficient 0 in the I=0 channel — the pn and np contributions
 * cancel with symmetric inter-baryon coupling. All I=0 signal comes
 * from quark-exchange topologies.
 *
 * Each exchange topology has 3 source quark pairings: 2 are Cg5
 * diquark contractions, 1 is the "open" (inter-baryon) pair where
 * each propagator is independently projected through T_unpol.
 *
 * Slot definitions (sink side, fixed for all topologies):
 *   slot 0: P.dq1 -> S_u   (proton diquark 1)
 *   slot 1: P.dq2 -> S_d   (proton diquark 2)
 *   slot 2: P.spec -> S_u  (proton spectator)
 *   slot 3: N.dq1 -> S_d   (neutron diquark 1)
 *   slot 4: N.dq2 -> S_u   (neutron diquark 2)
 *   slot 5: N.spec -> S_d  (neutron spectator)
 *
 * Sink epsilons always: eps_P = (0,1,2), eps_N = (3,4,5)
 * Sink Cg5: proton diquark (0,1), neutron diquark (3,4)
 * Open sink spins: s2 (proton spectator), s5 (neutron spectator)
 *
 * Reference: deuteron_wick_enum.py, dibaryon_contractions_derivation.py
 */

#ifndef __dibaryon_contract_w_h__
#define __dibaryon_contract_w_h__

#include "chromabase.h"

namespace Chroma
{
  namespace DibaryonContractions
  {
    //! Number of quark slots in the dibaryon system
    static const int NUM_SLOTS = 6;

    //! Number of exchange topologies (direct term handled separately)
    static const int NUM_EXCHANGE_TOPOS = 16;

    //! Topology descriptor for exchange terms
    struct ExchangeTopology {
      int coeff;       //! +1 or -1
      int eps_A[3];    //! Source epsilon group A (slot indices)
      int eps_B[3];    //! Source epsilon group B (slot indices)
      int cg5[3][2];   //! Source Cg5 pairings: 3 pairs of (slot_p, slot_q)
      int open_pair;   //! Which of the 3 pairs is inter-baryon (0,1,2)
    };

    //! Table of 16 exchange topologies
    extern const ExchangeTopology exchange_topos[NUM_EXCHANGE_TOPOS];

    //! Levi-Civita tensor for 3 colors: nonzero entries
    /*! Each entry is {i, j, k, sign} where eps_{ijk} = sign */
    struct EpsEntry {
      int i, j, k;
      int sign;
    };
    static const int NUM_EPS_ENTRIES = 6;
    extern const EpsEntry eps_table[NUM_EPS_ENTRIES];

    //! Compute exchange contribution to dibaryon correlator (spin-resolved)
    /*!
     * Returns a LatticeSpinMatrix with indices (s2_snk, s5_snk) —
     * the proton and neutron spectator sink spins — left open.
     * The caller applies the spin-1 projection / polarization.
     *
     * This sums all 16 exchange topologies using per-site explicit
     * color and spin loops. The computation is embarrassingly parallel
     * across lattice sites (QDP++ handles the parallelism).
     *
     * \param S_u     Up quark propagator (Read)
     * \param S_d     Down quark propagator (Read)
     * \param Cg5     Diquark spin matrix C*gamma_5 (Read)
     * \param T_unpol Unpolarized projection (1+gamma_4)/2 (Read)
     *
     * \return LatticeSpinMatrix B_{s2,s5}(x) exchange contribution
     */
    LatticeSpinMatrix dibaryonExchange(
        const LatticePropagator& S_u,
        const LatticePropagator& S_d,
        const SpinMatrix& Cg5,
        const SpinMatrix& T_unpol);

    //! Reference (slow) exchange contraction using lattice-wide QDP++ ops
    /*!
     * Identical algorithm to dibaryonExchange but using LatticeComplex
     * operations in nested scalar loops. Kept for verification only —
     * produces the same result but is ~1000x slower.
     */
    LatticeSpinMatrix dibaryonExchangeReference(
        const LatticePropagator& S_u,
        const LatticePropagator& S_d,
        const SpinMatrix& Cg5,
        const SpinMatrix& T_unpol);

    //! Compute FULL dibaryon correlator (spin-resolved)
    /*!
     * Returns a LatticeSpinMatrix from all 16 exchange topologies.
     * (Direct term has net coeff 0 in I=0 and is not included.)
     *
     * \param S_u     Up quark propagator (Read)
     * \param S_d     Down quark propagator (Read)
     * \param Cg5     Diquark spin matrix C*gamma_5 (Read)
     * \param T_unpol Unpolarized projection (1+gamma_4)/2 (Read)
     *
     * \return LatticeSpinMatrix B_{s2,s5}(x) full dibaryon contraction
     */
    LatticeSpinMatrix dibaryonFull(
        const LatticePropagator& S_u,
        const LatticePropagator& S_d,
        const SpinMatrix& Cg5,
        const SpinMatrix& T_unpol);

  } // namespace DibaryonContractions

} // namespace Chroma

#endif
