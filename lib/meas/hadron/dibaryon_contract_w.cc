/**
 * @file dibaryon_contract_w.cc
 * @brief Full 6-quark dibaryon contraction with exchange diagrams
 *
 * Optimized contraction using source diquark factorization:
 *   Diquark pairs: D(sp,sq) = (S[p] * Cg5 * S[q]^T)(sp,sq)
 *   Open (inter-baryon) pair: product of T_unpol-projected propagators
 *
 * The I=0 direct term has net coefficient 0 with the product source
 * (symmetric T_unpol coupling). All signal comes from exchange terms.
 *
 * See dibaryon_contract_w.h for full documentation.
 */

#include "meas/hadron/dibaryon_contract_w.h"

namespace Chroma
{
  namespace DibaryonContractions
  {
    // =================================================================
    // Levi-Civita tensor: 6 nonzero entries for 3 colors
    // =================================================================
    const EpsEntry eps_table[NUM_EPS_ENTRIES] = {
      {0, 1, 2, +1},
      {1, 2, 0, +1},
      {2, 0, 1, +1},
      {0, 2, 1, -1},
      {2, 1, 0, -1},
      {1, 0, 2, -1}
    };

    // =================================================================
    // 16 exchange topologies (from deuteron_wick_enum.py)
    // =================================================================
    // Coefficients recomputed for product source (T_unpol inter-baryon,
    // no Cg5 antisymmetry sign for the open pair).
    const ExchangeTopology exchange_topos[NUM_EXCHANGE_TOPOS] = {
      //                                                          open_pair
      { +1, {0,1,3}, {2,4,5}, {{0,1}, {2,3}, {4,5}}, 1 },  // #1
      { +1, {0,1,3}, {2,4,5}, {{0,1}, {2,5}, {3,4}}, 2 },  // #2
      { +1, {0,1,3}, {2,4,5}, {{0,3}, {1,2}, {4,5}}, 1 },  // #3
      { +1, {0,1,3}, {2,4,5}, {{0,3}, {1,4}, {2,5}}, 1 },  // #4
      { -1, {0,1,4}, {2,3,5}, {{0,1}, {2,3}, {4,5}}, 2 },  // #5
      { +1, {0,1,4}, {2,3,5}, {{0,5}, {1,4}, {2,3}}, 0 },  // #6
      { -1, {0,2,3}, {1,4,5}, {{0,3}, {1,2}, {4,5}}, 1 },  // #7
      { +1, {0,2,3}, {1,4,5}, {{0,3}, {1,4}, {2,5}}, 2 },  // #8
      { -1, {0,2,5}, {1,3,4}, {{0,5}, {1,2}, {3,4}}, 1 },  // #9
      { -1, {0,2,5}, {1,3,4}, {{0,5}, {1,4}, {2,3}}, 2 },  // #10
      { +1, {0,3,4}, {1,2,5}, {{0,1}, {2,5}, {3,4}}, 0 },  // #11
      { +1, {0,3,4}, {1,2,5}, {{0,3}, {1,2}, {4,5}}, 2 },  // #12
      { -1, {0,3,4}, {1,2,5}, {{0,3}, {1,4}, {2,5}}, 1 },  // #13
      { -1, {0,3,4}, {1,2,5}, {{0,5}, {1,2}, {3,4}}, 0 },  // #14
      { -1, {0,3,5}, {1,2,4}, {{0,3}, {1,2}, {4,5}}, 2 },  // #15
      { -1, {0,3,5}, {1,2,4}, {{0,3}, {1,4}, {2,5}}, 2 },  // #16
    };

    // Slot -> flavor index: 0=up, 1=down
    static const int slot_flavor[NUM_SLOTS] = {0, 1, 0, 1, 0, 1};


    // =================================================================
    // Exchange contraction: optimized with diquark factorization
    // =================================================================
    LatticeSpinMatrix dibaryonExchange(
        const LatticePropagator& S_u,
        const LatticePropagator& S_d,
        const SpinMatrix& Cg5,
        const SpinMatrix& T_unpol)
    {
#if QDP_NC == 3
      START_CODE();

      // ---------------------------------------------------------
      // Step 1: Extract all color blocks as LatticeSpinMatrix
      // S_color[flavor][a_snk][a_src]
      // flavor 0=up, 1=down; a,b = 0,1,2 (color indices)
      // ---------------------------------------------------------
      LatticeSpinMatrix S_color[2][Nc][Nc];
      for (int a = 0; a < Nc; a++)
        for (int b = 0; b < Nc; b++) {
          S_color[0][a][b] = peekColor(S_u, a, b);
          S_color[1][a][b] = peekColor(S_d, a, b);
        }

      // ---------------------------------------------------------
      // Step 2a: Precompute Cg5 * transpose(S) for diquark pairs
      // CgST[f][a][b] = Cg5 * transposeSpin(S_color[f][a][b])
      // ---------------------------------------------------------
      LatticeSpinMatrix CgST[2][Nc][Nc];
      for (int f = 0; f < 2; f++)
        for (int a = 0; a < Nc; a++)
          for (int b = 0; b < Nc; b++)
            CgST[f][a][b] = Cg5 * transposeSpin(S_color[f][a][b]);

      // ---------------------------------------------------------
      // Step 2b: Precompute T_unpol-projected propagator sums
      // for the inter-baryon (open) pair.
      //
      // Tv[f][a][b][s] = sum_beta S_color[f][a][b](s, beta) * T_unpol(beta, beta)
      //
      // Since T_unpol = diag(1,1,0,0) in DeGrand-Rossi, this sums
      // the first two source spin components for each sink spin.
      // ---------------------------------------------------------
      LatticeComplex Tv[2][Nc][Nc][Nd];
      for (int f = 0; f < 2; f++)
        for (int a = 0; a < Nc; a++)
          for (int b = 0; b < Nc; b++) {
            LatticeSpinMatrix ST = S_color[f][a][b] * T_unpol;
            for (int s = 0; s < Nd; s++) {
              Tv[f][a][b][s] = zero;
              for (int beta = 0; beta < Nd; beta++)
                Tv[f][a][b][s] += peekSpin(ST, s, beta);
            }
          }

      // ---------------------------------------------------------
      // Step 3: Precompute nonzero Cg5 entries for sink contraction
      // ---------------------------------------------------------
      struct Cg5Entry { int r, c; Complex val; };
      Cg5Entry cg5_nz[Nd*Nd];
      int n_cg5_nz = 0;
      for (int a = 0; a < Nd; a++)
        for (int b = 0; b < Nd; b++) {
          Complex v = peekSpin(Cg5, a, b);
          if (toBool(real(v) != 0.0) || toBool(imag(v) != 0.0)) {
            cg5_nz[n_cg5_nz].r = a;
            cg5_nz[n_cg5_nz].c = b;
            cg5_nz[n_cg5_nz].val = v;
            n_cg5_nz++;
          }
        }

      // ---------------------------------------------------------
      // Step 4: Accumulate results into LatticeComplex array
      // result_elem[s2][s5] accumulates the exchange contribution
      // ---------------------------------------------------------
      LatticeComplex result_elem[Nd][Nd];
      for (int s2 = 0; s2 < Nd; s2++)
        for (int s5 = 0; s5 < Nd; s5++)
          result_elem[s2][s5] = zero;

      // ---------------------------------------------------------
      // Step 5: Main loop over topologies and color assignments
      // ---------------------------------------------------------
      for (int t = 0; t < NUM_EXCHANGE_TOPOS; t++)
      {
        const ExchangeTopology& topo = exchange_topos[t];

        // Map slots to source epsilon positions
        int slot_eps_pos[NUM_SLOTS];
        int slot_eps_grp[NUM_SLOTS];  // 0=A, 1=B
        for (int s = 0; s < NUM_SLOTS; s++) {
          slot_eps_pos[s] = -1;
          slot_eps_grp[s] = -1;
        }
        for (int p = 0; p < 3; p++) {
          slot_eps_grp[topo.eps_A[p]] = 0;
          slot_eps_pos[topo.eps_A[p]] = p;
          slot_eps_grp[topo.eps_B[p]] = 1;
          slot_eps_pos[topo.eps_B[p]] = p;
        }

        // Source color loop (36 terms: eps_A x eps_B)
        for (int eA = 0; eA < NUM_EPS_ENTRIES; eA++)
        {
          int cA[3] = {eps_table[eA].i, eps_table[eA].j, eps_table[eA].k};
          int signA = eps_table[eA].sign;

          for (int eB = 0; eB < NUM_EPS_ENTRIES; eB++)
          {
            int cB[3] = {eps_table[eB].i, eps_table[eB].j, eps_table[eB].k};
            int signB = eps_table[eB].sign;

            // Source color for each slot
            int c_src[NUM_SLOTS];
            for (int s = 0; s < NUM_SLOTS; s++) {
              c_src[s] = (slot_eps_grp[s] == 0) ?
                cA[slot_eps_pos[s]] : cB[slot_eps_pos[s]];
            }

            // Sink color loop (36 terms: eps_P x eps_N)
            for (int eP = 0; eP < NUM_EPS_ENTRIES; eP++)
            {
              int aP[3] = {eps_table[eP].i, eps_table[eP].j, eps_table[eP].k};
              int signP = eps_table[eP].sign;

              for (int eN = 0; eN < NUM_EPS_ENTRIES; eN++)
              {
                int aN[3] = {eps_table[eN].i, eps_table[eN].j, eps_table[eN].k};
                int signN = eps_table[eN].sign;

                int a_snk[NUM_SLOTS];
                a_snk[0] = aP[0]; a_snk[1] = aP[1]; a_snk[2] = aP[2];
                a_snk[3] = aN[0]; a_snk[4] = aN[1]; a_snk[5] = aN[2];

                int color_sign = topo.coeff * signA * signB * signP * signN;

                // Build source pair elements as LatticeComplex
                // For diquark pairs: D(sp,sq) = (S[p] * Cg5 * S[q]^T)(sp,sq)
                // For open pair: De[k][r][c] = Tv_p[r] * Tv_q[c]
                //   (T_unpol projection on each propagator independently)
                LatticeComplex De[3][Nd][Nd];

                for (int k = 0; k < 3; k++) {
                  int pk = topo.cg5[k][0];
                  int qk = topo.cg5[k][1];

                  if (k == topo.open_pair) {
                    // Inter-baryon pair: T_unpol projection (no Cg5)
                    for (int r = 0; r < Nd; r++)
                      for (int c = 0; c < Nd; c++)
                        De[k][r][c] =
                          Tv[slot_flavor[pk]][a_snk[pk]][c_src[pk]][r] *
                          Tv[slot_flavor[qk]][a_snk[qk]][c_src[qk]][c];
                  } else {
                    // Diquark pair: Cg5 contraction
                    LatticeSpinMatrix D =
                      S_color[slot_flavor[pk]][a_snk[pk]][c_src[pk]]
                      * CgST[slot_flavor[qk]][a_snk[qk]][c_src[qk]];
                    for (int r = 0; r < Nd; r++)
                      for (int c = 0; c < Nd; c++)
                        De[k][r][c] = peekSpin(D, r, c);
                  }
                }

                // Contract sink Cg5 (proton diquark s0,s1; neutron diquark s3,s4)
                // and accumulate into result_elem[s2][s5]
                for (int s2 = 0; s2 < Nd; s2++) {
                  for (int s5 = 0; s5 < Nd; s5++) {

                    LatticeComplex acc;
                    acc = zero;

                    // Sum over sink Cg5 pairs using nonzero entries
                    for (int iP = 0; iP < n_cg5_nz; iP++) {
                      int s0 = cg5_nz[iP].r;
                      int s1 = cg5_nz[iP].c;

                      for (int iN = 0; iN < n_cg5_nz; iN++) {
                        int s3 = cg5_nz[iN].r;
                        int s4 = cg5_nz[iN].c;

                        // All 6 sink spins determined
                        int s[NUM_SLOTS] = {s0, s1, s2, s3, s4, s5};

                        // Product of 3 pair elements
                        acc += cg5_nz[iP].val * cg5_nz[iN].val *
                          De[0][s[topo.cg5[0][0]]][s[topo.cg5[0][1]]] *
                          De[1][s[topo.cg5[1][0]]][s[topo.cg5[1][1]]] *
                          De[2][s[topo.cg5[2][0]]][s[topo.cg5[2][1]]];
                      } // iN
                    } // iP

                    result_elem[s2][s5] +=
                      cmplx(Real(color_sign), Real(0)) * acc;

                  } // s5
                } // s2

              } // eN
            } // eP
          } // eB
        } // eA
      } // topology

      // ---------------------------------------------------------
      // Step 6: Pack result_elem into LatticeSpinMatrix
      // ---------------------------------------------------------
      LatticeSpinMatrix result;
      result = zero;
      for (int s2 = 0; s2 < Nd; s2++)
        for (int s5 = 0; s5 < Nd; s5++)
          pokeSpin(result, result_elem[s2][s5], s2, s5);

      END_CODE();
      return result;

#else
      LatticeSpinMatrix a;
      a = zero;
      return a;
#endif
    }


    // =================================================================
    // Full dibaryon: direct + exchange
    // =================================================================
    LatticeSpinMatrix dibaryonFull(
        const LatticePropagator& S_u,
        const LatticePropagator& S_d,
        const SpinMatrix& Cg5,
        const SpinMatrix& T_unpol)
    {
#if QDP_NC == 3
      START_CODE();

      // Direct term (topology 0) has net coefficient 0 in I=0 channel
      // with product source (T_unpol projection). The pn and np direct
      // terms cancel because the symmetric T_unpol coupling cannot
      // distinguish the two orderings. The I=0 signal comes entirely
      // from quark-exchange topologies.

      // All 16 exchange topologies
      LatticeSpinMatrix exchange = dibaryonExchange(S_u, S_d, Cg5, T_unpol);

      END_CODE();
      return exchange;

#else
      LatticeSpinMatrix a;
      a = zero;
      return a;
#endif
    }

  } // namespace DibaryonContractions

} // namespace Chroma
