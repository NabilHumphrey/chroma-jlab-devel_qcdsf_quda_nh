/**
 * @file dibaryon_contract_w.cc
 * @brief Full 6-quark dibaryon contraction with exchange diagrams
 *
 * Optimized contraction using source-Cg5 diquark factorization:
 *   D_{p,q}(sp,sq) = sum_{bp,bq} Cg5(bp,bq) * S[p](sp,bp) * S[q](sq,bq)
 *                  = (S[p] * Cg5 * S[q]^T)(sp,sq)
 *
 * This reduces the 6 source spin loops to 3 spin matrix multiplies.
 *
 * See dibaryon_contract_w.h for full documentation.
 */

#include "meas/hadron/dibaryon_contract_w.h"
#include "meas/hadron/nucleon_block_w.h"

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
    const ExchangeTopology exchange_topos[NUM_EXCHANGE_TOPOS] = {
      { -1, {0,1,3}, {2,4,5}, {{0,1}, {2,3}, {4,5}} },  // #1
      { +1, {0,1,3}, {2,4,5}, {{0,1}, {2,5}, {3,4}} },  // #2
      { +1, {0,1,3}, {2,4,5}, {{0,3}, {1,2}, {4,5}} },  // #3
      { +1, {0,1,3}, {2,4,5}, {{0,3}, {1,4}, {2,5}} },  // #4
      { -1, {0,1,4}, {2,3,5}, {{0,1}, {2,3}, {4,5}} },  // #5
      { +1, {0,1,4}, {2,3,5}, {{0,5}, {1,4}, {2,3}} },  // #6
      { +1, {0,2,3}, {1,4,5}, {{0,3}, {1,2}, {4,5}} },  // #7
      { +1, {0,2,3}, {1,4,5}, {{0,3}, {1,4}, {2,5}} },  // #8
      { -1, {0,2,5}, {1,3,4}, {{0,5}, {1,2}, {3,4}} },  // #9
      { +1, {0,2,5}, {1,3,4}, {{0,5}, {1,4}, {2,3}} },  // #10
      { +1, {0,3,4}, {1,2,5}, {{0,1}, {2,5}, {3,4}} },  // #11
      { +1, {0,3,4}, {1,2,5}, {{0,3}, {1,2}, {4,5}} },  // #12
      { +1, {0,3,4}, {1,2,5}, {{0,3}, {1,4}, {2,5}} },  // #13
      { -1, {0,3,4}, {1,2,5}, {{0,5}, {1,2}, {3,4}} },  // #14
      { +1, {0,3,5}, {1,2,4}, {{0,3}, {1,2}, {4,5}} },  // #15
      { +1, {0,3,5}, {1,2,4}, {{0,3}, {1,4}, {2,5}} },  // #16
    };

    // Slot -> flavor index: 0=up, 1=down
    static const int slot_flavor[NUM_SLOTS] = {0, 1, 0, 1, 0, 1};


    // =================================================================
    // Exchange contraction: optimized with diquark factorization
    // =================================================================
    LatticeSpinMatrix dibaryonExchange(
        const LatticePropagator& S_u,
        const LatticePropagator& S_d,
        const SpinMatrix& Cg5)
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
      // Step 2: Precompute Cg5 * transpose(S) for all color blocks
      //
      // Diquark: D(sp,sq) = sum_{bp,bq} S[p](sp,bp) * Cg5(bp,bq) * S[q](sq,bq)
      //        = sum_bp S[p](sp,bp) * (Cg5 * S[q]^T)(bp,sq)
      //        = (S[p] * Cg5 * S[q]^T)(sp,sq)
      //
      // CgST[f][a][b] = Cg5 * transposeSpin(S_color[f][a][b])
      // Then: D = S_color[fp][ap][cp] * CgST[fq][aq][cq]
      // ---------------------------------------------------------
      LatticeSpinMatrix CgST[2][Nc][Nc];
      for (int f = 0; f < 2; f++)
        for (int a = 0; a < Nc; a++)
          for (int b = 0; b < Nc; b++)
            CgST[f][a][b] = Cg5 * transposeSpin(S_color[f][a][b]);

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

                // Build 3 source-Cg5 diquark matrices
                // D[k](s_pk, s_qk) = (S[pk] * CgST[qk])(s_pk, s_qk)
                // = sum_{bp,bq} S[pk](s_pk,bp) * Cg5(bp,bq) * S[qk](s_qk,bq)
                LatticeSpinMatrix D[3];
                for (int k = 0; k < 3; k++) {
                  int pk = topo.cg5[k][0];
                  int qk = topo.cg5[k][1];
                  D[k] = S_color[slot_flavor[pk]][a_snk[pk]][c_src[pk]]
                       * CgST[slot_flavor[qk]][a_snk[qk]][c_src[qk]];
                }

                // Extract diquark elements as LatticeComplex
                LatticeComplex De[3][Nd][Nd];
                for (int k = 0; k < 3; k++)
                  for (int r = 0; r < Nd; r++)
                    for (int c = 0; c < Nd; c++)
                      De[k][r][c] = peekSpin(D[k], r, c);

                // Contract sink Cg5 (proton diquark s0,s1; neutron diquark s3,s4)
                // and accumulate into result_elem[s2][s5]
                //
                // result(s2,s5) += color_sign
                //   * sum_{s0,s1} Cg5(s0,s1) * sum_{s3,s4} Cg5(s3,s4)
                //   * D0(s_{p0},s_{q0}) * D1(s_{p1},s_{q1}) * D2(s_{p2},s_{q2})
                //
                // where D[k] row = s_{cg5[k][0]}, col = s_{cg5[k][1]}
                // and s0..s5 are the sink spins for slots 0..5

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

                        // Product of 3 diquark elements
                        // D[k] indexed by (s[cg5[k][0]], s[cg5[k][1]])
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
        const SpinMatrix& Cg5)
    {
#if QDP_NC == 3
      START_CODE();

      // Direct term (topology 0, coeff +2):
      //   B_P(s2,beta2) * Cg5(beta2,beta5) * B_N(s5,beta5)
      //   = (B_P * Cg5 * B_N^T)_{s2,s5}
      LatticeSpinMatrix B_P = NucleonBlocks::nucleonBlock(S_u, S_d, S_u, Cg5);
      LatticeSpinMatrix B_N = NucleonBlocks::nucleonBlock(S_d, S_u, S_d, Cg5);

      LatticeSpinMatrix direct = 2 * B_P * Cg5 * transposeSpin(B_N);

      // Exchange terms
      LatticeSpinMatrix exchange = dibaryonExchange(S_u, S_d, Cg5);

      END_CODE();
      return direct + exchange;

#else
      LatticeSpinMatrix a;
      a = zero;
      return a;
#endif
    }

  } // namespace DibaryonContractions

} // namespace Chroma
