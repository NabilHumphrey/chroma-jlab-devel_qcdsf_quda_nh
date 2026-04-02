/**
 * @file dibaryon_contract_w.cc
 * @brief Full 6-quark dibaryon contraction with exchange diagrams
 *
 * Optimized contraction using per-site scalar arithmetic.
 * All 16 exchange topology contractions are performed in a single
 * sweep over the lattice, extracting propagator elements at each
 * site and doing all color/spin index work with scalar Complex.
 *
 * The I=0 direct term has net coefficient 0 with the product source
 * (symmetric T_unpol coupling). All signal comes from exchange terms.
 *
 * See dibaryon_contract_w.h for full documentation.
 */

#include "meas/hadron/dibaryon_contract_w.h"
#include <complex>

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
    // Exchange contraction: hybrid GPU precompute + CPU per-site loop
    //
    // Phase 1: Use lattice-wide QDP++ ops to extract color blocks,
    //          precompute CgST and Tv (runs on GPU via QDP-JIT).
    // Phase 2: Extract spin elements from color blocks into flat
    //          CPU arrays via peekSpin (lattice-wide, GPU-friendly).
    // Phase 3: Per-site scalar contraction loop on CPU.
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
      // Precompute nonzero Cg5 entries (scalar, tiny)
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

      // Precompute nonzero T_unpol diagonal entries
      struct TEntry { int idx; Real val; };
      TEntry tunpol_nz[Nd];
      int n_tunpol_nz = 0;
      for (int b = 0; b < Nd; b++) {
        Complex v = peekSpin(T_unpol, b, b);
        if (toBool(real(v) != 0.0)) {
          tunpol_nz[n_tunpol_nz].idx = b;
          tunpol_nz[n_tunpol_nz].val = real(v);
          n_tunpol_nz++;
        }
      }

      // Precompute slot->eps mappings for all topologies
      struct TopoMap {
        int slot_eps_grp[NUM_SLOTS];
        int slot_eps_pos[NUM_SLOTS];
      };
      TopoMap topo_maps[NUM_EXCHANGE_TOPOS];
      for (int t = 0; t < NUM_EXCHANGE_TOPOS; t++) {
        const ExchangeTopology& topo = exchange_topos[t];
        for (int s = 0; s < NUM_SLOTS; s++) {
          topo_maps[t].slot_eps_grp[s] = -1;
          topo_maps[t].slot_eps_pos[s] = -1;
        }
        for (int p = 0; p < 3; p++) {
          topo_maps[t].slot_eps_grp[topo.eps_A[p]] = 0;
          topo_maps[t].slot_eps_pos[topo.eps_A[p]] = p;
          topo_maps[t].slot_eps_grp[topo.eps_B[p]] = 1;
          topo_maps[t].slot_eps_pos[topo.eps_B[p]] = p;
        }
      }

      // ---------------------------------------------------------
      // Phase 1: Lattice-wide precomputation (GPU-friendly)
      // Extract color blocks and precompute CgST, Tv
      // ---------------------------------------------------------
      LatticeSpinMatrix S_color[2][Nc][Nc];
      for (int a = 0; a < Nc; a++)
        for (int b = 0; b < Nc; b++) {
          S_color[0][a][b] = peekColor(S_u, a, b);
          S_color[1][a][b] = peekColor(S_d, a, b);
        }

      LatticeSpinMatrix CgST_lat[2][Nc][Nc];
      for (int f = 0; f < 2; f++)
        for (int a = 0; a < Nc; a++)
          for (int b = 0; b < Nc; b++)
            CgST_lat[f][a][b] = Cg5 * transposeSpin(S_color[f][a][b]);

      LatticeComplex Tv_lat[2][Nc][Nc][Nd];
      for (int f = 0; f < 2; f++)
        for (int a = 0; a < Nc; a++)
          for (int b = 0; b < Nc; b++) {
            LatticeSpinMatrix ST = S_color[f][a][b] * T_unpol;
            for (int s = 0; s < Nd; s++) {
              Tv_lat[f][a][b][s] = zero;
              for (int beta = 0; beta < Nd; beta++)
                Tv_lat[f][a][b][s] += peekSpin(ST, s, beta);
            }
          }

      // ---------------------------------------------------------
      // Phase 2: Extract spin elements to LatticeComplex arrays
      // (each peekSpin is one lattice-wide GPU op — fast)
      // ---------------------------------------------------------
      LatticeComplex SC[2][Nc][Nc][Nd][Nd];
      LatticeComplex CG[2][Nc][Nc][Nd][Nd];
      for (int f = 0; f < 2; f++)
        for (int a = 0; a < Nc; a++)
          for (int b = 0; b < Nc; b++)
            for (int s1 = 0; s1 < Nd; s1++)
              for (int s2 = 0; s2 < Nd; s2++) {
                SC[f][a][b][s1][s2] = peekSpin(S_color[f][a][b], s1, s2);
                CG[f][a][b][s1][s2] = peekSpin(CgST_lat[f][a][b], s1, s2);
              }

      // ---------------------------------------------------------
      // Phase 3: Per-site contraction loop
      // Access site data via .elem(site) on LatticeComplex objects.
      // LatticeComplex elements are small (1 complex per site) so
      // .elem(site) just reads from the host buffer after the
      // lattice-wide peekSpin ops above have synced data to host.
      // ---------------------------------------------------------
      LatticeSpinMatrix result;
      result = zero;

      int num_sites = Layout::sitesOnNode();

      for (int site = 0; site < num_sites; site++)
      {
        // Extract site-local scalar data from precomputed lattice arrays
        // These are OScalar Complex values (2 doubles each)
        typedef std::complex<double> cplx;
        cplx sc[2][Nc][Nc][Nd][Nd];
        cplx cg[2][Nc][Nc][Nd][Nd];
        cplx tv[2][Nc][Nc][Nd];

        for (int f = 0; f < 2; f++)
          for (int a = 0; a < Nc; a++)
            for (int b = 0; b < Nc; b++) {
              for (int s1 = 0; s1 < Nd; s1++)
                for (int s2 = 0; s2 < Nd; s2++) {
                  // LatticeComplex.elem(site) -> PScalar<PScalar<RComplex<REAL>>>
                  // .elem().elem() -> RComplex<REAL>
                  // .real()/.imag() -> REAL
                  sc[f][a][b][s1][s2] = cplx(
                    toDouble(SC[f][a][b][s1][s2].elem(site).elem().elem().real()),
                    toDouble(SC[f][a][b][s1][s2].elem(site).elem().elem().imag()));
                  cg[f][a][b][s1][s2] = cplx(
                    toDouble(CG[f][a][b][s1][s2].elem(site).elem().elem().real()),
                    toDouble(CG[f][a][b][s1][s2].elem(site).elem().elem().imag()));
                }
              for (int s = 0; s < Nd; s++) {
                tv[f][a][b][s] = cplx(
                  toDouble(Tv_lat[f][a][b][s].elem(site).elem().elem().real()),
                  toDouble(Tv_lat[f][a][b][s].elem(site).elem().elem().imag()));
              }
            }

        // Extract scalar Cg5 values (same for all sites)
        // Only needed on first iteration but cheap enough to redo
        cplx cg5_vals[Nd*Nd];
        for (int i = 0; i < n_cg5_nz; i++)
          cg5_vals[i] = cplx(toDouble(real(cg5_nz[i].val)),
                             toDouble(imag(cg5_nz[i].val)));

        // Accumulate result for this site
        cplx result_local[Nd][Nd];
        for (int s2 = 0; s2 < Nd; s2++)
          for (int s5 = 0; s5 < Nd; s5++)
            result_local[s2][s5] = cplx(0.0, 0.0);

        // Loop over all 16 exchange topologies
        for (int t = 0; t < NUM_EXCHANGE_TOPOS; t++)
        {
          const ExchangeTopology& topo = exchange_topos[t];
          const TopoMap& tm = topo_maps[t];

          for (int eA = 0; eA < NUM_EPS_ENTRIES; eA++)
          {
            int cA[3] = {eps_table[eA].i, eps_table[eA].j, eps_table[eA].k};
            int signA = eps_table[eA].sign;

            for (int eB = 0; eB < NUM_EPS_ENTRIES; eB++)
            {
              int cB[3] = {eps_table[eB].i, eps_table[eB].j, eps_table[eB].k};
              int signB = eps_table[eB].sign;

              int c_src[NUM_SLOTS];
              for (int s = 0; s < NUM_SLOTS; s++) {
                c_src[s] = (tm.slot_eps_grp[s] == 0) ?
                  cA[tm.slot_eps_pos[s]] : cB[tm.slot_eps_pos[s]];
              }

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

                  double color_sign = topo.coeff * signA * signB * signP * signN;

                  // Build 3 source pair elements: De[k][r][c]
                  cplx De[3][Nd][Nd];

                  for (int k = 0; k < 3; k++) {
                    int pk = topo.cg5[k][0];
                    int qk = topo.cg5[k][1];

                    if (k == topo.open_pair) {
                      for (int r = 0; r < Nd; r++)
                        for (int c = 0; c < Nd; c++)
                          De[k][r][c] =
                            tv[slot_flavor[pk]][a_snk[pk]][c_src[pk]][r] *
                            tv[slot_flavor[qk]][a_snk[qk]][c_src[qk]][c];
                    } else {
                      int fp = slot_flavor[pk], fq = slot_flavor[qk];
                      int ap = a_snk[pk], bp = c_src[pk];
                      int aq = a_snk[qk], bq = c_src[qk];
                      for (int r = 0; r < Nd; r++)
                        for (int c = 0; c < Nd; c++) {
                          cplx sum(0.0, 0.0);
                          for (int kk = 0; kk < Nd; kk++)
                            sum += sc[fp][ap][bp][r][kk] * cg[fq][aq][bq][kk][c];
                          De[k][r][c] = sum;
                        }
                    }
                  }

                  // Contract sink Cg5 and accumulate
                  for (int s2 = 0; s2 < Nd; s2++) {
                    for (int s5 = 0; s5 < Nd; s5++) {

                      cplx acc(0.0, 0.0);

                      for (int iP = 0; iP < n_cg5_nz; iP++) {
                        int s0 = cg5_nz[iP].r;
                        int s1 = cg5_nz[iP].c;

                        for (int iN = 0; iN < n_cg5_nz; iN++) {
                          int s3 = cg5_nz[iN].r;
                          int s4 = cg5_nz[iN].c;

                          int s[NUM_SLOTS] = {s0, s1, s2, s3, s4, s5};

                          acc += cg5_vals[iP] * cg5_vals[iN] *
                            De[0][s[topo.cg5[0][0]]][s[topo.cg5[0][1]]] *
                            De[1][s[topo.cg5[1][0]]][s[topo.cg5[1][1]]] *
                            De[2][s[topo.cg5[2][0]]][s[topo.cg5[2][1]]];
                        }
                      }

                      result_local[s2][s5] += color_sign * acc;

                    } // s5
                  } // s2

                } // eN
              } // eP
            } // eB
          } // eA
        } // topology

        // Store result into LatticeSpinMatrix at this site
        SpinMatrix res_site = zero;
        for (int s2 = 0; s2 < Nd; s2++)
          for (int s5 = 0; s5 < Nd; s5++) {
            Complex val = cmplx(Real(result_local[s2][s5].real()),
                                Real(result_local[s2][s5].imag()));
            pokeSpin(res_site, val, s2, s5);
          }
        result.elem(site) = res_site.elem();

      } // site loop

      END_CODE();
      return result;

#else
      LatticeSpinMatrix a;
      a = zero;
      return a;
#endif
    }


    // =================================================================
    // Reference (original) implementation using lattice-wide QDP++ ops.
    // Kept for verification - produces identical results but is O(23M)
    // lattice-wide operations instead of a single site sweep.
    // =================================================================
    LatticeSpinMatrix dibaryonExchangeReference(
        const LatticePropagator& S_u,
        const LatticePropagator& S_d,
        const SpinMatrix& Cg5,
        const SpinMatrix& T_unpol)
    {
#if QDP_NC == 3
      START_CODE();

      // ---------------------------------------------------------
      // Step 1: Extract all color blocks as LatticeSpinMatrix
      // ---------------------------------------------------------
      LatticeSpinMatrix S_color[2][Nc][Nc];
      for (int a = 0; a < Nc; a++)
        for (int b = 0; b < Nc; b++) {
          S_color[0][a][b] = peekColor(S_u, a, b);
          S_color[1][a][b] = peekColor(S_d, a, b);
        }

      // ---------------------------------------------------------
      // Step 2a: Precompute Cg5 * transpose(S)
      // ---------------------------------------------------------
      LatticeSpinMatrix CgST[2][Nc][Nc];
      for (int f = 0; f < 2; f++)
        for (int a = 0; a < Nc; a++)
          for (int b = 0; b < Nc; b++)
            CgST[f][a][b] = Cg5 * transposeSpin(S_color[f][a][b]);

      // ---------------------------------------------------------
      // Step 2b: Precompute T_unpol-projected propagator sums
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
      // Step 3: Precompute nonzero Cg5 entries
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
      // Step 4: Accumulate results
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

        int slot_eps_pos[NUM_SLOTS];
        int slot_eps_grp[NUM_SLOTS];
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

        for (int eA = 0; eA < NUM_EPS_ENTRIES; eA++)
        {
          int cA[3] = {eps_table[eA].i, eps_table[eA].j, eps_table[eA].k};
          int signA = eps_table[eA].sign;

          for (int eB = 0; eB < NUM_EPS_ENTRIES; eB++)
          {
            int cB[3] = {eps_table[eB].i, eps_table[eB].j, eps_table[eB].k};
            int signB = eps_table[eB].sign;

            int c_src[NUM_SLOTS];
            for (int s = 0; s < NUM_SLOTS; s++) {
              c_src[s] = (slot_eps_grp[s] == 0) ?
                cA[slot_eps_pos[s]] : cB[slot_eps_pos[s]];
            }

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

                LatticeComplex De[3][Nd][Nd];

                for (int k = 0; k < 3; k++) {
                  int pk = topo.cg5[k][0];
                  int qk = topo.cg5[k][1];

                  if (k == topo.open_pair) {
                    for (int r = 0; r < Nd; r++)
                      for (int c = 0; c < Nd; c++)
                        De[k][r][c] =
                          Tv[slot_flavor[pk]][a_snk[pk]][c_src[pk]][r] *
                          Tv[slot_flavor[qk]][a_snk[qk]][c_src[qk]][c];
                  } else {
                    LatticeSpinMatrix D =
                      S_color[slot_flavor[pk]][a_snk[pk]][c_src[pk]]
                      * CgST[slot_flavor[qk]][a_snk[qk]][c_src[qk]];
                    for (int r = 0; r < Nd; r++)
                      for (int c = 0; c < Nd; c++)
                        De[k][r][c] = peekSpin(D, r, c);
                  }
                }

                for (int s2 = 0; s2 < Nd; s2++) {
                  for (int s5 = 0; s5 < Nd; s5++) {

                    LatticeComplex acc;
                    acc = zero;

                    for (int iP = 0; iP < n_cg5_nz; iP++) {
                      int s0 = cg5_nz[iP].r;
                      int s1 = cg5_nz[iP].c;

                      for (int iN = 0; iN < n_cg5_nz; iN++) {
                        int s3 = cg5_nz[iN].r;
                        int s4 = cg5_nz[iN].c;

                        int s[NUM_SLOTS] = {s0, s1, s2, s3, s4, s5};

                        acc += cg5_nz[iP].val * cg5_nz[iN].val *
                          De[0][s[topo.cg5[0][0]]][s[topo.cg5[0][1]]] *
                          De[1][s[topo.cg5[1][0]]][s[topo.cg5[1][1]]] *
                          De[2][s[topo.cg5[2][0]]][s[topo.cg5[2][1]]];
                      }
                    }

                    result_elem[s2][s5] +=
                      cmplx(Real(color_sign), Real(0)) * acc;

                  }
                }

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
