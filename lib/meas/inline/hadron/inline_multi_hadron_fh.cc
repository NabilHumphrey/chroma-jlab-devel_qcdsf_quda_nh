/**
 * @file inline_multi_hadron_fh.cc
 * @brief Inline measurement for multi-hadron correlators with Feynman-Hellmann support
 *
 * Computes multi-hadron correlators (dibaryons, dimesons) from pre-computed
 * propagators using product ansatz: C_AB(t) = sum_p C_A(p,t) * C_B(-p,t)
 */

#include "inline_multi_hadron_fh.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"

#include <fstream>

namespace Chroma
{
  namespace InlineMultiHadronFHEnv
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                                              const std::string& path)
      {
        return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    // Measurement name
    const std::string measurement_name = "MULTI_HADRON_FH";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (! registered)
      {
        success &= TheInlineMeasurementFactory::Instance().registerObject(measurement_name, createMeasurement);
        registered = true;
      }
      return success;
    }


    //------------------------------------------------------------------
    // BoostedMomentum_t I/O
    //------------------------------------------------------------------
    void read(XMLReader& xml, const std::string& path, Params::BoostedMomentum_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "total_momentum", param.total_momentum);
      read(paramtop, "fiducial_mom1", param.fiducial_mom1);
      read(paramtop, "fiducial_mom2", param.fiducial_mom2);
      read(paramtop, "irrep", param.irrep);
    }

    void write(XMLWriter& xml, const std::string& path, const Params::BoostedMomentum_t& param)
    {
      push(xml, path);
      write(xml, "total_momentum", param.total_momentum);
      write(xml, "fiducial_mom1", param.fiducial_mom1);
      write(xml, "fiducial_mom2", param.fiducial_mom2);
      write(xml, "irrep", param.irrep);
      pop(xml);
    }


    //------------------------------------------------------------------
    // MultiHadronState_t I/O
    //------------------------------------------------------------------
    void read(XMLReader& xml, const std::string& path, Params::MultiHadronState_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "name", param.name);
      read(paramtop, "type", param.type);
      read(paramtop, "hadron1", param.hadron1);
      read(paramtop, "hadron2", param.hadron2);

      // Read boosted momentum parameters (optional)
      if (paramtop.count("is_boosted") == 1)
      {
        read(paramtop, "is_boosted", param.is_boosted);
        if (param.is_boosted)
        {
          read(paramtop, "boost", param.boost);
        }
      }
      else
      {
        param.is_boosted = false;
      }
    }

    void write(XMLWriter& xml, const std::string& path, const Params::MultiHadronState_t& param)
    {
      push(xml, path);
      write(xml, "name", param.name);
      write(xml, "type", param.type);
      write(xml, "hadron1", param.hadron1);
      write(xml, "hadron2", param.hadron2);
      write(xml, "is_boosted", param.is_boosted);
      if (param.is_boosted)
      {
        write(xml, "boost", param.boost);
      }
      pop(xml);
    }


    //------------------------------------------------------------------
    // FlavorCombination_t I/O
    //------------------------------------------------------------------
    void read(XMLReader& xml, const std::string& path, Params::FlavorCombination_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "up_is_perturbed", param.up_is_perturbed);
      read(paramtop, "down_is_perturbed", param.down_is_perturbed);
      read(paramtop, "label", param.label);
    }

    void write(XMLWriter& xml, const std::string& path, const Params::FlavorCombination_t& param)
    {
      push(xml, path);
      write(xml, "up_is_perturbed", param.up_is_perturbed);
      write(xml, "down_is_perturbed", param.down_is_perturbed);
      write(xml, "label", param.label);
      pop(xml);
    }


    //------------------------------------------------------------------
    // Params
    //------------------------------------------------------------------
    Params::Params()
    {
      frequency = 0;
      param.mom2_max = 4;
      param.avg_equiv_mom = true;
      param.t0 = 0;
    }

    Params::Params(XMLReader& xml_in, const std::string& path)
    {
      try
      {
        XMLReader paramtop(xml_in, path);

        if (paramtop.count("Frequency") == 1)
          read(paramtop, "Frequency", frequency);
        else
          frequency = 1;

        // Read Param section
        {
          XMLReader ptop(paramtop, "Param");

          read(ptop, "mom2_max", param.mom2_max);

          if (ptop.count("avg_equiv_mom") == 1)
            read(ptop, "avg_equiv_mom", param.avg_equiv_mom);
          else
            param.avg_equiv_mom = true;

          // Read source timeslice
          if (ptop.count("t0") == 1)
            read(ptop, "t0", param.t0);
          else
            param.t0 = 0;

          // Read multi-hadron states list
          if (ptop.count("states") == 1)
          {
            XMLReader stattop(ptop, "states");
            int num_states = stattop.count("elem");
            param.states.resize(num_states);
            for (int i = 0; i < num_states; ++i)
            {
              read(stattop, "elem[" + std::to_string(i+1) + "]", param.states[i]);
            }
          }
          else
          {
            // Default: deuteron (proton-neutron, spin-1)
            param.states.resize(1);
            param.states[0].name = "deuteron";
            param.states[0].type = "DIBARYON";
            param.states[0].hadron1 = "PROTON";
            param.states[0].hadron2 = "NEUTRON";
          }

          // Read contractions list
          if (ptop.count("contractions") == 1)
          {
            XMLReader contop(ptop, "contractions");
            int num_contractions = contop.count("elem");
            param.contractions.resize(num_contractions);
            for (int i = 0; i < num_contractions; ++i)
            {
              read(contop, "elem[" + std::to_string(i+1) + "]", param.contractions[i]);
            }
          }
          else
          {
            // Default: unperturbed
            param.contractions.resize(1);
            param.contractions[0].up_is_perturbed = false;
            param.contractions[0].down_is_perturbed = false;
            param.contractions[0].label = "unperturbed";
          }
        }

        // Read NamedObject section
        {
          XMLReader ntop(paramtop, "NamedObject");
          read(ntop, "gauge_id", named_obj.gauge_id);
          read(ntop, "perturbed_prop_id", named_obj.perturbed_prop_id);
          read(ntop, "unperturbed_prop_id", named_obj.unperturbed_prop_id);
          read(ntop, "output_file", named_obj.output_file);
        }
      }
      catch(const std::string& e)
      {
        QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
        QDP_abort(1);
      }
    }

    void Params::writeXML(XMLWriter& xml_out, const std::string& path) const
    {
      push(xml_out, path);

      write(xml_out, "Frequency", frequency);

      push(xml_out, "Param");
      write(xml_out, "mom2_max", param.mom2_max);
      write(xml_out, "avg_equiv_mom", param.avg_equiv_mom);
      write(xml_out, "t0", param.t0);

      push(xml_out, "states");
      for (int i = 0; i < param.states.size(); ++i)
      {
        write(xml_out, "elem", param.states[i]);
      }
      pop(xml_out);

      push(xml_out, "contractions");
      for (int i = 0; i < param.contractions.size(); ++i)
      {
        write(xml_out, "elem", param.contractions[i]);
      }
      pop(xml_out);

      pop(xml_out);

      push(xml_out, "NamedObject");
      write(xml_out, "gauge_id", named_obj.gauge_id);
      write(xml_out, "perturbed_prop_id", named_obj.perturbed_prop_id);
      write(xml_out, "unperturbed_prop_id", named_obj.unperturbed_prop_id);
      write(xml_out, "output_file", named_obj.output_file);
      pop(xml_out);

      pop(xml_out);
    }


    //------------------------------------------------------------------
    // Spin-resolved baryon correlator computation
    //------------------------------------------------------------------

    //! Compute spin-resolved proton correlators
    /*!
     * Returns separate spin-up and spin-down correlators for proper
     * spin projection in multi-hadron states.
     *
     * After positive parity projection, the Dirac structure is:
     *   | G_↑↑  G_↑↓  0  0 |
     *   | G_↓↑  G_↓↓  0  0 |
     *   |  0     0    0  0 |
     *   |  0     0    0  0 |
     *
     * At zero momentum, off-diagonal elements vanish: G_↑ = (0,0), G_↓ = (1,1)
     *
     * \param up_prop      Propagator for up quarks (Read)
     * \param down_prop    Propagator for down quarks (Read)
     * \param phases       Momentum phases (Read)
     * \param t0           Source timeslice for t_eff calculation (Read)
     * \param corr_spin_up   Spin-up correlator output (Write)
     * \param corr_spin_down Spin-down correlator output (Write)
     */
    void protonCorrelatorSpinResolved(
      const LatticePropagator& up_prop,
      const LatticePropagator& down_prop,
      const SftMom& phases,
      int t0,
      multi2d<DComplex>& corr_spin_up,
      multi2d<DComplex>& corr_spin_down)
    {
      START_CODE();

      int num_mom = phases.numMom();
      int Lt = phases.numSubsets();

      corr_spin_up.resize(num_mom, Lt);
      corr_spin_down.resize(num_mom, Lt);

      SpinMatrix g_one = 1.0;
      SpinMatrix Cg5 = Gamma(10) * (Gamma(15) * g_one);

      LatticePropagator diquark = quarkContract13(up_prop * Cg5,
                                                   Gamma(15) * down_prop);

      // Positive parity projection
      SpinMatrix parity_proj = 0.5 * (g_one + Gamma(8) * g_one);

      // Get the full spin matrix (color traced, with third quark contracted)
      LatticeSpinMatrix spin_matrix = parity_proj * traceColor(diquark * up_prop);

      // Extract spin-up component: Dirac indices (0,0)
      LatticeComplex spin_up_field = peekSpin(spin_matrix, 0, 0);

      // Extract spin-down component: Dirac indices (1,1)
      LatticeComplex spin_down_field = peekSpin(spin_matrix, 1, 1);

      // Momentum projection for spin-up
      multi2d<DComplex> hsum_up = phases.sft(spin_up_field);

      // Momentum projection for spin-down
      multi2d<DComplex> hsum_down = phases.sft(spin_down_field);

      // Apply t_eff shift
      for (int p = 0; p < num_mom; ++p)
        for (int t = 0; t < Lt; ++t)
        {
          int t_eff = (t - t0 + Lt) % Lt;
          corr_spin_up[p][t_eff] = hsum_up[p][t];
          corr_spin_down[p][t_eff] = hsum_down[p][t];
        }

      END_CODE();
    }

    //! Compute spin-resolved neutron correlators (u <-> d swap from proton)
    void neutronCorrelatorSpinResolved(
      const LatticePropagator& up_prop,
      const LatticePropagator& down_prop,
      const SftMom& phases,
      int t0,
      multi2d<DComplex>& corr_spin_up,
      multi2d<DComplex>& corr_spin_down)
    {
      protonCorrelatorSpinResolved(down_prop, up_prop, phases, t0, corr_spin_up, corr_spin_down);
    }

    //! Compute spin-traced proton correlator (for backward compatibility and output)
    void protonCorrelator(
      const LatticePropagator& up_prop,
      const LatticePropagator& down_prop,
      const SftMom& phases,
      int t0,
      multi2d<DComplex>& corr)
    {
      multi2d<DComplex> corr_up, corr_down;
      protonCorrelatorSpinResolved(up_prop, down_prop, phases, t0, corr_up, corr_down);

      int num_mom = phases.numMom();
      int Lt = phases.numSubsets();
      corr.resize(num_mom, Lt);

      for (int p = 0; p < num_mom; ++p)
        for (int t = 0; t < Lt; ++t)
          corr[p][t] = corr_up[p][t] + corr_down[p][t];
    }

    //! Compute spin-traced neutron correlator
    void neutronCorrelator(
      const LatticePropagator& up_prop,
      const LatticePropagator& down_prop,
      const SftMom& phases,
      int t0,
      multi2d<DComplex>& corr)
    {
      protonCorrelator(down_prop, up_prop, phases, t0, corr);
    }


    //------------------------------------------------------------------
    // Vector meson (rho) correlator computation
    //------------------------------------------------------------------

    //! Compute rho meson correlator for a specific polarization
    /*!
     * Rho meson correlator: C(t) = -Tr[Gamma * S_q * Gamma * gamma5 * S_qbar^dag * gamma5]
     * where the minus sign accounts for the fermion loop.
     *
     * Linear polarizations in Chroma convention:
     *   gamma = 1 (gamma_1) -> x-polarized (rho_x)
     *   gamma = 2 (gamma_2) -> y-polarized (rho_y)
     *   gamma = 4 (gamma_3) -> z-polarized (rho_z)
     *
     * For isovector rho (u dbar or d ubar), quark_prop and antiquark_prop
     * are the same propagator but represent different flavors.
     *
     * \param quark_prop     Quark propagator (Read)
     * \param antiquark_prop Antiquark propagator (Read)
     * \param gamma          Gamma matrix index for polarization (Read)
     * \param phases         Momentum phases (Read)
     * \param t0             Source timeslice for t_eff calculation (Read)
     * \param corr           Correlator output [momentum][time] (Write)
     */
    void rhoMesonCorrelator(
      const LatticePropagator& quark_prop,
      const LatticePropagator& antiquark_prop,
      int gamma,
      const SftMom& phases,
      int t0,
      multi2d<DComplex>& corr)
    {
      START_CODE();

      int num_mom = phases.numMom();
      int Lt = phases.numSubsets();

      corr.resize(num_mom, Lt);

      // Construct antiquark from propagator: Sbar = gamma5 * S^dag * gamma5
      LatticePropagator anti_quark = Gamma(15) * adj(antiquark_prop) * Gamma(15);

      // Meson correlator: -Tr[Gamma * S * Gamma * Sbar]
      // The minus sign is for the fermion loop
      LatticeComplex corr_fn = -trace(Gamma(gamma) * quark_prop * Gamma(gamma) * anti_quark);

      // Momentum projection
      multi2d<DComplex> hsum = phases.sft(corr_fn);

      // Apply t_eff shift
      for (int p = 0; p < num_mom; ++p)
        for (int t = 0; t < Lt; ++t)
        {
          int t_eff = (t - t0 + Lt) % Lt;
          corr[p][t_eff] = hsum[p][t];
        }

      END_CODE();
    }

    //! Compute rho meson correlators for all three linear polarizations
    /*!
     * Computes x, y, z polarized rho correlators simultaneously.
     *
     * \param quark_prop     Quark propagator (Read)
     * \param antiquark_prop Antiquark propagator (Read)
     * \param phases         Momentum phases (Read)
     * \param t0             Source timeslice for t_eff calculation (Read)
     * \param corr_x         x-polarized correlator output [momentum][time] (Write)
     * \param corr_y         y-polarized correlator output [momentum][time] (Write)
     * \param corr_z         z-polarized correlator output [momentum][time] (Write)
     */
    void rhoMesonCorrelatorAllPol(
      const LatticePropagator& quark_prop,
      const LatticePropagator& antiquark_prop,
      const SftMom& phases,
      int t0,
      multi2d<DComplex>& corr_x,
      multi2d<DComplex>& corr_y,
      multi2d<DComplex>& corr_z)
    {
      // gamma = 1 -> x, gamma = 2 -> y, gamma = 4 -> z (Chroma convention)
      rhoMesonCorrelator(quark_prop, antiquark_prop, 1, phases, t0, corr_x);
      rhoMesonCorrelator(quark_prop, antiquark_prop, 2, phases, t0, corr_y);
      rhoMesonCorrelator(quark_prop, antiquark_prop, 4, phases, t0, corr_z);
    }


    //------------------------------------------------------------------
    // Multi-hadron correlator computation with spin projection
    //------------------------------------------------------------------

    //! Find index of -p momentum in the phases list
    int findNegativeMomentum(const SftMom& phases, int ip)
    {
      multi1d<int> mom_p = phases.numToMom(ip);
      int num_mom = phases.numMom();

      for (int jp = 0; jp < num_mom; ++jp)
      {
        multi1d<int> mom_q = phases.numToMom(jp);

        bool match = true;
        for (int i = 0; i < Nd-1; ++i)
        {
          if (mom_q[i] != -mom_p[i])
          {
            match = false;
            break;
          }
        }
        if (match)
          return jp;
      }
      return -1;
    }

    //------------------------------------------------------------------
    // Utility functions for boosted (nonzero total momentum) operators
    //------------------------------------------------------------------

    //! Find SftMom index for a specific 3-momentum vector
    /*!
     * Searches the SftMom momentum list for a momentum matching the given 3-vector.
     * Returns the index, or -1 if not found (momentum exceeds mom2_max).
     */
    int findMomentumIndex(const SftMom& phases, const multi1d<int>& mom)
    {
      int num_mom = phases.numMom();
      for (int ip = 0; ip < num_mom; ++ip)
      {
        multi1d<int> mom_p = phases.numToMom(ip);
        bool match = true;
        for (int i = 0; i < Nd-1; ++i)
        {
          if (mom_p[i] != mom[i])
          {
            match = false;
            break;
          }
        }
        if (match)
          return ip;
      }
      return -1;
    }

    //! Apply C4x coset rotation to a 3-momentum vector
    /*!
     * For P = (1,0,0), the C_{4v} coset representatives are:
     *   coset 0: e          (identity)
     *   coset 1: C_{4x}     (90° about x: y -> -z, z -> y)
     *   coset 2: C_{2x}     (180° about x: y -> -y, z -> -z)
     *   coset 3: C_{4x}^{-1} (270° about x: y -> z, z -> -y)
     */
    void applyC4xCosetRotation(const multi1d<int>& mom, int coset_idx, multi1d<int>& result)
    {
      result.resize(3);
      switch(coset_idx)
      {
        case 0: // Identity
          result[0] = mom[0];
          result[1] = mom[1];
          result[2] = mom[2];
          break;
        case 1: // C_{4x}: (x,y,z) -> (x,-z,y)
          result[0] = mom[0];
          result[1] = -mom[2];
          result[2] = mom[1];
          break;
        case 2: // C_{2x}: (x,y,z) -> (x,-y,-z)
          result[0] = mom[0];
          result[1] = -mom[1];
          result[2] = -mom[2];
          break;
        case 3: // C_{4x}^{-1}: (x,y,z) -> (x,z,-y)
          result[0] = mom[0];
          result[1] = mom[2];
          result[2] = -mom[1];
          break;
        default:
          QDPIO::cerr << "applyC4xCosetRotation: invalid coset index " << coset_idx << std::endl;
          QDP_abort(1);
      }
    }

    //! Build the 4-element C4x orbit from fiducial momenta
    /*!
     * Constructs all 4 momentum configurations in the orbit by applying
     * coset representatives {e, C4x, C2x, C4x^{-1}} to the fiducial pair.
     *
     * \param fid_mom1    Fiducial momentum 1 (Read)
     * \param fid_mom2    Fiducial momentum 2 (Read)
     * \param orbit_mom1  Output: 4 rotated momentum-1 vectors (Write)
     * \param orbit_mom2  Output: 4 rotated momentum-2 vectors (Write)
     */
    void buildC4xOrbit(const multi1d<int>& fid_mom1, const multi1d<int>& fid_mom2,
                       multi1d<multi1d<int>>& orbit_mom1, multi1d<multi1d<int>>& orbit_mom2)
    {
      orbit_mom1.resize(4);
      orbit_mom2.resize(4);

      for (int m = 0; m < 4; ++m)
      {
        applyC4xCosetRotation(fid_mom1, m, orbit_mom1[m]);
        applyC4xCosetRotation(fid_mom2, m, orbit_mom2[m]);
      }
    }


    //------------------------------------------------------------------
    // Boosted dibaryon correlator (nonzero total momentum)
    //------------------------------------------------------------------

    //! Compute boosted dibaryon correlators with C4v irrep projection
    /*!
     * For spin-1 dibaryon at P = (1,0,0), projects onto A1 (longitudinal)
     * and E (transverse) irreps of C_{4v}.
     *
     * At each orbit configuration m, computes:
     *   C^SS(m) = 0.5*(G_up1(n1_m)*G_up2(n2_m) + G_down1(n1_m)*G_down2(n2_m))
     *   C^mix(m) = 0.5*(G_up1(n1_m)*G_down2(n2_m) + G_down1(n1_m)*G_up2(n2_m))
     *
     * Then accumulates with rotation weights:
     *   A1:  += C^SS(m)    [uniform sum of x-polarization]
     *   Ey:  += Ryy(m)*C^SS(m) + Ryz(m)*C^mix(m)
     *   Ez:  += Rzy(m)*C^SS(m) + Rzz(m)*C^mix(m)
     *
     * Rotation coefficients for coset reps {e, C4x, C2x, C4x^{-1}}:
     *   Ryy = {1, 0, -1, 0}, Ryz = {0, -1, 0, 1}
     *   Rzy = {0, 1, 0, -1}, Rzz = {1, 0, -1, 0}
     */
    void dibaryonCorrelatorBoostedCartesian(
      const multi2d<DComplex>& h1_up,
      const multi2d<DComplex>& h1_down,
      const multi2d<DComplex>& h2_up,
      const multi2d<DComplex>& h2_down,
      const SftMom& phases,
      const multi1d<int>& fid_mom1,
      const multi1d<int>& fid_mom2,
      multi1d<DComplex>& corr_A1,
      multi1d<DComplex>& corr_Ey,
      multi1d<DComplex>& corr_Ez)
    {
      START_CODE();

      int Lt = h1_up.size2();

      corr_A1.resize(Lt);
      corr_Ey.resize(Lt);
      corr_Ez.resize(Lt);
      corr_A1 = zero;
      corr_Ey = zero;
      corr_Ez = zero;

      // Build the 4-element orbit
      multi1d<multi1d<int>> orbit_mom1, orbit_mom2;
      buildC4xOrbit(fid_mom1, fid_mom2, orbit_mom1, orbit_mom2);

      // Rotation coefficients for yz-block of coset reps
      // {e, C4x, C2x, C4x^{-1}}
      int Ryy[4] = { 1,  0, -1,  0};
      int Ryz[4] = { 0, -1,  0,  1};
      int Rzy[4] = { 0,  1,  0, -1};
      int Rzz[4] = { 1,  0, -1,  0};

      for (int m = 0; m < 4; ++m)
      {
        // Find momentum indices for this orbit configuration
        int ip1 = findMomentumIndex(phases, orbit_mom1[m]);
        int ip2 = findMomentumIndex(phases, orbit_mom2[m]);

        if (ip1 < 0 || ip2 < 0)
        {
          QDPIO::cerr << "dibaryonCorrelatorBoostedCartesian: momentum not found for orbit element "
                      << m << ". Check that mom2_max >= 2." << std::endl;
          QDPIO::cerr << "  n1 = (" << orbit_mom1[m][0] << "," << orbit_mom1[m][1] << "," << orbit_mom1[m][2] << ")"
                      << "  n2 = (" << orbit_mom2[m][0] << "," << orbit_mom2[m][1] << "," << orbit_mom2[m][2] << ")" << std::endl;
          QDP_abort(1);
        }

        for (int t = 0; t < Lt; ++t)
        {
          // C^SS(m) = 0.5*(up1*up2 + down1*down2) — corresponds to C_x = C_y
          DComplex CSS = 0.5 * (h1_up[ip1][t] * h2_up[ip2][t] +
                                h1_down[ip1][t] * h2_down[ip2][t]);

          // C^mix(m) = 0.5*(up1*down2 + down1*up2) — corresponds to C_z
          DComplex Cmix = 0.5 * (h1_up[ip1][t] * h2_down[ip2][t] +
                                 h1_down[ip1][t] * h2_up[ip2][t]);

          // A1 (longitudinal): uniform sum of x-polarization = C^SS
          corr_A1[t] += CSS;

          // E_y: Ryy*C^SS + Ryz*C^mix
          corr_Ey[t] += Real(Ryy[m]) * CSS + Real(Ryz[m]) * Cmix;

          // E_z: Rzy*C^SS + Rzz*C^mix
          corr_Ez[t] += Real(Rzy[m]) * CSS + Real(Rzz[m]) * Cmix;
        }
      }

      END_CODE();
    }


    //------------------------------------------------------------------
    // Boosted di-rho correlator (nonzero total momentum)
    //------------------------------------------------------------------

    //! Compute boosted di-rho correlators with C4v irrep projection
    /*!
     * Same orbit structure as boosted dibaryon, but uses Levi-Civita coupling
     * for spin-1 from 1x1.
     *
     * At each orbit config m, first computes S=1 Cartesian components:
     *   C_x^{S=1}(m) = 0.5*[C_y^rho1(n1_m)*C_z^rho2(n2_m) + C_z^rho1(n1_m)*C_y^rho2(n2_m)]
     *   C_y^{S=1}(m) = 0.5*[C_z^rho1(n1_m)*C_x^rho2(n2_m) + C_x^rho1(n1_m)*C_z^rho2(n2_m)]
     *   C_z^{S=1}(m) = 0.5*[C_x^rho1(n1_m)*C_y^rho2(n2_m) + C_y^rho1(n1_m)*C_x^rho2(n2_m)]
     *
     * Then projects:
     *   A1:  += C_x^{S=1}(m)
     *   Ey:  += Ryy(m)*C_y^{S=1}(m) + Ryz(m)*C_z^{S=1}(m)
     *   Ez:  += Rzy(m)*C_y^{S=1}(m) + Rzz(m)*C_z^{S=1}(m)
     */
    void diRhoCorrelatorBoostedCartesian(
      const multi2d<DComplex>& rho1_x,
      const multi2d<DComplex>& rho1_y,
      const multi2d<DComplex>& rho1_z,
      const multi2d<DComplex>& rho2_x,
      const multi2d<DComplex>& rho2_y,
      const multi2d<DComplex>& rho2_z,
      const SftMom& phases,
      const multi1d<int>& fid_mom1,
      const multi1d<int>& fid_mom2,
      multi1d<DComplex>& corr_A1,
      multi1d<DComplex>& corr_Ey,
      multi1d<DComplex>& corr_Ez)
    {
      START_CODE();

      int Lt = rho1_x.size2();

      corr_A1.resize(Lt);
      corr_Ey.resize(Lt);
      corr_Ez.resize(Lt);
      corr_A1 = zero;
      corr_Ey = zero;
      corr_Ez = zero;

      // Build the 4-element orbit
      multi1d<multi1d<int>> orbit_mom1, orbit_mom2;
      buildC4xOrbit(fid_mom1, fid_mom2, orbit_mom1, orbit_mom2);

      // Rotation coefficients for yz-block of coset reps
      int Ryy[4] = { 1,  0, -1,  0};
      int Ryz[4] = { 0, -1,  0,  1};
      int Rzy[4] = { 0,  1,  0, -1};
      int Rzz[4] = { 1,  0, -1,  0};

      for (int m = 0; m < 4; ++m)
      {
        int ip1 = findMomentumIndex(phases, orbit_mom1[m]);
        int ip2 = findMomentumIndex(phases, orbit_mom2[m]);

        if (ip1 < 0 || ip2 < 0)
        {
          QDPIO::cerr << "diRhoCorrelatorBoostedCartesian: momentum not found for orbit element "
                      << m << ". Check that mom2_max >= 2." << std::endl;
          QDP_abort(1);
        }

        for (int t = 0; t < Lt; ++t)
        {
          // Levi-Civita coupled S=1 components at this orbit config
          // C_x^{S=1} = 0.5*(rho1_y * rho2_z + rho1_z * rho2_y)
          DComplex Cx_S1 = 0.5 * (rho1_y[ip1][t] * rho2_z[ip2][t] +
                                   rho1_z[ip1][t] * rho2_y[ip2][t]);

          // C_y^{S=1} = 0.5*(rho1_z * rho2_x + rho1_x * rho2_z)
          DComplex Cy_S1 = 0.5 * (rho1_z[ip1][t] * rho2_x[ip2][t] +
                                   rho1_x[ip1][t] * rho2_z[ip2][t]);

          // C_z^{S=1} = 0.5*(rho1_x * rho2_y + rho1_y * rho2_x)
          DComplex Cz_S1 = 0.5 * (rho1_x[ip1][t] * rho2_y[ip2][t] +
                                   rho1_y[ip1][t] * rho2_x[ip2][t]);

          // A1 (longitudinal): sum of x-component
          corr_A1[t] += Cx_S1;

          // E_y: Ryy*Cy + Ryz*Cz
          corr_Ey[t] += Real(Ryy[m]) * Cy_S1 + Real(Ryz[m]) * Cz_S1;

          // E_z: Rzy*Cy + Rzz*Cz
          corr_Ez[t] += Real(Rzy[m]) * Cy_S1 + Real(Rzz[m]) * Cz_S1;
        }
      }

      END_CODE();
    }


    //! Compute dibaryon correlators in Cartesian polarization basis
    /*!
     * For spin-1 (triplet) dibaryon states in the Cartesian basis:
     *   C_x(t) = 0.5 * [G↑₁G↑₂ + G↓₁G↓₂]
     *   C_y(t) = 0.5 * [G↑₁G↑₂ + G↓₁G↓₂]    (identical to C_x for product ansatz)
     *   C_z(t) = 0.5 * [G↑₁G↓₂ + G↓₁G↑₂]
     *
     * These are needed for studying polarized structure functions:
     *   - Tensor polarization: b1 ~ C_z - 0.5*(C_x + C_y)
     *
     * The sum C_x + C_y + C_z gives the unpolarized spin-1 correlator.
     * Note: C_x = C_y identically (bit-for-bit) but both are output for uniformity.
     *
     * \param h1_up     Hadron 1 spin-up correlator [momentum][time]
     * \param h1_down   Hadron 1 spin-down correlator [momentum][time]
     * \param h2_up     Hadron 2 spin-up correlator [momentum][time]
     * \param h2_down   Hadron 2 spin-down correlator [momentum][time]
     * \param phases    Momentum phases
     * \param corr_x    x-polarization correlator output [time]
     * \param corr_y    y-polarization correlator output [time]
     * \param corr_z    z-polarization correlator output [time]
     */
    void dibaryonCorrelatorCartesian(
      const multi2d<DComplex>& h1_up,
      const multi2d<DComplex>& h1_down,
      const multi2d<DComplex>& h2_up,
      const multi2d<DComplex>& h2_down,
      const SftMom& phases,
      multi1d<DComplex>& corr_x,
      multi1d<DComplex>& corr_y,
      multi1d<DComplex>& corr_z)
    {
      START_CODE();

      int num_mom = phases.numMom();
      int Lt = h1_up.size2();

      corr_x.resize(Lt);
      corr_y.resize(Lt);
      corr_z.resize(Lt);
      corr_x = zero;
      corr_y = zero;
      corr_z = zero;

      // For each momentum p, find -p and accumulate spin products
      for (int ip = 0; ip < num_mom; ++ip)
      {
        int ip_neg = findNegativeMomentum(phases, ip);
        if (ip_neg < 0) continue;

        for (int t = 0; t < Lt; ++t)
        {
          // C_x = C_y = 0.5 * (G↑₁G↑₂ + G↓₁G↓₂)
          DComplex upup   = h1_up[ip][t]   * h2_up[ip_neg][t];
          DComplex downdown = h1_down[ip][t] * h2_down[ip_neg][t];
          corr_x[t] += 0.5 * (upup + downdown);

          // C_z = 0.5 * (G↑₁G↓₂ + G↓₁G↑₂)
          DComplex updown = h1_up[ip][t]   * h2_down[ip_neg][t];
          DComplex downup = h1_down[ip][t] * h2_up[ip_neg][t];
          corr_z[t] += 0.5 * (updown + downup);
        }
      }

      // Normalize (factor of 1/2 for double counting in momentum sum)
      // and set C_y = C_x (identical for product ansatz)
      for (int t = 0; t < Lt; ++t)
      {
        corr_x[t] *= 0.5;
        corr_z[t] *= 0.5;
        corr_y[t] = corr_x[t];
      }

      END_CODE();
    }


    //! Compute spin-1 (triplet) dibaryon correlator as sum of Cartesian components
    /*!
     * Spin-1 (triplet, e.g., deuteron):
     *   C^{S=1}(t) = C_x + C_y + C_z
     *
     * \param corr_x    x-polarization correlator [time]
     * \param corr_y    y-polarization correlator [time]
     * \param corr_z    z-polarization correlator [time]
     * \param dibaryon_corr  Output correlator [time]
     */
    void dibaryonCorrelatorSpin1(
      const multi1d<DComplex>& corr_x,
      const multi1d<DComplex>& corr_y,
      const multi1d<DComplex>& corr_z,
      multi1d<DComplex>& dibaryon_corr)
    {
      int Lt = corr_x.size();
      dibaryon_corr.resize(Lt);

      for (int t = 0; t < Lt; ++t)
      {
        dibaryon_corr[t] = corr_x[t] + corr_y[t] + corr_z[t];
      }
    }


    //------------------------------------------------------------------
    // Vector-vector (di-rho) dimeson correlators with spin projection
    //------------------------------------------------------------------

    //! Compute di-rho correlators in Cartesian basis using Levi-Civita (cross product) coupling
    /*!
     * For two spin-1 particles coupling to total spin S=1 (antisymmetric),
     * the Cartesian basis uses the Levi-Civita tensor epsilon_{ijk}:
     *
     *   C_x(t) = 0.5 * sum_p [C_y^rho1(p) * C_z^rho2(-p) + C_z^rho1(p) * C_y^rho2(-p)]
     *   C_y(t) = 0.5 * sum_p [C_z^rho1(p) * C_x^rho2(-p) + C_x^rho1(p) * C_z^rho2(-p)]
     *   C_z(t) = 0.5 * sum_p [C_x^rho1(p) * C_y^rho2(-p) + C_y^rho1(p) * C_x^rho2(-p)]
     *
     * The 0.5 CG prefactor comes from |1/sqrt(2)|^2.
     * An additional 0.5 momentum normalization for double counting is applied.
     *
     * \param rho1_x    Rho 1 x-polarized correlator [momentum][time]
     * \param rho1_y    Rho 1 y-polarized correlator [momentum][time]
     * \param rho1_z    Rho 1 z-polarized correlator [momentum][time]
     * \param rho2_x    Rho 2 x-polarized correlator [momentum][time]
     * \param rho2_y    Rho 2 y-polarized correlator [momentum][time]
     * \param rho2_z    Rho 2 z-polarized correlator [momentum][time]
     * \param phases    Momentum phases
     * \param corr_x    x-polarization correlator output [time]
     * \param corr_y    y-polarization correlator output [time]
     * \param corr_z    z-polarization correlator output [time]
     */
    void diRhoCorrelatorCartesian(
      const multi2d<DComplex>& rho1_x,
      const multi2d<DComplex>& rho1_y,
      const multi2d<DComplex>& rho1_z,
      const multi2d<DComplex>& rho2_x,
      const multi2d<DComplex>& rho2_y,
      const multi2d<DComplex>& rho2_z,
      const SftMom& phases,
      multi1d<DComplex>& corr_x,
      multi1d<DComplex>& corr_y,
      multi1d<DComplex>& corr_z)
    {
      START_CODE();

      int num_mom = phases.numMom();
      int Lt = rho1_x.size2();

      corr_x.resize(Lt);
      corr_y.resize(Lt);
      corr_z.resize(Lt);
      corr_x = zero;
      corr_y = zero;
      corr_z = zero;

      // For each momentum p, find -p and accumulate Levi-Civita products
      for (int ip = 0; ip < num_mom; ++ip)
      {
        int ip_neg = findNegativeMomentum(phases, ip);
        if (ip_neg < 0) continue;

        for (int t = 0; t < Lt; ++t)
        {
          // C_x = 0.5 * [C_y^1 * C_z^2 + C_z^1 * C_y^2]
          corr_x[t] += 0.5 * (rho1_y[ip][t] * rho2_z[ip_neg][t] +
                               rho1_z[ip][t] * rho2_y[ip_neg][t]);

          // C_y = 0.5 * [C_z^1 * C_x^2 + C_x^1 * C_z^2]
          corr_y[t] += 0.5 * (rho1_z[ip][t] * rho2_x[ip_neg][t] +
                               rho1_x[ip][t] * rho2_z[ip_neg][t]);

          // C_z = 0.5 * [C_x^1 * C_y^2 + C_y^1 * C_x^2]
          corr_z[t] += 0.5 * (rho1_x[ip][t] * rho2_y[ip_neg][t] +
                               rho1_y[ip][t] * rho2_x[ip_neg][t]);
        }
      }

      // Normalize (factor of 1/2 for double counting in momentum sum)
      for (int t = 0; t < Lt; ++t)
      {
        corr_x[t] *= 0.5;
        corr_y[t] *= 0.5;
        corr_z[t] *= 0.5;
      }

      END_CODE();
    }

    //! Compute total spin-1 di-rho correlator (sum of Cartesian components)
    void diRhoCorrelatorSpin1(
      const multi2d<DComplex>& rho1_x,
      const multi2d<DComplex>& rho1_y,
      const multi2d<DComplex>& rho1_z,
      const multi2d<DComplex>& rho2_x,
      const multi2d<DComplex>& rho2_y,
      const multi2d<DComplex>& rho2_z,
      const SftMom& phases,
      multi1d<DComplex>& dirho_corr)
    {
      multi1d<DComplex> corr_x, corr_y, corr_z;

      diRhoCorrelatorCartesian(rho1_x, rho1_y, rho1_z,
                               rho2_x, rho2_y, rho2_z,
                               phases, corr_x, corr_y, corr_z);

      int Lt = corr_x.size();
      dirho_corr.resize(Lt);

      for (int t = 0; t < Lt; ++t)
      {
        dirho_corr[t] = corr_x[t] + corr_y[t] + corr_z[t];
      }
    }


    //------------------------------------------------------------------
    // Main measurement
    //------------------------------------------------------------------
    void InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      push(xml_out, "MultiHadronFH");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << "MULTI_HADRON_FH: Multi-hadron correlator measurement with FH support" << std::endl;

      // Grab gauge configuration (for plaquette measurement)
      multi1d<LatticeColorMatrix> u =
        TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.gauge_id);

      // Measure plaquette
      MesPlq(xml_out, "Observables", u);

      // Load propagators
      QDPIO::cout << "Loading perturbed propagator: " << params.named_obj.perturbed_prop_id << std::endl;
      LatticePropagator perturbed_prop;
      try
      {
        perturbed_prop = TheNamedObjMap::Instance().getData<LatticePropagator>(
          params.named_obj.perturbed_prop_id);
      }
      catch (std::bad_cast)
      {
        QDPIO::cerr << "MULTI_HADRON_FH: perturbed_prop_id cast failed" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "MULTI_HADRON_FH: perturbed_prop_id not found: " << e << std::endl;
        QDP_abort(1);
      }

      QDPIO::cout << "Loading unperturbed propagator: " << params.named_obj.unperturbed_prop_id << std::endl;
      LatticePropagator unperturbed_prop;
      try
      {
        unperturbed_prop = TheNamedObjMap::Instance().getData<LatticePropagator>(
          params.named_obj.unperturbed_prop_id);
      }
      catch (std::bad_cast)
      {
        QDPIO::cerr << "MULTI_HADRON_FH: unperturbed_prop_id cast failed" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "MULTI_HADRON_FH: unperturbed_prop_id not found: " << e << std::endl;
        QDP_abort(1);
      }

      // Setup momentum projection
      SftMom phases(params.param.mom2_max, false, Nd-1);
      int num_mom = phases.numMom();
      int Lt = Layout::lattSize()[Nd-1];

      QDPIO::cout << "Number of momenta: " << num_mom << std::endl;
      QDPIO::cout << "Lattice time extent: " << Lt << std::endl;
      QDPIO::cout << "Source timeslice t0: " << params.param.t0 << std::endl;
      QDPIO::cout << "Number of multi-hadron states: " << params.param.states.size() << std::endl;
      QDPIO::cout << "Number of contractions: " << params.param.contractions.size() << std::endl;

      // Open output file
      std::ofstream ofs;
      if (Layout::primaryNode())
      {
        ofs.open(params.named_obj.output_file.c_str());
        ofs << "# Multi-hadron correlator output (Feynman-Hellmann)\n";
        ofs << "# t_eff = (t - t0 + Lt) % Lt, with t0 = " << params.param.t0 << "\n";
        ofs << "# Columns: state flavor t_eff re im\n";
        ofs << "#\n";
        ofs << "# Spin-1 states use Cartesian polarization basis (x, y, z):\n";
        ofs << "#   state_x = x-polarization\n";
        ofs << "#   state_y = y-polarization\n";
        ofs << "#   state_z = z-polarization\n";
        ofs << "#\n";
        ofs << "# DIBARYON (1/2 x 1/2 -> 1): e.g., deuteron\n";
        ofs << "#   C_x = C_y = 0.5*(G_up1*G_up2 + G_down1*G_down2)\n";
        ofs << "#   C_z = 0.5*(G_up1*G_down2 + G_down1*G_up2)\n";
        ofs << "#\n";
        ofs << "# DIMESON RHO-RHO (1 x 1 -> 1): e.g., di-rho\n";
        ofs << "#   Uses Levi-Civita (cross product) coupling for antisymmetric spin-1\n";
        ofs << "#\n";
        ofs << "# Sum relation: C(total) = C_x + C_y + C_z\n";
        ofs << "# Tensor polarization: b1 ~ C_z - 0.5*(C_x + C_y)\n";
        ofs << "#\n";
        ofs << "# BOOSTED states (nonzero total momentum P):\n";
        ofs << "#   state_A1 = A1 irrep (longitudinal, along P)\n";
        ofs << "#   state_Ey = E irrep, y-component (transverse)\n";
        ofs << "#   state_Ez = E irrep, z-component (transverse)\n";
      }

      // Loop over flavor combinations
      push(xml_out, "Correlators");

      for (int ic = 0; ic < params.param.contractions.size(); ++ic)
      {
        const Params::FlavorCombination_t& combo = params.param.contractions[ic];

        QDPIO::cout << "Computing contraction: " << combo.label
                    << " (up=" << (combo.up_is_perturbed ? "perturbed" : "unperturbed")
                    << ", down=" << (combo.down_is_perturbed ? "perturbed" : "unperturbed")
                    << ")" << std::endl;

        // Select propagators based on flavor assignment
        const LatticePropagator& up_prop = combo.up_is_perturbed ? perturbed_prop : unperturbed_prop;
        const LatticePropagator& down_prop = combo.down_is_perturbed ? perturbed_prop : unperturbed_prop;

        push(xml_out, "flavor_combination");
        write(xml_out, "label", combo.label);
        write(xml_out, "up_is_perturbed", combo.up_is_perturbed);
        write(xml_out, "down_is_perturbed", combo.down_is_perturbed);
        write(xml_out, "t0", params.param.t0);

        // Loop over multi-hadron states
        for (int is = 0; is < params.param.states.size(); ++is)
        {
          const Params::MultiHadronState_t& state = params.param.states[is];

          QDPIO::cout << "  Computing state: " << state.name
                      << " (" << state.hadron1 << " x " << state.hadron2
                      << ", S=1, Cartesian basis"
                      << (state.is_boosted ? ", BOOSTED" : "")
                      << ")" << std::endl;

          // Compute spin-projected multi-hadron correlator
          multi1d<DComplex> multi_hadron_corr;

          // Cartesian polarization correlators (x, y, z) for P=0
          multi1d<DComplex> corr_x, corr_y, corr_z;
          bool have_cartesian_states = false;

          // Irrep-projected correlators for boosted states
          multi1d<DComplex> corr_A1, corr_Ey, corr_Ez;
          bool have_boosted_irreps = false;

          // Reference single-hadron correlators for output
          multi2d<DComplex> hadron1_corr, hadron2_corr;

          if (state.type == "DIBARYON")
          {
            // Compute spin-resolved constituent baryon correlators
            multi2d<DComplex> h1_up, h1_down;
            multi2d<DComplex> h2_up, h2_down;

            if (state.hadron1 == "PROTON")
            {
              protonCorrelatorSpinResolved(up_prop, down_prop, phases, params.param.t0, h1_up, h1_down);
            }
            else if (state.hadron1 == "NEUTRON")
            {
              neutronCorrelatorSpinResolved(up_prop, down_prop, phases, params.param.t0, h1_up, h1_down);
            }
            else
            {
              QDPIO::cerr << "MULTI_HADRON_FH: Unknown baryon type for hadron1: " << state.hadron1 << std::endl;
              QDP_abort(1);
            }

            if (state.hadron2 == "PROTON")
            {
              protonCorrelatorSpinResolved(up_prop, down_prop, phases, params.param.t0, h2_up, h2_down);
            }
            else if (state.hadron2 == "NEUTRON")
            {
              neutronCorrelatorSpinResolved(up_prop, down_prop, phases, params.param.t0, h2_up, h2_down);
            }
            else
            {
              QDPIO::cerr << "MULTI_HADRON_FH: Unknown baryon type for hadron2: " << state.hadron2 << std::endl;
              QDP_abort(1);
            }

            if (state.is_boosted)
            {
              // Boosted dibaryon: C4v irrep projection
              dibaryonCorrelatorBoostedCartesian(h1_up, h1_down, h2_up, h2_down,
                                                 phases,
                                                 state.boost.fiducial_mom1,
                                                 state.boost.fiducial_mom2,
                                                 corr_A1, corr_Ey, corr_Ez);

              // Total spin-1 = A1 + Ey + Ez
              multi_hadron_corr.resize(Lt);
              for (int t = 0; t < Lt; ++t)
                multi_hadron_corr[t] = corr_A1[t] + corr_Ey[t] + corr_Ez[t];

              have_boosted_irreps = true;
            }
            else
            {
              // P=0 dibaryon: Cartesian polarization correlators
              dibaryonCorrelatorCartesian(h1_up, h1_down, h2_up, h2_down,
                                          phases, corr_x, corr_y, corr_z);

              // Total spin-1 correlator = sum of Cartesian components
              dibaryonCorrelatorSpin1(corr_x, corr_y, corr_z, multi_hadron_corr);
              have_cartesian_states = true;
            }

            // Compute spin-traced single-hadron correlators for reference output
            hadron1_corr.resize(h1_up.size1(), h1_up.size2());
            hadron2_corr.resize(h2_up.size1(), h2_up.size2());
            for (int p = 0; p < h1_up.size1(); ++p)
              for (int t = 0; t < h1_up.size2(); ++t)
              {
                hadron1_corr[p][t] = h1_up[p][t] + h1_down[p][t];
                hadron2_corr[p][t] = h2_up[p][t] + h2_down[p][t];
              }
          }
          else if (state.type == "DIMESON")
          {
            // Check if this is a vector-vector (rho-rho) dimeson
            bool is_rho1 = (state.hadron1 == "RHO" || state.hadron1 == "RHO_X" ||
                            state.hadron1 == "RHO_Y" || state.hadron1 == "RHO_Z");
            bool is_rho2 = (state.hadron2 == "RHO" || state.hadron2 == "RHO_X" ||
                            state.hadron2 == "RHO_Y" || state.hadron2 == "RHO_Z");

            if (is_rho1 && is_rho2)
            {
              // Vector-vector dimeson (rho-rho)
              // Compute rho correlators for all polarizations
              multi2d<DComplex> rho1_x, rho1_y, rho1_z;
              multi2d<DComplex> rho2_x, rho2_y, rho2_z;

              // For rho+ (u dbar): quark = up, antiquark = down
              // For rho- (d ubar): quark = down, antiquark = up
              // For rho0: (u ubar - d dbar)/sqrt(2) - we use u ubar for simplicity
              rhoMesonCorrelatorAllPol(up_prop, down_prop, phases, params.param.t0,
                                       rho1_x, rho1_y, rho1_z);
              rhoMesonCorrelatorAllPol(up_prop, down_prop, phases, params.param.t0,
                                       rho2_x, rho2_y, rho2_z);

              if (state.is_boosted)
              {
                // Boosted di-rho: C4v irrep projection
                diRhoCorrelatorBoostedCartesian(rho1_x, rho1_y, rho1_z,
                                                 rho2_x, rho2_y, rho2_z,
                                                 phases,
                                                 state.boost.fiducial_mom1,
                                                 state.boost.fiducial_mom2,
                                                 corr_A1, corr_Ey, corr_Ez);

                // Total spin-1 = A1 + Ey + Ez
                multi_hadron_corr.resize(Lt);
                for (int t = 0; t < Lt; ++t)
                  multi_hadron_corr[t] = corr_A1[t] + corr_Ey[t] + corr_Ez[t];

                have_boosted_irreps = true;
              }
              else
              {
                // P=0 di-rho: Cartesian polarization correlators
                diRhoCorrelatorCartesian(rho1_x, rho1_y, rho1_z,
                                         rho2_x, rho2_y, rho2_z,
                                         phases, corr_x, corr_y, corr_z);

                // Total spin-1 correlator = sum of Cartesian components
                diRhoCorrelatorSpin1(rho1_x, rho1_y, rho1_z,
                                     rho2_x, rho2_y, rho2_z,
                                     phases, multi_hadron_corr);
                have_cartesian_states = true;
              }

              // Use z-polarized rho as reference single-hadron correlator
              hadron1_corr = rho1_z;
              hadron2_corr = rho2_z;
            }
            else
            {
              QDPIO::cerr << "MULTI_HADRON_FH: Unsupported DIMESON hadron types: "
                          << state.hadron1 << " x " << state.hadron2 << std::endl;
              QDP_abort(1);
            }
          }
          else
          {
            QDPIO::cerr << "MULTI_HADRON_FH: Unknown state type: " << state.type << std::endl;
            QDP_abort(1);
          }

          // Write XML output
          push(xml_out, "state");
          write(xml_out, "name", state.name);
          write(xml_out, "type", state.type);
          write(xml_out, "hadron1", state.hadron1);
          write(xml_out, "hadron2", state.hadron2);
          write(xml_out, "is_boosted", state.is_boosted);

          // Write single-hadron correlators at zero momentum for reference
          push(xml_out, "hadron1_p0");
          for (int t = 0; t < Lt; ++t)
          {
            push(xml_out, "elem");
            write(xml_out, "t_eff", t);
            write(xml_out, "re", Real(real(hadron1_corr[0][t])));
            write(xml_out, "im", Real(imag(hadron1_corr[0][t])));
            pop(xml_out);
          }
          pop(xml_out);

          push(xml_out, "hadron2_p0");
          for (int t = 0; t < Lt; ++t)
          {
            push(xml_out, "elem");
            write(xml_out, "t_eff", t);
            write(xml_out, "re", Real(real(hadron2_corr[0][t])));
            write(xml_out, "im", Real(imag(hadron2_corr[0][t])));
            pop(xml_out);
          }
          pop(xml_out);

          // Write multi-hadron correlator (summed over all irreps/polarizations)
          push(xml_out, "multi_hadron");
          for (int t = 0; t < Lt; ++t)
          {
            push(xml_out, "elem");
            write(xml_out, "t_eff", t);
            write(xml_out, "re", Real(real(multi_hadron_corr[t])));
            write(xml_out, "im", Real(imag(multi_hadron_corr[t])));
            pop(xml_out);
          }
          pop(xml_out);

          // Write irrep-projected correlators for boosted states
          if (have_boosted_irreps)
          {
            // A1 (longitudinal) irrep correlator
            push(xml_out, "irrep_A1");
            for (int t = 0; t < Lt; ++t)
            {
              push(xml_out, "elem");
              write(xml_out, "t_eff", t);
              write(xml_out, "re", Real(real(corr_A1[t])));
              write(xml_out, "im", Real(imag(corr_A1[t])));
              pop(xml_out);
            }
            pop(xml_out);

            // E irrep y-component (transverse)
            push(xml_out, "irrep_Ey");
            for (int t = 0; t < Lt; ++t)
            {
              push(xml_out, "elem");
              write(xml_out, "t_eff", t);
              write(xml_out, "re", Real(real(corr_Ey[t])));
              write(xml_out, "im", Real(imag(corr_Ey[t])));
              pop(xml_out);
            }
            pop(xml_out);

            // E irrep z-component (transverse)
            push(xml_out, "irrep_Ez");
            for (int t = 0; t < Lt; ++t)
            {
              push(xml_out, "elem");
              write(xml_out, "t_eff", t);
              write(xml_out, "re", Real(real(corr_Ez[t])));
              write(xml_out, "im", Real(imag(corr_Ez[t])));
              pop(xml_out);
            }
            pop(xml_out);
          }

          // Write Cartesian polarization correlators for P=0 spin-1 states
          if (have_cartesian_states)
          {
            // x-polarization correlator
            push(xml_out, "pol_x");
            for (int t = 0; t < Lt; ++t)
            {
              push(xml_out, "elem");
              write(xml_out, "t_eff", t);
              write(xml_out, "re", Real(real(corr_x[t])));
              write(xml_out, "im", Real(imag(corr_x[t])));
              pop(xml_out);
            }
            pop(xml_out);

            // y-polarization correlator
            push(xml_out, "pol_y");
            for (int t = 0; t < Lt; ++t)
            {
              push(xml_out, "elem");
              write(xml_out, "t_eff", t);
              write(xml_out, "re", Real(real(corr_y[t])));
              write(xml_out, "im", Real(imag(corr_y[t])));
              pop(xml_out);
            }
            pop(xml_out);

            // z-polarization correlator
            push(xml_out, "pol_z");
            for (int t = 0; t < Lt; ++t)
            {
              push(xml_out, "elem");
              write(xml_out, "t_eff", t);
              write(xml_out, "re", Real(real(corr_z[t])));
              write(xml_out, "im", Real(imag(corr_z[t])));
              pop(xml_out);
            }
            pop(xml_out);
          }

          pop(xml_out);  // state

          // Write plain text output
          if (Layout::primaryNode())
          {
            // Summed correlator
            for (int t = 0; t < Lt; ++t)
            {
              ofs << state.name << " " << combo.label << " " << t
                  << " " << real(multi_hadron_corr[t]) << " " << imag(multi_hadron_corr[t])
                  << "\n";
            }

            // Irrep-projected correlators for boosted states
            if (have_boosted_irreps)
            {
              for (int t = 0; t < Lt; ++t)
              {
                ofs << state.name << "_A1 " << combo.label << " " << t
                    << " " << real(corr_A1[t]) << " " << imag(corr_A1[t])
                    << "\n";
              }
              for (int t = 0; t < Lt; ++t)
              {
                ofs << state.name << "_Ey " << combo.label << " " << t
                    << " " << real(corr_Ey[t]) << " " << imag(corr_Ey[t])
                    << "\n";
              }
              for (int t = 0; t < Lt; ++t)
              {
                ofs << state.name << "_Ez " << combo.label << " " << t
                    << " " << real(corr_Ez[t]) << " " << imag(corr_Ez[t])
                    << "\n";
              }
            }

            // Individual Cartesian polarization correlators for P=0
            if (have_cartesian_states)
            {
              for (int t = 0; t < Lt; ++t)
              {
                ofs << state.name << "_x " << combo.label << " " << t
                    << " " << real(corr_x[t]) << " " << imag(corr_x[t])
                    << "\n";
              }
              for (int t = 0; t < Lt; ++t)
              {
                ofs << state.name << "_y " << combo.label << " " << t
                    << " " << real(corr_y[t]) << " " << imag(corr_y[t])
                    << "\n";
              }
              for (int t = 0; t < Lt; ++t)
              {
                ofs << state.name << "_z " << combo.label << " " << t
                    << " " << real(corr_z[t]) << " " << imag(corr_z[t])
                    << "\n";
              }
            }
          }
        }

        pop(xml_out);  // flavor_combination
      }

      pop(xml_out);  // Correlators

      if (Layout::primaryNode())
      {
        ofs.close();
        QDPIO::cout << "Correlators written to " << params.named_obj.output_file << std::endl;
      }

      // Timing summary
      snoop.stop();
      write(xml_out, "total_time_sec", snoop.getTimeInSeconds());

      QDPIO::cout << "MULTI_HADRON_FH: total time = "
                  << snoop.getTimeInSeconds() << " sec" << std::endl;

      pop(xml_out);  // MultiHadronFH

      END_CODE();
    }

  } // namespace InlineMultiHadronFHEnv

} // namespace Chroma
