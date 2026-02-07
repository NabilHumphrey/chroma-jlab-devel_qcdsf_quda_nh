/**
 * @file inline_meson_spectrum_fh.cc
 * @brief Inline measurement for meson correlators with Feynman-Hellmann support
 *
 * Computes meson correlators from pre-computed propagators with
 * configurable gamma structures and flavor combinations for FH analysis.
 */

#include "inline_meson_spectrum_fh.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"

#include <fstream>

namespace Chroma
{
  namespace InlineMesonSpectrumFHEnv
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
    const std::string measurement_name = "MESON_SPECTRUM_FH";

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
    // MesonInterpolator_t I/O
    //------------------------------------------------------------------
    void read(XMLReader& xml, const std::string& path, Params::MesonInterpolator_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "name", param.name);
      read(paramtop, "gamma", param.gamma);
    }

    void write(XMLWriter& xml, const std::string& path, const Params::MesonInterpolator_t& param)
    {
      push(xml, path);
      write(xml, "name", param.name);
      write(xml, "gamma", param.gamma);
      pop(xml);
    }


    //------------------------------------------------------------------
    // FlavorCombination_t I/O
    //------------------------------------------------------------------
    void read(XMLReader& xml, const std::string& path, Params::FlavorCombination_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "quark_is_perturbed", param.quark_is_perturbed);
      read(paramtop, "antiquark_is_perturbed", param.antiquark_is_perturbed);
      read(paramtop, "label", param.label);
    }

    void write(XMLWriter& xml, const std::string& path, const Params::FlavorCombination_t& param)
    {
      push(xml, path);
      write(xml, "quark_is_perturbed", param.quark_is_perturbed);
      write(xml, "antiquark_is_perturbed", param.antiquark_is_perturbed);
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

          // Read meson interpolators list
          if (ptop.count("mesons") == 1)
          {
            XMLReader mestop(ptop, "mesons");
            int num_mesons = mestop.count("elem");
            param.mesons.resize(num_mesons);
            for (int i = 0; i < num_mesons; ++i)
            {
              read(mestop, "elem[" + std::to_string(i+1) + "]", param.mesons[i]);
            }
          }
          else
          {
            // Default: rho meson in all three polarizations
            param.mesons.resize(3);
            param.mesons[0].name = "rho_x";
            param.mesons[0].gamma = 1;
            param.mesons[1].name = "rho_y";
            param.mesons[1].gamma = 2;
            param.mesons[2].name = "rho_z";
            param.mesons[2].gamma = 4;
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
            // Default: single contraction with quark=perturbed, antiquark=unperturbed
            param.contractions.resize(1);
            param.contractions[0].quark_is_perturbed = true;
            param.contractions[0].antiquark_is_perturbed = false;
            param.contractions[0].label = "quark_FH";
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

      push(xml_out, "mesons");
      for (int i = 0; i < param.mesons.size(); ++i)
      {
        write(xml_out, "elem", param.mesons[i]);
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
    // Meson correlator computation
    //------------------------------------------------------------------

    //! Compute meson two-point function
    /*!
     * Meson correlator: Tr[ Γ S_q(x,0) Γ γ₅ S_qbar†(x,0) γ₅ ]
     *
     * For a meson with interpolator qbar Γ q, the correlator is:
     * C(t) = sum_x < 0 | (qbar Γ q)(x,t) (qbar Γ q)†(0,0) | 0 >
     *
     * \param quark_prop     Propagator for quark (Read)
     * \param antiquark_prop Propagator for antiquark (Read)
     * \param gamma          Gamma matrix index for interpolator (Read)
     * \param phases         Momentum phases (Read)
     * \param t0             Source timeslice for t_eff calculation (Read)
     * \param corr           Correlator output indexed by t_eff (Write)
     */
    void mesonCorrelator(
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

      // Construct anti-quark propagator: γ₅ S† γ₅
      int G5 = Ns*Ns - 1;  // γ₅ = 15
      LatticePropagator anti_quark = Gamma(G5) * adj(antiquark_prop) * Gamma(G5);

      // Meson correlation function: Tr[ Γ S_q Γ (γ₅ S_qbar† γ₅) ]
      // = Tr[ Γ S_q Γ anti_quark ]
      LatticeComplex corr_fn = trace(Gamma(gamma) * quark_prop * Gamma(gamma) * anti_quark);

      // Momentum projection
      multi2d<DComplex> hsum = phases.sft(corr_fn);

      // Apply t_eff shift: t_eff = (t - t0 + Lt) % Lt
      for (int p = 0; p < num_mom; ++p)
        for (int t = 0; t < Lt; ++t)
        {
          int t_eff = (t - t0 + Lt) % Lt;
          corr[p][t_eff] = hsum[p][t];
        }

      END_CODE();
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

      push(xml_out, "MesonSpectrumFH");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << "MESON_SPECTRUM_FH: Meson correlator measurement with FH support" << std::endl;

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
        QDPIO::cerr << "MESON_SPECTRUM_FH: perturbed_prop_id cast failed" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "MESON_SPECTRUM_FH: perturbed_prop_id not found: " << e << std::endl;
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
        QDPIO::cerr << "MESON_SPECTRUM_FH: unperturbed_prop_id cast failed" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "MESON_SPECTRUM_FH: unperturbed_prop_id not found: " << e << std::endl;
        QDP_abort(1);
      }

      // Setup momentum projection
      SftMom phases(params.param.mom2_max, false, Nd-1);
      int num_mom = phases.numMom();
      int Lt = Layout::lattSize()[Nd-1];

      QDPIO::cout << "Number of momenta: " << num_mom << std::endl;
      QDPIO::cout << "Lattice time extent: " << Lt << std::endl;
      QDPIO::cout << "Source timeslice t0: " << params.param.t0 << std::endl;
      QDPIO::cout << "Number of meson types: " << params.param.mesons.size() << std::endl;
      QDPIO::cout << "Number of contractions: " << params.param.contractions.size() << std::endl;

      // Open output file
      std::ofstream ofs;
      if (Layout::primaryNode())
      {
        ofs.open(params.named_obj.output_file.c_str());
        ofs << "# Meson correlator output (Feynman-Hellmann)\n";
        ofs << "# t_eff = (t - t0 + Lt) % Lt, with t0 = " << params.param.t0 << "\n";
        ofs << "# Columns: meson flavor t_eff px py pz re im\n";
      }

      // Loop over flavor combinations
      push(xml_out, "Correlators");

      for (int ic = 0; ic < params.param.contractions.size(); ++ic)
      {
        const Params::FlavorCombination_t& combo = params.param.contractions[ic];

        QDPIO::cout << "Computing contraction: " << combo.label
                    << " (quark=" << (combo.quark_is_perturbed ? "perturbed" : "unperturbed")
                    << ", antiquark=" << (combo.antiquark_is_perturbed ? "perturbed" : "unperturbed")
                    << ")" << std::endl;

        // Select propagators based on flavor assignment
        const LatticePropagator& quark_prop = combo.quark_is_perturbed ? perturbed_prop : unperturbed_prop;
        const LatticePropagator& antiquark_prop = combo.antiquark_is_perturbed ? perturbed_prop : unperturbed_prop;

        push(xml_out, "flavor_combination");
        write(xml_out, "label", combo.label);
        write(xml_out, "quark_is_perturbed", combo.quark_is_perturbed);
        write(xml_out, "antiquark_is_perturbed", combo.antiquark_is_perturbed);
        write(xml_out, "t0", params.param.t0);

        // Loop over meson interpolators
        for (int im = 0; im < params.param.mesons.size(); ++im)
        {
          const Params::MesonInterpolator_t& meson = params.param.mesons[im];

          QDPIO::cout << "  Computing meson: " << meson.name << " (gamma=" << meson.gamma << ")" << std::endl;

          // Compute meson correlator (indexed by t_eff)
          multi2d<DComplex> meson_corr;
          mesonCorrelator(quark_prop, antiquark_prop, meson.gamma, phases, params.param.t0, meson_corr);

          // Write XML output
          push(xml_out, "meson");
          write(xml_out, "name", meson.name);
          write(xml_out, "gamma", meson.gamma);

          push(xml_out, "correlators");
          for (int p = 0; p < num_mom; ++p)
          {
            push(xml_out, "momentum");
            multi1d<int> mom = phases.numToMom(p);
            write(xml_out, "mom_num", p);
            write(xml_out, "px", mom[0]);
            write(xml_out, "py", mom[1]);
            write(xml_out, "pz", mom[2]);

            push(xml_out, "correlator");
            for (int t = 0; t < Lt; ++t)
            {
              push(xml_out, "elem");
              write(xml_out, "t_eff", t);
              write(xml_out, "re", Real(real(meson_corr[p][t])));
              write(xml_out, "im", Real(imag(meson_corr[p][t])));
              pop(xml_out);
            }
            pop(xml_out);  // correlator

            pop(xml_out);  // momentum
          }
          pop(xml_out);  // correlators

          pop(xml_out);  // meson

          // Write plain text output for all momenta
          if (Layout::primaryNode())
          {
            for (int p = 0; p < num_mom; ++p)
            {
              multi1d<int> mom = phases.numToMom(p);
              for (int t = 0; t < Lt; ++t)
              {
                ofs << meson.name << " " << combo.label << " " << t
                    << " " << mom[0] << " " << mom[1] << " " << mom[2]
                    << " " << real(meson_corr[p][t]) << " " << imag(meson_corr[p][t])
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

      QDPIO::cout << "MESON_SPECTRUM_FH: total time = "
                  << snoop.getTimeInSeconds() << " sec" << std::endl;

      pop(xml_out);  // MesonSpectrumFH

      END_CODE();
    }

  } // namespace InlineMesonSpectrumFHEnv

} // namespace Chroma
