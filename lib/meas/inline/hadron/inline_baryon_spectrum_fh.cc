/**
 * @file inline_baryon_spectrum_fh.cc
 * @brief Inline measurement for baryon correlators with Feynman-Hellmann support
 *
 * Computes baryon correlators from pre-computed propagators with
 * configurable interpolators and flavor combinations for FH analysis.
 */

#include "inline_baryon_spectrum_fh.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"

#include <fstream>

namespace Chroma
{
  namespace InlineBaryonSpectrumFHEnv
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
    const std::string measurement_name = "BARYON_SPECTRUM_FH";

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
    // BaryonInterpolator_t I/O
    //------------------------------------------------------------------
    void read(XMLReader& xml, const std::string& path, Params::BaryonInterpolator_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "name", param.name);
      read(paramtop, "type", param.type);

      if (paramtop.count("parity") == 1)
        read(paramtop, "parity", param.parity);
      else
        param.parity = "PLUS";
    }

    void write(XMLWriter& xml, const std::string& path, const Params::BaryonInterpolator_t& param)
    {
      push(xml, path);
      write(xml, "name", param.name);
      write(xml, "type", param.type);
      write(xml, "parity", param.parity);
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

          // Read baryon interpolators list
          if (ptop.count("baryons") == 1)
          {
            XMLReader bartop(ptop, "baryons");
            int num_baryons = bartop.count("elem");
            param.baryons.resize(num_baryons);
            for (int i = 0; i < num_baryons; ++i)
            {
              read(bartop, "elem[" + std::to_string(i+1) + "]", param.baryons[i]);
            }
          }
          else
          {
            // Default: proton with positive parity
            param.baryons.resize(1);
            param.baryons[0].name = "proton";
            param.baryons[0].type = "PROTON";
            param.baryons[0].parity = "PLUS";
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
            // Default: single contraction with up=perturbed, down=unperturbed
            param.contractions.resize(1);
            param.contractions[0].up_is_perturbed = true;
            param.contractions[0].down_is_perturbed = false;
            param.contractions[0].label = "up_FH";
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

      push(xml_out, "baryons");
      for (int i = 0; i < param.baryons.size(); ++i)
      {
        write(xml_out, "elem", param.baryons[i]);
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
    // Baryon correlator computation
    //------------------------------------------------------------------

    //! Compute proton two-point function with separate up/down propagators
    /*!
     * Proton interpolator: epsilon_abc (u^a C gamma_5 d^b) u^c
     *
     * \param up_prop    Propagator for up quarks (Read)
     * \param down_prop  Propagator for down quarks (Read)
     * \param parity     Parity projection: "PLUS", "MINUS", or "UNPROJECTED" (Read)
     * \param phases     Momentum phases (Read)
     * \param t0         Source timeslice for t_eff calculation (Read)
     * \param corr       Correlator output indexed by t_eff (Write)
     */
    void protonCorrelator(
      const LatticePropagator& up_prop,
      const LatticePropagator& down_prop,
      const std::string& parity,
      const SftMom& phases,
      int t0,
      multi2d<DComplex>& corr)
    {
      START_CODE();

      int num_mom = phases.numMom();
      int Lt = phases.numSubsets();

      corr.resize(num_mom, Lt);

      // C gamma_5 matrix for diquark
      SpinMatrix g_one = 1.0;
      SpinMatrix Cg5 = Gamma(10) * (Gamma(15) * g_one);  // C = gamma_4 * gamma_2, Cg5 = C * gamma_5

      // Proton = epsilon_abc (u^a C gamma_5 d^b) u^c
      // Diquark: (u C gamma_5 d) - contracts color indices a,b
      // Third quark: u^c
      LatticePropagator diquark = quarkContract13(up_prop * Cg5,
                                                   Gamma(15) * down_prop);

      // Parity projection
      SpinMatrix parity_proj;
      if (parity == "PLUS")
      {
        parity_proj = 0.5 * (g_one + Gamma(8) * g_one);  // (1 + gamma_4)/2
      }
      else if (parity == "MINUS")
      {
        parity_proj = 0.5 * (g_one - Gamma(8) * g_one);  // (1 - gamma_4)/2
      }
      else  // UNPROJECTED
      {
        parity_proj = g_one;
      }

      // Contract diquark with third up quark and trace
      LatticeComplex baryon_field = trace(parity_proj * traceColor(diquark * up_prop));

      // Momentum projection
      multi2d<DComplex> hsum = phases.sft(baryon_field);

      // Apply t_eff shift: t_eff = (t - t0 + Lt) % Lt
      for (int p = 0; p < num_mom; ++p)
        for (int t = 0; t < Lt; ++t)
        {
          int t_eff = (t - t0 + Lt) % Lt;
          corr[p][t_eff] = hsum[p][t];
        }

      END_CODE();
    }


    //! Compute neutron two-point function with separate up/down propagators
    /*!
     * Neutron interpolator: epsilon_abc (d^a C gamma_5 u^b) d^c
     * (just swap u <-> d from proton)
     */
    void neutronCorrelator(
      const LatticePropagator& up_prop,
      const LatticePropagator& down_prop,
      const std::string& parity,
      const SftMom& phases,
      int t0,
      multi2d<DComplex>& corr)
    {
      // Neutron is proton with u <-> d swap
      protonCorrelator(down_prop, up_prop, parity, phases, t0, corr);
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

      push(xml_out, "BaryonSpectrumFH");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << "BARYON_SPECTRUM_FH: Baryon correlator measurement with FH support" << std::endl;

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
        QDPIO::cerr << "BARYON_SPECTRUM_FH: perturbed_prop_id cast failed" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "BARYON_SPECTRUM_FH: perturbed_prop_id not found: " << e << std::endl;
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
        QDPIO::cerr << "BARYON_SPECTRUM_FH: unperturbed_prop_id cast failed" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "BARYON_SPECTRUM_FH: unperturbed_prop_id not found: " << e << std::endl;
        QDP_abort(1);
      }

      // Setup momentum projection
      SftMom phases(params.param.mom2_max, false, Nd-1);
      int num_mom = phases.numMom();
      int Lt = Layout::lattSize()[Nd-1];

      QDPIO::cout << "Number of momenta: " << num_mom << std::endl;
      QDPIO::cout << "Lattice time extent: " << Lt << std::endl;
      QDPIO::cout << "Source timeslice t0: " << params.param.t0 << std::endl;
      QDPIO::cout << "Number of baryon types: " << params.param.baryons.size() << std::endl;
      QDPIO::cout << "Number of contractions: " << params.param.contractions.size() << std::endl;

      // Open output file
      std::ofstream ofs;
      if (Layout::primaryNode())
      {
        ofs.open(params.named_obj.output_file.c_str());
        ofs << "# Baryon correlator output (Feynman-Hellmann)\n";
        ofs << "# t_eff = (t - t0 + Lt) % Lt, with t0 = " << params.param.t0 << "\n";
        ofs << "# Columns: baryon flavor t_eff px py pz re im\n";
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

        // Loop over baryon interpolators
        for (int ib = 0; ib < params.param.baryons.size(); ++ib)
        {
          const Params::BaryonInterpolator_t& baryon = params.param.baryons[ib];

          QDPIO::cout << "  Computing baryon: " << baryon.name
                      << " (type=" << baryon.type << ", parity=" << baryon.parity << ")" << std::endl;

          // Compute baryon correlator (indexed by t_eff)
          multi2d<DComplex> baryon_corr;

          if (baryon.type == "PROTON")
          {
            protonCorrelator(up_prop, down_prop, baryon.parity, phases, params.param.t0, baryon_corr);
          }
          else if (baryon.type == "NEUTRON")
          {
            neutronCorrelator(up_prop, down_prop, baryon.parity, phases, params.param.t0, baryon_corr);
          }
          else
          {
            QDPIO::cerr << "BARYON_SPECTRUM_FH: Unknown baryon type: " << baryon.type << std::endl;
            QDP_abort(1);
          }

          // Write XML output
          push(xml_out, "baryon");
          write(xml_out, "name", baryon.name);
          write(xml_out, "type", baryon.type);
          write(xml_out, "parity", baryon.parity);

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
              write(xml_out, "re", Real(real(baryon_corr[p][t])));
              write(xml_out, "im", Real(imag(baryon_corr[p][t])));
              pop(xml_out);
            }
            pop(xml_out);  // correlator

            pop(xml_out);  // momentum
          }
          pop(xml_out);  // correlators

          pop(xml_out);  // baryon

          // Write plain text output for all momenta
          if (Layout::primaryNode())
          {
            for (int p = 0; p < num_mom; ++p)
            {
              multi1d<int> mom = phases.numToMom(p);
              for (int t = 0; t < Lt; ++t)
              {
                ofs << baryon.name << " " << combo.label << " " << t
                    << " " << mom[0] << " " << mom[1] << " " << mom[2]
                    << " " << real(baryon_corr[p][t]) << " " << imag(baryon_corr[p][t])
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

      QDPIO::cout << "BARYON_SPECTRUM_FH: total time = "
                  << snoop.getTimeInSeconds() << " sec" << std::endl;

      pop(xml_out);  // BaryonSpectrumFH

      END_CODE();
    }

  } // namespace InlineBaryonSpectrumFHEnv

} // namespace Chroma
