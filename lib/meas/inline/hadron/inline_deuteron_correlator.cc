/**
 * @file inline_deuteron_correlator.cc
 * @brief Inline measurement for deuteron correlators (contraction-only)
 *
 * Computes proton, neutron, deuteron, and dineutron correlators
 * from pre-computed propagators with specified flavor combinations.
 */

#include "inline_deuteron_correlator.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"

#include <fstream>

namespace Chroma
{
  namespace InlineDeuteronCorrelatorEnv
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
    const std::string measurement_name = "DEUTERON_CORRELATOR";

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
            param.contractions[0].label = "up_perturbed";
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
    // Proton correlator computation
    //------------------------------------------------------------------

    //! Compute proton two-point function with separate up/down propagators
    /*!
     * Proton interpolator: epsilon_abc (u^a C gamma_5 d^b) u^c
     *
     * \param up_prop    Propagator for up quarks (Read)
     * \param down_prop  Propagator for down quarks (Read)
     * \param phases     Momentum phases (Read)
     * \param t0         Source timeslice for t_eff calculation (Read)
     * \param corr       Correlator output indexed by t_eff (Write)
     */
    void protonCorrelator(
      const LatticePropagator& up_prop,
      const LatticePropagator& down_prop,
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
      SpinMatrix Cg5 = Gamma(10) * (Gamma(15) * g_one);  // C = gamma_4 * gamma_2

      // Proton = epsilon_abc (u^a C gamma_5 d^b) u^c
      // Diquark: (u C gamma_5 d) - contracts color indices a,b
      // Third quark: u^c
      LatticePropagator diquark = quarkContract13(up_prop * Cg5,
                                                   Gamma(15) * down_prop);

      // Positive parity projection
      SpinMatrix parity_proj = 0.5 * (g_one + Gamma(8) * g_one);

      // Contract diquark with third up quark and trace
      LatticeComplex proton_field = trace(parity_proj * traceColor(diquark * up_prop));

      // Momentum projection
      multi2d<DComplex> hsum = phases.sft(proton_field);

      // Apply t_eff shift: t_eff = (t - t0 + Lt) % Lt
      for (int p = 0; p < num_mom; ++p)
        for (int t = 0; t < Lt; ++t)
        {
          int t_eff = (t - t0 + Lt) % Lt;
          corr[p][t_eff] = hsum[p][t];
        }

      END_CODE();
    }


    //! Compute two-nucleon (deuteron/dineutron) correlators
    void twoNucleonCorrelator(
      const multi2d<DComplex>& nucleon_corr,
      const SftMom& phases,
      multi1d<DComplex>& deuteron_corr)
    {
      START_CODE();

      int num_mom = phases.numMom();
      int Lt = nucleon_corr.size2();

      deuteron_corr.resize(Lt);
      deuteron_corr = zero;

      // For each momentum p, find -p and accumulate N(p) * N(-p)
      for (int ip = 0; ip < num_mom; ++ip)
      {
        multi1d<int> mom_p = phases.numToMom(ip);

        int ip_neg = -1;
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
          {
            ip_neg = jp;
            break;
          }
        }

        if (ip_neg < 0) continue;

        for (int t = 0; t < Lt; ++t)
        {
          deuteron_corr[t] += nucleon_corr[ip][t] * nucleon_corr[ip_neg][t];
        }
      }

      // Normalize
      for (int t = 0; t < Lt; ++t)
      {
        deuteron_corr[t] *= 0.5;
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

      push(xml_out, "DeuteronCorrelator");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << "DEUTERON_CORRELATOR: Deuteron correlator measurement (contraction-only)" << std::endl;

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
        QDPIO::cerr << "DEUTERON_CORRELATOR: perturbed_prop_id cast failed" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "DEUTERON_CORRELATOR: perturbed_prop_id not found: " << e << std::endl;
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
        QDPIO::cerr << "DEUTERON_CORRELATOR: unperturbed_prop_id cast failed" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "DEUTERON_CORRELATOR: unperturbed_prop_id not found: " << e << std::endl;
        QDP_abort(1);
      }

      // Setup momentum projection
      SftMom phases(params.param.mom2_max, false, Nd-1);
      int num_mom = phases.numMom();
      int Lt = Layout::lattSize()[Nd-1];

      QDPIO::cout << "Number of momenta: " << num_mom << std::endl;
      QDPIO::cout << "Lattice time extent: " << Lt << std::endl;
      QDPIO::cout << "Source timeslice t0: " << params.param.t0 << std::endl;
      QDPIO::cout << "Number of contractions: " << params.param.contractions.size() << std::endl;

      // Open output file
      std::ofstream ofs;
      if (Layout::primaryNode())
      {
        ofs.open(params.named_obj.output_file.c_str());
        ofs << "# Deuteron correlator output (contraction-only)\n";
        ofs << "# t_eff = (t - t0 + Lt) % Lt, with t0 = " << params.param.t0 << "\n";
        ofs << "# Columns: label t_eff px py pz proton_re proton_im deuteron_re deuteron_im\n";
      }

      // Loop over flavor combinations
      push(xml_out, "Contractions");

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

        // Compute proton correlator (indexed by t_eff)
        multi2d<DComplex> proton_corr;
        protonCorrelator(up_prop, down_prop, phases, params.param.t0, proton_corr);

        // Compute deuteron correlator
        multi1d<DComplex> deuteron_corr;
        twoNucleonCorrelator(proton_corr, phases, deuteron_corr);

        // Write XML output
        push(xml_out, "contraction");
        write(xml_out, "label", combo.label);
        write(xml_out, "up_is_perturbed", combo.up_is_perturbed);
        write(xml_out, "down_is_perturbed", combo.down_is_perturbed);
        write(xml_out, "t0", params.param.t0);

        // Write proton correlators for all momenta
        push(xml_out, "proton_correlators");
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
            write(xml_out, "re", Real(real(proton_corr[p][t])));
            write(xml_out, "im", Real(imag(proton_corr[p][t])));
            pop(xml_out);
          }
          pop(xml_out);  // correlator

          pop(xml_out);  // momentum
        }
        pop(xml_out);  // proton_correlators

        // Write deuteron correlator (summed over momentum pairs)
        push(xml_out, "deuteron");
        for (int t = 0; t < Lt; ++t)
        {
          push(xml_out, "elem");
          write(xml_out, "t_eff", t);
          write(xml_out, "re", Real(real(deuteron_corr[t])));
          write(xml_out, "im", Real(imag(deuteron_corr[t])));
          pop(xml_out);
        }
        pop(xml_out);

        pop(xml_out);  // contraction

        // Write plain text output for all momenta
        if (Layout::primaryNode())
        {
          for (int p = 0; p < num_mom; ++p)
          {
            multi1d<int> mom = phases.numToMom(p);
            for (int t = 0; t < Lt; ++t)
            {
              // Output: label t_eff px py pz proton_re proton_im deuteron_re deuteron_im
              // Note: deuteron is the same for all momenta (summed over p pairs)
              ofs << combo.label << " " << t
                  << " " << mom[0] << " " << mom[1] << " " << mom[2]
                  << " " << real(proton_corr[p][t]) << " " << imag(proton_corr[p][t])
                  << " " << real(deuteron_corr[t]) << " " << imag(deuteron_corr[t])
                  << "\n";
            }
          }
        }
      }

      pop(xml_out);  // Contractions

      if (Layout::primaryNode())
      {
        ofs.close();
        QDPIO::cout << "Correlators written to " << params.named_obj.output_file << std::endl;
      }

      // Timing summary
      snoop.stop();
      write(xml_out, "total_time_sec", snoop.getTimeInSeconds());

      QDPIO::cout << "DEUTERON_CORRELATOR: total time = "
                  << snoop.getTimeInSeconds() << " sec" << std::endl;

      pop(xml_out);  // DeuteronCorrelator

      END_CODE();
    }

  } // namespace InlineDeuteronCorrelatorEnv

} // namespace Chroma
