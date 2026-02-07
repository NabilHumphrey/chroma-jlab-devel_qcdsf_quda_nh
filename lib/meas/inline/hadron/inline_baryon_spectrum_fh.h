// -*- C++ -*-
/**
 * @file inline_baryon_spectrum_fh.h
 * @brief Inline measurement for baryon correlators with Feynman-Hellmann support
 *
 * Computes baryon correlators from pre-computed propagators with
 * configurable interpolators and flavor combinations for FH analysis.
 */

#ifndef __inline_baryon_spectrum_fh_h__
#define __inline_baryon_spectrum_fh_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma
{
  namespace InlineBaryonSpectrumFHEnv
  {
    //! Registration
    bool registerAll();

    //! Parameters
    struct Params
    {
      //! Baryon interpolator specification
      struct BaryonInterpolator_t
      {
        std::string name;             //!< Label for this baryon (e.g., "proton")
        std::string type;             //!< Baryon type: PROTON, NEUTRON
        std::string parity;           //!< Parity projection: PLUS, MINUS, UNPROJECTED
      };

      //! Flavor combination for contractions
      struct FlavorCombination_t
      {
        bool up_is_perturbed;         //!< Use perturbed prop for up quarks
        bool down_is_perturbed;       //!< Use perturbed prop for down quarks
        std::string label;            //!< Label for this combination in output
      };

      //! General parameters
      struct Param_t
      {
        int mom2_max;                 //!< Maximum momentum squared
        bool avg_equiv_mom;           //!< Average equivalent momenta
        int t0;                       //!< Source timeslice for t_eff calculation

        //! List of baryon interpolators to compute
        multi1d<BaryonInterpolator_t> baryons;

        //! List of flavor combinations to compute
        multi1d<FlavorCombination_t> contractions;
      };

      //! Named object parameters
      struct NamedObject_t
      {
        std::string gauge_id;             //!< Gauge field ID
        std::string perturbed_prop_id;    //!< Perturbed (FH) propagator
        std::string unperturbed_prop_id;  //!< Unperturbed propagator
        std::string output_file;          //!< Output filename
      };

      unsigned long frequency;
      Param_t param;
      NamedObject_t named_obj;

      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path) const;
    };

    //! Inline measurement
    class InlineMeas : public AbsInlineMeasurement
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const { return params.frequency; }

      void operator()(unsigned long update_no, XMLWriter& xml_out);

    private:
      Params params;
    };

  } // namespace InlineBaryonSpectrumFHEnv

} // namespace Chroma

#endif
