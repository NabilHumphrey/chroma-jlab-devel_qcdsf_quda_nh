// -*- C++ -*-
/**
 * @file inline_meson_spectrum_fh.h
 * @brief Inline measurement for meson correlators with Feynman-Hellmann support
 *
 * Computes meson correlators from pre-computed propagators with
 * configurable gamma structures and flavor combinations for FH analysis.
 */

#ifndef __inline_meson_spectrum_fh_h__
#define __inline_meson_spectrum_fh_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma
{
  namespace InlineMesonSpectrumFHEnv
  {
    //! Registration
    bool registerAll();

    //! Parameters
    struct Params
    {
      //! Meson interpolator specification
      struct MesonInterpolator_t
      {
        std::string name;             //!< Label for this meson (e.g., "rho_x")
        int gamma;                    //!< Gamma matrix index (1=γ₁, 2=γ₂, 4=γ₃, etc.)
      };

      //! Flavor combination for contractions
      struct FlavorCombination_t
      {
        bool quark_is_perturbed;      //!< Use perturbed prop for quark
        bool antiquark_is_perturbed;  //!< Use perturbed prop for antiquark
        std::string label;            //!< Label for this combination in output
      };

      //! General parameters
      struct Param_t
      {
        int mom2_max;                 //!< Maximum momentum squared
        bool avg_equiv_mom;           //!< Average equivalent momenta
        int t0;                       //!< Source timeslice for t_eff calculation

        //! List of meson interpolators to compute
        multi1d<MesonInterpolator_t> mesons;

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

  } // namespace InlineMesonSpectrumFHEnv

} // namespace Chroma

#endif
