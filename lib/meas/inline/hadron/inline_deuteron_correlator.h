// -*- C++ -*-
/**
 * @file inline_deuteron_correlator.h
 * @brief Inline measurement for deuteron correlators (contraction-only)
 *
 * Computes proton, neutron, deuteron, and dineutron correlators
 * from pre-computed propagators with specified flavor combinations.
 */

#ifndef __inline_deuteron_correlator_h__
#define __inline_deuteron_correlator_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  namespace InlineDeuteronCorrelatorEnv
  {
    //! Registration
    bool registerAll();

    //! Parameters
    struct Params
    {
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

  } // namespace InlineDeuteronCorrelatorEnv

} // namespace Chroma

#endif
