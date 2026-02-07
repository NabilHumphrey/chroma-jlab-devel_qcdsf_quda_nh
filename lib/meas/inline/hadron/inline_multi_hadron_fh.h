// -*- C++ -*-
/**
 * @file inline_multi_hadron_fh.h
 * @brief Inline measurement for multi-hadron correlators with Feynman-Hellmann support
 *
 * Computes multi-hadron correlators (dibaryons, dimesons) from pre-computed
 * propagators using product ansatz: C_AB(t) = sum_p C_A(p,t) * C_B(-p,t)
 */

#ifndef __inline_multi_hadron_fh_h__
#define __inline_multi_hadron_fh_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma
{
  namespace InlineMultiHadronFHEnv
  {
    //! Registration
    bool registerAll();

    //! Parameters
    struct Params
    {
      //! Boosted momentum specification for nonzero total momentum
      struct BoostedMomentum_t
      {
        multi1d<int> total_momentum;  //!< Total momentum P (3-vector)
        multi1d<int> fiducial_mom1;   //!< Fiducial momentum for hadron 1
        multi1d<int> fiducial_mom2;   //!< Fiducial momentum for hadron 2
        std::string irrep;            //!< Irrep to project onto (e.g., "A1+E")
      };

      //! Multi-hadron state specification
      struct MultiHadronState_t
      {
        std::string name;             //!< Label for this state (e.g., "deuteron")
        std::string type;             //!< State type: DIBARYON, DIMESON
        std::string hadron1;          //!< First hadron type (PROTON, NEUTRON, PION, etc.)
        std::string hadron2;          //!< Second hadron type
        bool is_boosted;              //!< Whether this uses nonzero total momentum
        BoostedMomentum_t boost;      //!< Boosted momentum parameters (if is_boosted)
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

        //! List of multi-hadron states to compute
        multi1d<MultiHadronState_t> states;

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

  } // namespace InlineMultiHadronFHEnv

} // namespace Chroma

#endif
