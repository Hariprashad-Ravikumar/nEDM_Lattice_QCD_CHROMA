/*! \file
 *  \brief Inline Weinberg three-gluon operator term
 *
 */

#include "meas/inline/glue/inline_ggg.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/ggg.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"


namespace Chroma 
{ 

  namespace InlineGGGEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
          const std::string& path)
      {
        Params p(xml_in, path);
        return new InlineMeas(p);
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "WEINBERG_GGG";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
        success &= CreateGaugeStateEnv::registerAll();
        success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
        registered = true;
      }
      return success;
    }


    //! GGG input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
        case 1:
          if (paramtop.count("GaugeState") != 0)
            param.cgs = readXMLGroup(paramtop, "GaugeState", "Name");
          else
            param.cgs = CreateGaugeStateEnv::nullXMLGroup();
          break;

        default:
          QDPIO::cerr << "Params::Param_t: " << version 
            << " unsupported." << std::endl;
          QDP_abort(1);
      }

    }

    //! GGG output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      xml << param.cgs.xml;

      pop(xml);
    }


    //! GGG input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
    }

    //! GGG output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);

      pop(xml);
    }


    // Params
    Params::Params()
    { 
      frequency = 0; 
      param.cgs = CreateGaugeStateEnv::nullXMLGroup();
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

        // Params
        read(paramtop, "Param", param);

        // Ids
        read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
        QDPIO::cerr << "Caught Exception reading XML: " << e << std::endl;
        QDP_abort(1);
      }
    }

    void 
      InlineMeas::operator()(unsigned long update_no,
          XMLWriter& xml_out) 
      {
        START_CODE();

        StopWatch snoop;
        snoop.reset();
        snoop.start();

        // Grab the object
        const multi1d<LatticeColorMatrix>& u = 
          TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

        push(xml_out, "Weinberg_GGG");
        write(xml_out, "update_no", update_no);

        multi1d<DComplex> ggg_ts, qtop_ts;
        DComplex ggg, qtop;
        weinberg_ggg(u, ggg_ts, ggg, qtop_ts, qtop);

        write(xml_out, "ggg", real(ggg));
        write(xml_out, "qtop", real(qtop));
        write(xml_out, "ggg_ts", ggg_ts);
        write(xml_out, "qtop_ts", qtop_ts);

        pop(xml_out); // pop("Weinberg_GGG");

      snoop.stop();
      QDPIO::cout << InlineGGGEnv::name << ": total time = "
      << snoop.getTimeInSeconds()
      << " secs" << std::endl;
        END_CODE();
      } 

  }
}
