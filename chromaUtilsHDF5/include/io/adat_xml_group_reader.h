// -*- C++ -*-
/*! \file
 *  \brief Read an XML group as a string
 */

#ifndef __adat_xml_group_reader_h__
#define __adat_xml_group_reader_h__

#include "adat_xmlio.h"

namespace ADATXML
{

  //! Hold group xml and type id
  /*! \ingroup io */
  struct GroupXML_t
  {
    std::string  xml;     /*!< xml holding group */
    std::string  id;      /*!< typeid within group */
    std::string  path;    /*!< pathname of root */
  };


  //! Read group and return in a struct
  /*! \ingroup io */
  GroupXML_t readXMLGroup(XMLReader& xml, 
			  const std::string& path, const std::string& type_name);

  //! Read group and return in a struct
  /*! \ingroup io */
  Array<GroupXML_t> readXMLArrayGroup(XMLReader& xml, 
				      const std::string& path, 
				      const std::string& type_name);

  //! Read group and return in a struct
  /*! \ingroup io */
  std::vector<GroupXML_t> readXMLVectorGroup(XMLReader& xml, 
					     const std::string& path, 
					     const std::string& type_name);

} //end namespace ADATXML

#endif
