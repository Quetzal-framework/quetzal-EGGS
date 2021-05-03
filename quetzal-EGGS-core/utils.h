// Copyright 2020 Arnaud Becheler    <arnaud.becheler@gmail.com>

/***********************************************************************                                                                         *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
************************************************************************/

#ifndef __DECRYPT_UTILS_H_INCLUDED__
#define __DECRYPT_UTILS_H_INCLUDED__

#include "quetzal/quetzal.h"

namespace decrypt::utils{

  // used in dispersal kernel definition
  template<typename T>
  auto make_neighboring_cells_functor(T const& landscape){

    using coord_type = quetzal::geography::GeographicCoordinates;
    return [landscape](coord_type const& x0){
      auto res = landscape.get().resolution();

      std::vector<coord_type> v;

      coord_type x1(x0);
      x1.lon() += res.lon();
      if(landscape.get().is_in_spatial_extent(x1)){
        v.push_back(landscape.get().reproject_to_centroid(x1));
      }

      coord_type x2(x0);
      x2.lat() += res.lat();
      if(landscape.get().is_in_spatial_extent(x2)){
        v.push_back(landscape.get().reproject_to_centroid(x2));
      }

      coord_type x3(x0);
      x3.lon() -= res.lon();
      if(landscape.get().is_in_spatial_extent(x3)){
        v.push_back(landscape.get().reproject_to_centroid(x3));
      }

      coord_type x4(x0);
      x4.lat() -= res.lat() ;
      if(landscape.get().is_in_spatial_extent(x4)){
        v.push_back(landscape.get().reproject_to_centroid(x4));
      }
      return v;
    };
  };

  // Class for representing sampled individual gene copies
  class GeneCopy
  {
  private:
    using coord_type = quetzal::geography::GeographicCoordinates;
    std::string m_id;
    coord_type m_x;
  public:
    GeneCopy(std::string const& id, coord_type const& x): m_id(id), m_x(x){}
    GeneCopy(coord_type const& x): m_id(""), m_x(x){}
    auto id() const {return m_id;}
    auto x() const {return m_x;}
  };
}

#include <string.h>
#include <iostream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/any.hpp"

namespace po = boost::program_options;

inline void PrintVariableMap(const boost::program_options::variables_map vm) {
    for (po::variables_map::const_iterator it = vm.begin(); it != vm.end(); it++) {
        std::cout << "> " << it->first;
        if (((boost::any)it->second.value()).empty()) {
            std::cout << "(empty)";
        }
        if (vm[it->first].defaulted() || it->second.defaulted()) {
            std::cout << "(default)";
        }
        std::cout << "=";

        bool is_char;
        try {
            boost::any_cast<const char *>(it->second.value());
            is_char = true;
        } catch (const boost::bad_any_cast &) {
            is_char = false;
        }
        bool is_str;
        try {
            boost::any_cast<std::string>(it->second.value());
            is_str = true;
        } catch (const boost::bad_any_cast &) {
            is_str = false;
        }

        if (((boost::any)it->second.value()).type() == typeid(int)) {
            std::cout << vm[it->first].as<int>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(bool)) {
            std::cout << vm[it->first].as<bool>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(double)) {
            std::cout << vm[it->first].as<double>() << std::endl;
        } else if (is_char) {
            std::cout << vm[it->first].as<const char * >() << std::endl;
        } else if (is_str) {
            std::string temp = vm[it->first].as<std::string>();
            if (temp.size()) {
                std::cout << temp << std::endl;
            } else {
                std::cout << "true" << std::endl;
            }
        } else { // Assumes that the only remainder is vector<string>
            try {
                std::vector<std::string> vect = vm[it->first].as<std::vector<std::string> >();
                uint i = 0;
                for (std::vector<std::string>::iterator oit=vect.begin();
                     oit != vect.end(); oit++, ++i) {
                    std::cout << "\r> " << it->first << "[" << i << "]=" << (*oit) << std::endl;
                }
            } catch (const boost::bad_any_cast &) {
                std::cout << "UnknownType(" << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
            }
        }
    }
}

#endif
