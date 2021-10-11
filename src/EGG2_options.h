
// Copyright 2020 Arnaud Becheler    <arnaud.becheler@gmail.com>

/***********************************************************************                                                                         *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
************************************************************************/

#ifndef __EGG2_OPTIONS_H_INCLUDED__
#define __EGG2_OPTIONS_H_INCLUDED__

#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>
namespace bpo = boost::program_options;

/// @brief Returns a map with the program options
auto handle_options(int argc, char* argv[])
{
  // Declare a group of options that will be allowed only on command line
  bpo::options_description generic_options{"Command-line-only options"};
  generic_options.add_options()
  ("help,h", "help screen")
  ("verbose,v", "verbose mode")
  ("version", "software version")
  ("config", bpo::value<std::string>()->required(), "configuration file")
  ;
  // Allowed both on command line and in config file
  bpo::options_description general_options{"Input/output options"};
  general_options.add_options()
  ("suitability", bpo::value<std::string>()->required(), "input suitability map in GeoTiFF (.tiff) format")
  ("tips",  bpo::value<std::string>()->required(), "input CSV file listing <ID,latitude,longitude> for every node")
  ("output", bpo::value<std::string>()->default_value("out.db"), "SQLite 3 database output")
  ;
  // Declare a simpler way to call on command line
  bpo::positional_options_description positional_options;
  positional_options.add("config", 1);
  positional_options.add("suitability", 1);
  positional_options.add("tips", 1);
  positional_options.add("output", 1);
  // Allowed both on command line and in config file
  bpo::options_description model_options{"Demogenetic model parameters"};
  model_options.add_options()
  ("n_loci", bpo::value<int>(), "number of loci")
  ("lon_0", bpo::value<double>(), "Origin point longitude")
  ("lat_0", bpo::value<double>(), "Origin point latitude")
  ("N_0", bpo::value<int>(), "Number of gene copies at introduction point")
  ("duration", bpo::value<int>(), "Number of generations to simulate")
  ("K_suit", bpo::value<int>(), "Carrying capacity in suitable areas")
  ("K_max_ocean", bpo::value<int>(), "Highest carrying capacity in areas with NA suitability")
  ("K_min_ocean", bpo::value<int>(), "Lowest carrying capacity in areas with NA suitability")
  ("p_K_ocean", bpo::value<double>(), "Probability to have highest carrying capacity in areas with NA suitability")
  ("K_max_matrix", bpo::value<int>(), "Highest carrying capacity in areas with 0 suitability")
  ("K_min_matrix", bpo::value<int>(), "Lowest carrying capacity in areas with 0 suitability")
  ("p_K_matrix", bpo::value<double>(), "Probability to have highest carrying capacity in areas with 0 suitability")
  ("r", bpo::value<double>(), "Growth rate")
  ("emigrant_rate", bpo::value<double>(), "Emigrant rate between the four neighboring cells")
  ;
  // Allowed both on command line and in config file
  bpo::options_description other_options{"Other options"};
  other_options.add_options()
  ("reuse", bpo::value<int>()->default_value(1), "number of pseudo-observed data to be simulated under one demographic history")
  ("log-history",  bpo::value<std::string>(), "output history in GeoTiff format if option specified")
  ;
  bpo::options_description command_line_options;
  command_line_options.add(generic_options).add(general_options).add(model_options).add(other_options);

  bpo::options_description file_options{"General options (command line values will overwrite congif file values)"};
  file_options.add(general_options).add(model_options).add(other_options);

  bpo::variables_map vm;
  try
  {
    bpo::store(
      bpo::command_line_parser(argc, argv).
      options(command_line_options).
      positional(positional_options).
      run(), vm); // can throw
    // --help option
    if (vm.count("help"))
    {
      std::cout << "--------------------------------------------------------------------------------------" << std::endl;
      std::cout << "| This is Quetzal-EGG-1 coalescence simulator.                                        |" << std::endl;
      std::cout << "|   - Purpose: simulate gene trees in an heterogeneous landscape.                     |" << std::endl;
      std::cout << "|   - Author: Arnaud Becheler, 2021.                                                  |" << std::endl;
      std::cout << "|   - Usage: " << argv[0] << " [options] <config> <suitability> <tips> <output> ...    " << std::endl;
      std::cout << "--------------------------------------------------------------------------------------|" << std::endl;
      std::cout << "\n" << generic_options << std::endl;
      std::cout << "\n" << file_options << std::endl;
      return vm;
      // SUCCESS
    }
    // --version option
    if (vm.count("version"))
    {
      std::cout << "quetzal-EGG-1 version 0.1" << std::endl;
      return vm;
      // SUCCESS
    }
    bpo::notify(vm); // throws on error, so do after help in case there are any problems
  }
  catch(boost::program_options::required_option& e)
  {
    throw;
  }
  catch(boost::program_options::error& e)
  {
    throw;
  }
  if (vm.count("config"))
  {
    std::ifstream ifs{vm["config"].as<std::string>().c_str()};
    if (ifs){
      store(parse_config_file(ifs, file_options), vm);
    }
  }
  notify(vm);
  return vm;
} // end of handle_options

#endif
