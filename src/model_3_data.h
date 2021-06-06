
// Copyright 2020 Arnaud Becheler    <arnaud.becheler@gmail.com>

/***********************************************************************                                                                         *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
************************************************************************/

#ifndef __M3_REALITY_CHECK_H_INCLUDED__
#define __M3_REALITY_CHECK_H_INCLUDED__

#include "include/quetzal.h"

#include <boost/program_options.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "sqlite3pp.h"

#include <random>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

namespace coal = quetzal::coalescence;
namespace geo = quetzal::geography;
namespace demo = quetzal::demography;
namespace sim = quetzal::simulator;
namespace schemes = quetzal::sampling_scheme;
namespace bpo = boost::program_options;

/*
* Return a map with the program options
*/
auto handle_options(int argc, char* argv[]){
  bpo::variables_map vm;
  try
  {
    bpo::options_description generalOptions{"General"};
    generalOptions.add_options()
    ("help,h", "Help screen")
    ("config", bpo::value<std::string>(), "Config file")
    ("landscape", bpo::value<std::string>()->required(), "Geospatial file in tiff format giving the friction map");

    bpo::options_description fileOptions{"File"};
    fileOptions.add_options()
    ("help", "produce help message")
    ("n_loci", bpo::value<unsigned int>()->required(), "Number of loci to simulate")
    ("sample",  bpo::value<std::string>(), "File name for the lon/lat of sampled genetic material")
    ("lon_0", bpo::value<double>()->required(), "Introduction point longitude")
    ("lat_0", bpo::value<double>()->required(), "Introduction point latitude")
    ("N_0", bpo::value<unsigned int>()->required(), "Number of gene copies at introduction point")
    ("duration", bpo::value<unsigned int>()->required(), "Number of generations to simulate")
    ("K", bpo::value<unsigned int>()->required(), "Carrying capacity")
    ("r", bpo::value<double>()->required(), "Growth rate")
    ("emigrant_rate", bpo::value<double>()->required(), "Emigrant rate between the four neighboring cells")
    ("demography_out",  bpo::value<std::string>(), "File name for the simulated demography output")
    ("database", bpo::value<std::string>(), "Filename database storing the output");

    store(parse_command_line(argc, argv, generalOptions), vm);
    if (vm.count("config"))
    {
      std::ifstream ifs{vm["config"].as<std::string>().c_str()};
      if (ifs){
        store(parse_config_file(ifs, fileOptions), vm);
      }
    }
    notify(vm);
    if (vm.count("help"))
    {
      std::cout << generalOptions << '\n';
    }
  }
  catch (const bpo::error &ex)
  {
    std::cerr << ex.what() << '\n';
  }
  return vm;
}

/*
 * Class for representing sampled individual gene copies
 */
class Ind {
private:
    using coord_type = geo::GeographicCoordinates;
    unsigned int m_id;
    coord_type m_x;
    static unsigned int m_next_available_id;
    auto next_id(){
      auto id = m_next_available_id;
      ++m_next_available_id;
      return id;
    }
public:
    Ind(coord_type const& x):
    m_id(next_id()),
    m_x(x)
    {}
    auto id() const {return m_id;}
    auto x() const {return m_x;}
};

/*
 * Initialization of the first ID
 */
unsigned int Ind::m_next_available_id = 0;

class SimulationContext{
public:
  static auto run(bpo::variables_map const& vm)
  {
    /******************************
    * Geospatial dataset
    *****************************/
    using time_type = int;
    using landscape_type = quetzal::geography::DiscreteLandscape<std::string,time_type>;
    using coord_type = landscape_type::coord_type;
    const std::string file = vm["landscape"].as<std::string>();
    landscape_type env( {{"friction", file}}, {time_type(0)} );
    std::cout << "Landscape initialized" << std::endl;

    /******************************
    * Genetic dataset
    *****************************/
    std::string datafile = vm["sample"].as<std::string>();
  	quetzal::genetics::Loader<coord_type, quetzal::genetics::microsatellite> reader;
  	auto dataset = reader.read(datafile);
    std::cout << "Dataset initialized:\n\n" << dataset << std::endl;
    dataset.reproject(env);
    std::cout << "Dataset reprojected:\n\n" << dataset << std::endl;

    /******************************
    * Simulator configuration
    *****************************/
    using demographic_policy = quetzal::demography::strategy::mass_based;
    using coalescence_policy = coal::policies::distance_to_parent_leaf_name<coord_type, time_type>;
    using simulator_type = sim::SpatiallyExplicit<coord_type, time_type, demographic_policy, coalescence_policy>;

    coord_type x_0(vm["lat_0"].as<double>(), vm["lon_0"].as<double>());
    x_0 = env.reproject_to_centroid(x_0);
    time_type t_0 = 0;
    time_type sample_time = vm["duration"].as<unsigned int>();
    unsigned int N_0 = vm["N_0"].as<unsigned int>();
    std::random_device rd;
    std::mt19937 gen(rd());

    // Initialize the simulator
    simulator_type simulator(x_0, t_0, N_0);
    std::cout << "Simulator initialized" << std::endl;

    /******************************
    * Niche and growth functions
    *****************************/
    using quetzal::expressive::literal_factory;
    using quetzal::expressive::use;
    literal_factory<coord_type, time_type> lit;

    auto K = lit( vm["K"].as<unsigned int>() );
    auto r = lit( vm["r"].as<double>() );

    // Retrieve population size reference to define a logistic growth process
    auto pop_sizes = simulator.pop_size_history();
    auto N = use( [pop_sizes](coord_type x, time_type t){ return pop_sizes(x,t);} );

    auto g = N * ( lit(1) + r ) / ( lit(1) + ( (r * N)/K ));
    auto sim_children = [g](auto& gen, coord_type const&x, time_type t){
      std::poisson_distribution<unsigned int> poisson(g(x,t));
      return poisson(gen);
    };

    /******************
    * Dispersal
    ******************/
    auto f = env["friction"];
    auto friction = [f](coord_type x){return f(x, 0);};
    double emigrant_rate = vm["emigrant_rate"].as<double>();
    auto env_ref = std::cref(env);
    auto get_neighbors = [env_ref](coord_type const& x){return env_ref.get().four_nearest_defined_cells(x);};
    auto const& space = env.geographic_definition_space();
    auto dispersal = demographic_policy::make_neighboring_migration(space, emigrant_rate, friction, get_neighbors);

    auto distance = [](coord_type const& a, coord_type const& b) -> double {return memoized_great_circle_distance(a,b);};
    std::cout << "Dispersal initialized" << std::endl;

    /**************************
    * Demographic expansion
    **************************/
    std::cout << "--- Expanding demography" << std::endl;
    simulator.expand_demography(sample_time, sim_children, dispersal, gen);

    // Visualization
    if(vm.count("demography_out"))
    {
      std::string filename = vm["demography_out"].as<std::string>();
      env.export_to_geotiff(N, t_0, sample_time, [&pop_sizes](time_type const& t){return pop_sizes.get().definition_space(t);}, filename);
    }

    /**********************
    *  Coalescence
    **********************/


    std::cout << "--- Simulating coalescents" << std::endl;
    std::vector<Ind> v;

    for(auto const& it1 : dataset.get_sampling_points())
    {
      for(unsigned int i = 0; i < dataset.individuals_at(it1).size(); ++ i)
      {
        v.emplace_back(it1);
      }
    }

    auto get_name = [](auto const& ind, time_type){return std::to_string(ind.id());};
    auto get_position = [](auto const& ind, time_type){return ind.x();};

    std::string genealogies;
    unsigned int n_loci = vm["n_loci"].as<unsigned int>();
    for(unsigned int locus = 0; locus < n_loci ; ++locus){
      genealogies.append(simulator.coalesce_to_mrca<>(v, sample_time, get_position, get_name, gen));
      genealogies.append("\n\n");
    }
    genealogies.pop_back();
    genealogies.pop_back();

    if(vm.count("database")){
      try{
        std::string file = vm["database"].as<std::string>();
        sqlite3pp::database db(file.c_str());
        sqlite3pp::command cmd(db, "CREATE TABLE IF NOT EXISTS results(id INTEGER PRIMARY KEY AUTOINCREMENT UNIQUE, genealogies TEXT)");
        cmd.execute();

        sqlite3pp::command cmd2(db, "INSERT INTO results (genealogies) VALUES (?)");
        cmd2.binder() << genealogies;
        cmd2.execute();
      }
      catch(const std::exception& e){
        std::cout << e.what() << std::endl;
      }
    }else{
      std::cout    << "--- Genealogies in Newick format:\n\n"
                   << genealogies
                   << std::endl;
    }
    return simulator;

  }
};

#endif
