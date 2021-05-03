
// Copyright 2020 Arnaud Becheler    <arnaud.becheler@gmail.com>

/***********************************************************************                                                                         *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
************************************************************************/

#ifndef __M4_REALITY_CHECK_H_INCLUDED__
#define __M4_REALITY_CHECK_H_INCLUDED__

#include "quetzal/include/quetzal/geography/DiscreteLandscape.h"

#include "utils.h"
#include "rapidcsv.h"

#include <boost/program_options.hpp>

#include "sqlite3pp.h"

#include <random>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <functional>

namespace coal = quetzal::coalescence;
namespace geo = quetzal::geography;
namespace demography = quetzal::demography;
namespace genet = quetzal::genetics;
namespace sim = quetzal::simulator;
namespace expr = quetzal::expressive;
namespace bpo = boost::program_options;


// Returns a map with the program options
auto handle_options(int argc, char* argv[])
{
  bpo::variables_map vm;
  try
  {
    bpo::options_description generalOptions{"General"};
    generalOptions.add_options()
    ("help,h", "Help screen")
    ("v", "verbose mode")
    ("config", bpo::value<std::string>()->required(), "Config file")
    ("landscape", bpo::value<std::string>()->required(), "Geospatial file in tiff format giving the friction map")
    ("reuse", bpo::value<int>()->required(), "How many times pseudo-observed should be generated under the same spatial history")
    ("n_loci", bpo::value<int>()->required(), "Number of loci to simulate")
    ("sample",  bpo::value<std::string>()->required(), "File name for the lon/lat of sampled genetic material")
    ("lon_0", bpo::value<double>()->required(), "Introduction point longitude")
    ("lat_0", bpo::value<double>()->required(), "Introduction point latitude")
    ("N_0", bpo::value<int>()->required(), "Number of gene copies at introduction point")
    ("duration", bpo::value<int>()->required(), "Number of generations to simulate")
    ("K_suit", bpo::value<int>()->required(), "Carrying capacity in suitable areas")
    ("K_max", bpo::value<int>()->required(), "Highest carrying capacity in areas with null suitability")
    ("K_min", bpo::value<int>()->required(), "Lowest carrying capacity in areas with null suitability")
    ("p_K", bpo::value<double>()->required(), "Probability to have highest carrying capacity in areas with 0 suitability")
    ("r", bpo::value<double>()->required(), "Growth rate")
    ("emigrant_rate", bpo::value<double>()->required(), "Emigrant rate between the four neighboring cells")
    ("demography_out",  bpo::value<std::string>(), "File name for the simulated demography output")
    ("database", bpo::value<std::string>()->required(), "Filename database storing the output");

    bpo::options_description fileOptions{"File"};
    fileOptions.add_options()
    ("landscape", bpo::value<std::string>()->required(), "Geospatial file in tiff format giving the friction map")
    ("reuse", bpo::value<int>()->required(), "How many times pseudo-observed should be generated under the same spatial history")
    ("n_loci", bpo::value<int>()->required(), "Number of loci to simulate")
    ("sample",  bpo::value<std::string>()->required(), "File name for the lon/lat of sampled genetic material")
    ("lon_0", bpo::value<double>()->required(), "Introduction point longitude")
    ("lat_0", bpo::value<double>()->required(), "Introduction point latitude")
    ("N_0", bpo::value<int>()->required(), "Number of gene copies at introduction point")
    ("duration", bpo::value<int>()->required(), "Number of generations to simulate")
    ("K_suit", bpo::value<int>()->required(), "Carrying capacity in suitable areas")
    ("K_max", bpo::value<int>()->required(), "Highest carrying capacity in areas with null suitability")
    ("K_min", bpo::value<int>()->required(), "Lowest carrying capacity in areas with null suitability")
    ("p_K", bpo::value<double>()->required(), "Probability to have highest carrying capacity in areas with 0 suitability")
    ("r", bpo::value<double>()->required(), "Growth rate")
    ("emigrant_rate", bpo::value<double>()->required(), "Emigrant rate between the four neighboring cells")
    ("demography_out",  bpo::value<std::string>(), "File name for the simulated demography output")
    ("database", bpo::value<std::string>()->required(), "Filename database storing the output");

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
} // end of handle_options

//Class for wrapping sqlite3pp code
class database_type
{
public:

  database_type(std::string const& filename){
  this->m_database = sqlite3pp::database(filename.c_str());
  create_table();
  }

  auto insert_params_results_and_get_rowid(bpo::variables_map vm, std::string const& newicks)
  {
    sqlite3pp::command cmd(
      this->m_database,
      "INSERT INTO core4_pods (lon_0, lat_0, N_0, duration, K_suit, K_max, K_min, p_K, r, emigrant_rate, newicks) VALUES (?,?,?,?,?,?,?,?,?,?,?)"
    );
    cmd.binder() << vm["lon_0"].as<double>()
                  << vm["lat_0"].as<double>()
                  << vm["N_0"].as<int>()
                  << vm["duration"].as<int>()
                  << vm["K_suit"].as<int>()
                  << vm["K_max"].as<int>()
                  << vm["K_min"].as<int>()
                  << vm["p_K"].as<double>()
                  << vm["r"].as<double>()
                  << vm["emigrant_rate"].as<double>()
                  << newicks;
    cmd.execute();
    return this->m_database.last_insert_rowid();
  }

private:
  sqlite3pp::database m_database;

  void create_table()
  {
    sqlite3pp::command cmd(
      this->m_database,
      "CREATE TABLE IF NOT EXISTS core4_pods(lon_0 DOUBLE, lat_0 DOUBLE, N_0 INTEGER, duration INTEGER, K_suit INTEGER, K_max INTEGER, K_min INTEGER, p_K DOUBLE, r DOUBLE, emigrant_rate DOUBLE, newicks TEXT)"
    );
    cmd.execute();
  }

};

class SimulationContext
{
public:
  SimulationContext(bpo::variables_map opts, std::mt19937& gen, bool verbose):
  verbose(verbose),
  m_vm(opts),
  m_database(build_database()),
  m_landscape(build_landscape()),
  m_sample(build_sample()),
  m_core(build_simulation_core()),
  m_reproduction_expr(build_reproduction_function(gen)),
  m_dispersal_kernel(build_dispersal_kernel())
  {
    if(verbose){
      show_reprojected_sample();
    }
  };

  void run(std::mt19937& gen)
  {
    expand_demography(gen);
    maybe_save_demography();
    int nb_reuse = m_vm["reuse"].as<int>();
    for(int i=1; i <= nb_reuse; ++i){
      simulate_coalescence(gen);
      save_genealogies();
    }
  }

private:

  using time_type = int;
  using landscape_type = geo::DiscreteLandscape<std::string,time_type>;
  using coord_type = landscape_type::coord_type;
  using sample_type = std::vector<decrypt::utils::GeneCopy>;
  using demographic_policy = demography::strategy::mass_based;
  using coal_policy = coal::policies::distance_to_parent_leaf_name<coord_type, time_type>;
  using core_type = sim::SpatiallyExplicit<coord_type, time_type, demographic_policy, coal_policy>;
  using options_type = bpo::variables_map;
  // eww, but it works
  using dispersal_type = demography::strategy::mass_based::light_neighboring_migration
  <
    coord_type,
    std::function<double(coord_type)>,
    std::function<std::vector<coord_type>(coord_type)>
  >;
  using reproduction_type = std::function<unsigned int(std::mt19937&, coord_type, time_type)>;

  bool verbose;
  bpo::variables_map m_vm;
  database_type m_database;
  landscape_type m_landscape;
  sample_type m_sample;
  core_type m_core;
  reproduction_type m_reproduction_expr;
  time_type m_t_0;
  time_type m_sample_time;
  dispersal_type m_dispersal_kernel;
  std::string m_newicks;

  database_type build_database()
  {
    if(verbose){std::cout << "Database initialization" << std::endl;}
    try
    {
      std::string filename = m_vm["database"].as<std::string>();
      return database_type(filename);
    }
    catch(const std::exception& e)
    {
      throw std::runtime_error(std::string("In SimulationContext: Error when building database. ") + e.what());
    }
  }

  // TODO: pas sûr de la syntaxe du try/catch
  landscape_type build_landscape()
  {
    std::string filename = m_vm["landscape"].as<std::string>();
    if(verbose){std::cout << "Landscape initialization" << std::endl;}
    try
    {
      return landscape_type({{"suitability", filename}}, {time_type(0)});
    }
    catch(const std::exception& e)
    {
      throw std::runtime_error(std::string("In SimulationContext: Error when building landscape. ") + e.what());
    }
  }

  // TODO haploid versus diploid. Plus, format changed: ID/coordinates/no genetics
  sample_type build_sample()
  {
    if(verbose){std::cout << "Sample initialization" << std::endl;}
    std::string datafile = m_vm["sample"].as<std::string>();
    using decrypt::utils::GeneCopy;
    rapidcsv::Document doc(datafile);
    std::vector<std::string> IDs = doc.GetColumn<std::string>("ID");
    std::vector<double> v_lat = doc.GetColumn<double>("latitude");
    std::vector<double> v_lon = doc.GetColumn<double>("longitude");
    unsigned int i = 0;
    std::vector<GeneCopy> gene_copies;
    for(auto const& id : IDs)
    {
      coord_type x;
      x.lat(v_lat.at(i));
      x.lon(v_lon.at(i));
      x = m_landscape.reproject_to_centroid(x);
      gene_copies.push_back(GeneCopy(id, x));
      ++i;
    }
    return gene_copies;
  }

  void show_reprojected_sample()
  {
    std::cout << "Reprojected sample:\n\n";
    for(auto const& gene_copy: m_sample)
    {
      std::cout << gene_copy.id() << "\t" << gene_copy.x() << std::endl;
    }
  }

  core_type build_simulation_core()
  {
    if(verbose){std::cout << "Simulation core initialization" << std::endl;}
    coord_type x_0(m_vm["lat_0"].as<double>(), m_vm["lon_0"].as<double>());
    x_0 = m_landscape.reproject_to_centroid(x_0);
    m_t_0 = 0;
    m_sample_time = m_vm["duration"].as<int>();
    int N_0 = m_vm["N_0"].as<int>();
    core_type core(x_0, m_t_0, N_0);
    core.ancestral_Wright_Fisher_N(N_0);
    return core;
  }

  reproduction_type build_reproduction_function(std::mt19937 & gen)
  {
    if(verbose){std::cout << "Reproduction expression initialization" << std::endl;}
    using expr::literal_factory;
    using expr::use;

    // growth rate
    literal_factory<coord_type, time_type> lit;
    auto r = lit( m_vm["r"].as<double>() );

    // carrying capacity
    auto suitability = m_landscape["suitability"];
    int K_suit = m_vm["K_suit"].as<int>();
    int K_min = m_vm["K_min"].as<int>();
    int K_max = m_vm["K_max"].as<int>();
    double p_K = m_vm["p_K"].as<double>();
    auto K = [K_suit, K_min, K_max, p_K, &gen, suitability](coord_type const& x, time_type)
    {
      if( suitability(x,0) <= 0.1)
      { //ocean cell
        return std::bernoulli_distribution(p_K)(gen) ? K_max : K_min;
      }else{
        return static_cast<int>(static_cast<double>(K_suit)*suitability(x,0));
      }
    };

    // Retrieve population size reference to define a logistic growth process
    auto pop_sizes = m_core.pop_size_history();
    auto N = use( [pop_sizes](coord_type x, time_type t){ return pop_sizes(x,t);} );
    auto g = N * ( lit(1) + r ) / ( lit(1) + ( (r * N)/K ));
    auto reproduction = [g](auto& gen, coord_type const&x, time_type t){
      std::poisson_distribution<unsigned int> poisson(g(x,t));
      return poisson(gen);
    };
    return reproduction;
  }

  dispersal_type build_dispersal_kernel()
  {
    if(verbose){std::cout << "Dispersal kernel initialization" << std::endl;}
    auto suitability = m_landscape["suitability"];
    std::function<double(coord_type)> friction = [suitability](coord_type const& x){
      if(suitability(x,0) <= 0.5) {return 0.99;} //ocean cell
      else return 1.0 - suitability(x, 0);
    };
    double emigrant_rate = m_vm["emigrant_rate"].as<double>();
    auto env_ref = std::cref(m_landscape);
    std::function<std::vector<coord_type>(coord_type)> get_neighbors = decrypt::utils::make_neighboring_cells_functor(env_ref);
    return demographic_policy::make_light_neighboring_migration(coord_type(), emigrant_rate, friction, get_neighbors);
  }

  void expand_demography(std::mt19937 & gen)
  {
    m_core.expand_demography(m_sample_time, m_reproduction_expr, m_dispersal_kernel, gen);
  }

  void maybe_save_demography()
  {
    if(m_vm.count("demography_out"))
    {
      std::string filename = m_vm["demography_out"].as<std::string>();
      if(std::filesystem::exists(filename))
      {
        std::string message("Unable to save demography: file " +filename+ " already exists.");
        throw(std::runtime_error(message));
      }
      using expr::use;
      auto pop_sizes = m_core.pop_size_history();
      auto N = use( [pop_sizes](coord_type x, time_type t){ return pop_sizes(x,t);} );
      m_landscape.export_to_geotiff(N, m_t_0, m_sample_time, [&pop_sizes](time_type const& t){return pop_sizes.get().definition_space(t);}, filename);
    }
  }

  void simulate_coalescence(std::mt19937 & gen)
  {
    auto get_name = [](auto const& gene_copy, time_type){return gene_copy.id();};
    auto get_position = [](auto const& gene_copy, time_type){return gene_copy.x();};
    std::string genealogies;
    int n_loci = m_vm["n_loci"].as<int>();
    for(unsigned int locus = 0; locus < n_loci ; ++locus)
    {
      genealogies.append(m_core.coalesce_to_mrca<>(m_sample, m_sample_time, get_position, get_name, gen));
      genealogies.append("\n\n");
    }
    genealogies.pop_back();
    genealogies.pop_back();
    m_newicks = genealogies;
  }

  void save_genealogies()
  {
    auto ID = m_database.insert_params_results_and_get_rowid(m_vm, m_newicks);
    std::cout << ID << std::endl;
  }
};

#endif
