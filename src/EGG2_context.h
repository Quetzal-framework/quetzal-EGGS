
// Copyright 2020 Arnaud Becheler    <arnaud.becheler@gmail.com>

/***********************************************************************                                                                         *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
************************************************************************/

#ifndef __EGG2_CONTEXT_H_INCLUDED__
#define __EGG2_CONTEXT_H_INCLUDED__

#include "EGG2_database.h"
#include "quetzal/quetzal.h"

#include "utils.h"
#include "rapidcsv.h"
#include "sqlite3pp.h"

#include <boost/program_options.hpp>

#include <random>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <functional>
#include <cmath>

namespace geo = quetzal::geography;
namespace demography = quetzal::demography;
namespace genet = quetzal::genetics;
namespace expr = quetzal::expressive;
namespace bpo = boost::program_options;

/// @brief Principal simulation class gathering simulation context elements.
/// @details
/// Gather context elements such as program options, database access, landscape representation,
/// sampled nodes, reproduction/disersal expressions and forward them to a simulation core.
class SimulationContext
{
private:
  // Type declaration
  using time_type = int;
  using landscape_type = geo::DiscreteLandscape<std::string,time_type>;
  using coord_type = landscape_type::coord_type;
  using sample_type = std::vector<EGGS::utils::GeneCopy>;
  using dispersal_policy = quetzal::demography::dispersal_policy::mass_based;
  using coalescence_policy = quetzal::coalescence::newick_with_distance_to_parent_and_leaf_name<coord_type, time_type>;
  using memory_policy = quetzal::memory::on_demand;
  using core_type = quetzal::ForwardBackwardSpatiallyExplicit<coord_type, dispersal_policy, coalescence_policy, memory_policy>;
  using options_type = bpo::variables_map;
  using dispersal_type = dispersal_policy::neighboring_migration
  <
    coord_type,
    std::function<double(coord_type)>,
    std::function<std::vector<coord_type>(coord_type)>
  >;
  using reproduction_type = std::function<unsigned int(std::mt19937&, coord_type, time_type)>;
  // Members
  bool verbose;
  bpo::variables_map m_vm;
  database_type m_database;
  landscape_type m_landscape;
  sample_type m_sample;
  core_type m_core;
  reproduction_type m_reproduction_expr;
  time_type m_duration;
  dispersal_type m_dispersal_kernel;
  std::string m_newicks;
public:
  SimulationContext(bpo::variables_map opts, std::mt19937& gen, bool verbose):
  verbose(verbose),
  m_vm(opts),
  m_database(build_database()),
  m_landscape(build_landscape()),
  m_sample(build_sample()),
  m_core(build_simulation_core()),
  m_reproduction_expr(build_reproduction_function(gen)),
  m_duration(m_vm["duration"].as<int>()),
  m_dispersal_kernel(build_dispersal_kernel())
  {
    if(verbose) show_reprojected_sample();
  };
  /// @brief Entry point after contruction
  void run(std::mt19937& gen){
    simulate_forward_demography(gen);
    maybe_save_demography();
    int nb_reuse = m_vm["reuse"].as<int>();
    for(int i=1; i <= nb_reuse; ++i)
    {
      try
      {
        simulate_coalescence(gen);
        record_params_and_genealogies();
      }
      catch(const std::domain_error &)
      {
        record_params_and_failure();
      }
    }
  }
private:
  database_type build_database()
  {
    if(verbose){std::cout << "Database initialization" << std::endl;}
    try
    {
      std::string filename = m_vm["output"].as<std::string>();
      return database_type(filename);
    }
    catch(const std::exception& e)
    {
      throw std::runtime_error(std::string("In SimulationContext: Error when building database. ") + e.what());
    }
  }

  landscape_type build_landscape()
  {
    std::string filename = m_vm["suitability"].as<std::string>();
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

  sample_type build_sample()
  {
    if(verbose){std::cout << "Tip nodes initialization" << std::endl;}
    std::string datafile = m_vm["tips"].as<std::string>();
    using EGGS::utils::GeneCopy;
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
    int N_0 = m_vm["N_0"].as<int>();
    int duration = m_vm["duration"].as<int>();
    core_type core(x_0, N_0, duration);
    core.ancestral_Wright_Fisher_N(N_0);
    return core;
  }

  reproduction_type build_reproduction_function(std::mt19937 & gen)
  {
    if(verbose){std::cout << "Reproduction expression initialization" << std::endl;}
    using expr::literal_factory;
    using expr::use;
    // growth rate:
    literal_factory<coord_type, time_type> lit;
    auto r = lit( m_vm["r"].as<double>() );
    // carrying capacity:
    auto suitability = m_landscape["suitability"];
    int K_suit = m_vm["K_suit"].as<int>();
    int K_min_ocean = m_vm["K_min_ocean"].as<int>();
    int K_max_ocean = m_vm["K_max_ocean"].as<int>();
    double p_K_ocean = m_vm["p_K_ocean"].as<double>();
    int K_min_matrix = m_vm["K_min_matrix"].as<int>();
    int K_max_matrix = m_vm["K_max_matrix"].as<int>();
    double p_K_matrix = m_vm["p_K_matrix"].as<double>();
    auto K = [K_suit, K_min_ocean, K_max_ocean, p_K_ocean, K_min_matrix, K_max_matrix, p_K_matrix, &gen, suitability](coord_type const& x, time_type t)
    {
      if ( suitability(x,0) < 0.0) // NA: ocean cell
      {
        return std::bernoulli_distribution(p_K_ocean)(gen) ? K_max_ocean : K_min_ocean;
      }
      else if (suitability(x,0) == 0.0) // 0: matrix cell
      {
        double value = sin(static_cast<double>(t)*2.3/100.0);
        if(value < 0.0) value = 0.0;
        return std::bernoulli_distribution(p_K_matrix)(gen) ? static_cast<int>(static_cast<double>(K_max_matrix)*value) : static_cast<int>(static_cast<double>(K_min_matrix)*value);
      }
      else // non-zero, non-na: sky islands
      {
        return static_cast<int>(static_cast<double>(K_suit)*suitability(x,0));
      }
    };
    // Retrieve population size reference to define a logistic growth process
    auto N = m_core.get_functor_N();
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
    double emigrant_rate = m_vm["emigrant_rate"].as<double>();
    auto suitability = m_landscape["suitability"];
    // std::function for knowing type to store as member
    std::function<double(coord_type)> friction = [suitability](coord_type const& x)
    {
      if(suitability(x,0) <= 0.0) {return 0.01;} //ocean cell or matrix
      else return 1.0 - suitability(x, 0);
    };
    auto env_ref = std::cref(m_landscape);
    // std::function for knowing type to store as member
    std::function<std::vector<coord_type>(coord_type)> get_neighbors = [env_ref](coord_type const& x)
    {
      return env_ref.get().direct_neighbors(x);
    };
    dispersal_type d =  dispersal_policy::make_neighboring_migration(coord_type(), emigrant_rate, friction, get_neighbors);
    return d;
  }

  void simulate_forward_demography(std::mt19937 & gen)
  {
    m_core.simulate_forward_demography(m_reproduction_expr, m_dispersal_kernel, gen);
  }

  void maybe_save_demography()
  {
    if(m_vm.count("log-history"))
    {
      std::string filename = m_vm["log-history"].as<std::string>();
      if(std::filesystem::exists(filename))
      {
        std::string message("Unable to save demography: file " +filename+ " already exists.");
        throw(std::runtime_error(message));
      }
      using expr::use;
      auto core_ref = std::cref(m_core);
      auto N = m_core.get_functor_N();
      auto space = [core_ref](time_type t){return core_ref.get().distribution_area(t);};
      m_landscape.export_to_geotiff(N, 0, m_duration - 1, space, filename);
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
      genealogies.append(m_core.coalesce_to_mrca<>(m_sample, m_duration, get_position, get_name, gen));
      genealogies.append("\n\n");
    }
    genealogies.pop_back();
    genealogies.pop_back();
    m_newicks = genealogies;
  }

  void record_params_and_genealogies()
  {
    auto ID = m_database.insert_params_results_and_get_rowid(m_vm, m_newicks);
    std::cout << ID << "\tRECORDED" <<std::endl;
  }

  void record_params_and_failure()
  {
    auto ID = m_database.insert_params_failure_and_get_rowid(m_vm);
    std::cerr << ID << "\tFAILED" <<std::endl;
  }
}; // end SimulationContext

#endif
