
// Copyright 2020 Arnaud Becheler    <arnaud.becheler@gmail.com>

/***********************************************************************                                                                         *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
************************************************************************/

#ifndef __EGG2_DATABASE_H_INCLUDED__
#define __EGG2_DATABASE_H_INCLUDED__

#include "sqlite3pp.h"
#include <boost/program_options.hpp>
#include <string>

namespace bpo = boost::program_options;

/// @brief Class for wrapping sqlite3pp code
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
      "INSERT INTO quetzal_EGG_2 (lon_0, lat_0, N_0, duration, K_suit, K_max_ocean, K_min_ocean, p_K_ocean, K_max_matrix, K_min_matrix, p_K_matrix, r, emigrant_rate, newicks) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    );
    cmd.binder() << vm["lon_0"].as<double>()
                  << vm["lat_0"].as<double>()
                  << vm["N_0"].as<int>()
                  << vm["duration"].as<int>()
                  << vm["K_suit"].as<int>()
                  << vm["K_max_ocean"].as<int>()
                  << vm["K_min_ocean"].as<int>()
                  << vm["p_K_ocean"].as<double>()
                  << vm["K_max_matrix"].as<int>()
                  << vm["K_min_matrix"].as<int>()
                  << vm["p_K_matrix"].as<double>()
                  << vm["r"].as<double>()
                  << vm["emigrant_rate"].as<double>()
                  << newicks;
    cmd.execute();
    return this->m_database.last_insert_rowid();
  }

  auto insert_params_failure_and_get_rowid(bpo::variables_map vm)
  {
    sqlite3pp::command cmd(
      this->m_database,
      "INSERT INTO quetzal_EGG_2 (lon_0, lat_0, N_0, duration, K_suit, K_max_ocean, K_min_ocean, p_K_ocean, K_max_matrix, K_min_matrix, p_K_matrix, r, emigrant_rate, newicks) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    );
    cmd.binder() << vm["lon_0"].as<double>()
                  << vm["lat_0"].as<double>()
                  << vm["N_0"].as<int>()
                  << vm["duration"].as<int>()
                  << vm["K_suit"].as<int>()
                  << vm["K_max_ocean"].as<int>()
                  << vm["K_min_ocean"].as<int>()
                  << vm["p_K_ocean"].as<double>()
                  << vm["K_max_matrix"].as<int>()
                  << vm["K_min_matrix"].as<int>()
                  << vm["p_K_matrix"].as<double>()
                  << vm["r"].as<double>()
                  << vm["emigrant_rate"].as<double>()
                  << "";
    cmd.execute();
    return this->m_database.last_insert_rowid();
  }

private:
  sqlite3pp::database m_database;

  void create_table()
  {
    sqlite3pp::command cmd(
      this->m_database,
      "CREATE TABLE IF NOT EXISTS quetzal_EGG_2(lon_0 DOUBLE, lat_0 DOUBLE, N_0 INTEGER, duration INTEGER, K_suit INTEGER, K_max_ocean INTEGER, K_min_ocean INTEGER, p_K_ocean DOUBLE, K_max_matrix INTEGER, K_min_matrix INTEGER, p_K_matrix DOUBLE, r DOUBLE, emigrant_rate DOUBLE, newicks TEXT)"
    );
    cmd.execute();
  }

};

#endif
