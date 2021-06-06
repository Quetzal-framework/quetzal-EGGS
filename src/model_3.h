#include "include/quetzal.h"
#include <boost/program_options.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/progress.hpp>

#include "sqlite3pp.h"

#include <random>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string>

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
    ("landscape", bpo::value<std::string>()->required(), "Geospatial file in tiff format");

    bpo::options_description fileOptions{"File"};
    fileOptions.add_options()
    ("help", "produce help message")
    ("n_sim_gen", bpo::value<unsigned int>()->required(), "Number of genetic simulations")
    ("n_loci", bpo::value<unsigned int>()->required(), "Number of loci to simulate")
    ("lon_0", bpo::value<double>()->required(), "Introduction point longitude")
    ("lat_0", bpo::value<double>()->required(), "Introduction point latitude")
    ("N_0", bpo::value<unsigned int>()->required(), "Number of gene copies at introduction point")
    ("duration", bpo::value<unsigned int>()->required(), "Number of generations to simulate")
    ("lon_1", bpo::value<double>()->required(), "Fixed sampling point longitude")
    ("lat_1", bpo::value<double>()->required(), "Fixed point latitude")
    ("n_sample_1", bpo::value<unsigned int>()->required(), "Number of gene copies to sample around center 1")
    ("radius_sample_1", bpo::value<double>()->required(), "Sampling radius around center 1")
    ("n_sample_2", bpo::value<unsigned int>()->required(), "Number of gene copies to sample around center 2")
    ("radius_sample_2", bpo::value<double>()->required(), "Sampling radius around center 2")
    ("sampling_threshold",  bpo::value<unsigned int>()->required(), "Population size under which random sampling is not considered")
    ("suitability_threshold", bpo::value<double>()->required(), "Environmental threshold delimiting suitable (above threshold) and unsuitable (under threshold) areas")
    ("K_max", bpo::value<unsigned int>()->required(), "Carrying capacity in suitable areas")
    ("p", bpo::value<double>()->required(), "Probability of transition between two carrying capacity values in unsuitable areas")
    ("K_min_a", bpo::value<unsigned int>()->required(), "Carrying capacity in unsuitable areas with probability p")
    ("K_min_b", bpo::value<unsigned int>()->required(), "Carrying capacity in unsuitable areas with probability 1-p")
    ("r", bpo::value<double>()->required(), "Constant growth rate")
    ("emigrant_rate", bpo::value<double>()->required(), "Emigrant rate between the four neighboring cells")
    ("friction_suitable", bpo::value<double>()->required(), "Friction coefficient in suitable areas")
    ("friction_unsuitable", bpo::value<double>()->required(), "Friction coefficient in unsuitable areas")
    ("demography_out",  bpo::value<std::string>(), "File name for the simulated demography output")
    ("last_layer_out",  bpo::value<std::string>(), "File name for last demographic layer at sampling time")
    ("distribution_area_out", bpo::value<std::string>(), "File name for coordinates at which population size is positive at sampling time")
    ("mask_out", bpo::value<std::string>(), "File name for coordinates at which random sampling of individuals is operable")
    ("sample_out",  bpo::value<std::string>(), "File name for the simulated sampling scheme output")
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
    std::string m_pop;
    static unsigned int m_next_available_id;
    auto next_id(){
      auto id = m_next_available_id;
      ++m_next_available_id;
      return id;
    }
public:
    Ind(coord_type const& x, std::string pop):
    m_id(next_id()),
    m_x(x),
    m_pop(pop)
    {}
    auto id() const {return m_id;}
    auto x() const {return m_x;}
    auto pop() const {return m_pop;}
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
    landscape_type env( {{"var", file}}, {time_type(0)} );

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

    /******************************
    * Niche and growth functions
    *****************************/
    using quetzal::expressive::literal_factory;
    using quetzal::expressive::use;
    literal_factory<coord_type, time_type> lit;

    // Carrying capacity depends on environment heterogeneity
    auto rain = env["var"];
    double threshold = vm["suitability_threshold"].as<double>();
    auto is_rainfall_high = [rain, threshold](coord_type x, time_type){return (rain(x, 0) >= threshold) ? true : false;};
    double K_max = static_cast<double>(vm["K_max"].as<unsigned int>());
    double K_min_a = static_cast<double>(vm["K_min_a"].as<unsigned int>());
    double K_min_b = static_cast<double>(vm["K_min_b"].as<unsigned int>());
    double p = vm["p"].as<double>();
    auto K = [&gen, is_rainfall_high, K_max, K_min_a, K_min_b, p](coord_type const&x, time_type t){
      if(is_rainfall_high(x,t)){
        return K_max;
      }else{
        return std::bernoulli_distribution(p)(gen) ? K_min_a : K_min_b;
      }
    };

    // Growth rate is constant across the landscape
    auto r = lit(vm["r"].as<double>());

    // Retrieve population size reference to define a logistic growth process
    auto pop_sizes = simulator.pop_size_history();
    auto N = use( [pop_sizes](coord_type x, time_type t){ return pop_sizes(x,t);} );

    auto g = N * ( lit(1) + r ) / ( lit(1) + ( (r * N)/use(K) ));
    auto sim_children = [g](auto& gen, coord_type const&x, time_type t){
      std::poisson_distribution<unsigned int> poisson(g(x,t));
      return poisson(gen);
    };

    /******************
    * Dispersal
    ******************/
    double f_a = vm["friction_suitable"].as<double>();
    double f_b = vm["friction_unsuitable"].as<double>();
    auto friction = [is_rainfall_high, f_a, f_b](coord_type const& x){return (is_rainfall_high(x,0)) ? f_a : f_b ;} ;
    double emigrant_rate = vm["emigrant_rate"].as<double>();
    auto env_ref = std::cref(env);
    auto get_neighbors = [env_ref](coord_type const& x){return env_ref.get().four_nearest_defined_cells(x);};
    auto const& space = env.geographic_definition_space();
    auto dispersal = demographic_policy::make_neighboring_migration(space, emigrant_rate, friction, get_neighbors);

    auto distance = [](coord_type const& a, coord_type const& b) -> double {return memoized_great_circle_distance(a,b);};

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

    if(vm.count("last_layer_out"))
    {
      std::string filename = vm["last_layer_out"].as<std::string>();
      env.export_to_geotiff(N, sample_time -1, sample_time, [&pop_sizes](time_type const& t){return pop_sizes.get().definition_space(t);}, filename);
    }

    /**********************
    *  Sampling schemes
    **********************/

    // Sampling schemes are constrained to be in the last simulated distribution area
    auto distribution_area = pop_sizes.get().definition_space(sample_time);
    // Retrieve the last demographic size distribution
    auto last_N = [N, sample_time](coord_type const& x){return N(x, sample_time);};
    // Filter distribution area to discard less populated demes from sampling space
    auto mask = distribution_area;
    auto sampling_threshold = vm["sampling_threshold"].as<unsigned int>();
    mask.erase(std::remove_if(mask.begin(), mask.end(),
        [last_N, sampling_threshold](coord_type const& x) -> bool { return last_N(x) < sampling_threshold;}
    ));

    if(vm.count("distribution_area_out"))
    {
      std::string filename = vm["distribution_area_out"].as<std::string>();
      env.export_to_shapefile(distribution_area, filename);
    }

    if(vm.count("mask_out"))
    {
      std::string filename = vm["mask_out"].as<std::string>();
      env.export_to_shapefile(mask, filename);
    }

    std::cout << "--- Simulating coalescents" << std::endl;
    unsigned int n_sim_gen = vm["n_sim_gen"].as<unsigned int>();
    boost::progress_display show_progress( n_sim_gen );
    for(unsigned int i_sim_gen = 0; i_sim_gen < n_sim_gen; ++i_sim_gen)
    {
      /************************************************************************
      * Sample 1: Fixed point
      *************************************************************************/
      unsigned int n_1 = vm["n_sample_1"].as<unsigned int>();
      coord_type x_1(vm["lat_1"].as<double>(), vm["lon_1"].as<double>());
      coord_type::km radius_1 = vm["radius_sample_1"].as<double>();
      // Sample points within the sampling radius
      auto density_1 = [radius_1, distance, x_1](coord_type const& x)
      {
        if(distance(x_1, x) < radius_1){return 1.0;} else {return 0.0;};
      };
      auto sampler_1 = schemes::make_constrained_sampler(mask, density_1, last_N, n_1);
      auto sample_1 = sampler_1(gen);

      /************************************************************************
      * Sample 2: Uniformely at random across distribution area
      ************************************************************************/
      unsigned int n_2 = vm["n_sample_2"].as<unsigned int>();
      coord_type::km radius_2 = vm["radius_sample_2"].as<double>();
      std::uniform_int_distribution<int> dist(0, mask.size() - 1);
      auto x_2 = mask.at(dist(gen));
      if(simulator.pop_size_history().get().is_defined(x_2, sample_time) != true)
      {
        throw std::string("invalid sampling seed 2, this error should not be raised, please contact A. Becheler.");
      }
      auto density_2 = [radius_2, distance, x_2](coord_type const& x)
      {
        if(distance(x_2, x) < radius_2){return 1.0;} else {return 0.0;};
      };
      auto sampler_2 = schemes::make_constrained_sampler(mask, density_2, last_N, n_2);
      auto sample_2 = sampler_2(gen);

      /**********************
       *  Coalescence
       **********************/

       std::vector<Ind> v;

       for(auto const& it1 : sample_1){
         for(unsigned int i = 1; i <= it1.second; ++ i){
           v.emplace_back(it1.first, "Pop1");
         }
       }

       for(auto const& it1 : sample_2){
         for(unsigned int i = 1; i <= it1.second; ++ i){
           v.emplace_back(it1.first, "Pop2");
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

      std::string imap;
      for(auto const& it: v){
        imap.append(std::to_string(it.id()) + "\t" + it.pop() + "\n");
      }


      if(vm.count("database")){
        try{
          std::string file = vm["database"].as<std::string>();
          sqlite3pp::database db(file.c_str());
          sqlite3pp::command cmd(db, "CREATE TABLE IF NOT EXISTS results(id INTEGER PRIMARY KEY AUTOINCREMENT UNIQUE, genealogies TEXT, imap TEXT, lat REAL, lon REAL)");
          cmd.execute();

          sqlite3pp::command cmd2(db, "INSERT INTO results (genealogies, imap, lat, lon) VALUES (?,?,?,?)");
          cmd2.binder()<< genealogies
                       << imap
                       << std::to_string(x_2.lat())
                       << std::to_string(x_2.lon());
          cmd2.execute();
          ++show_progress;
        }
        catch(const std::exception& e){
          std::cout << e.what() << std::endl;
        }
      }else{
        std::cout    << "--- Genealogies in Newick format:\n\n"
                     << genealogies << "\n\n"
                     << "--- Mapping:\n\n"
                     << imap << "\n\n"
                     << "--- Coordinates of the center of the variable sampling area:\t"
                     << x_2 << "\n\n" << std::endl;
      }
    } // n_sim_gen
    return simulator;

  }
};
