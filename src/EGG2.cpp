//
//  Copyright © 2020 Arnaud Becheler. All rights reserved.
//
#include "EGG2_options.h"
#include "EGG2_context.h"
#include "utils.h"

namespace bpo = boost::program_options;

namespace
{
  const size_t ERROR_IN_COMMAND_LINE = 1;
  const size_t SUCCESS = 0;
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;
}

int main(int argc, char* argv[])
{
  bpo::variables_map vm;
  bool verbose = false;
  try{
    vm = handle_options(argc, argv);
    // --help option
    if (vm.count("help") || vm.count("version") )
    {
      return SUCCESS;
    }
  }
  catch(boost::program_options::required_option& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    return ERROR_IN_COMMAND_LINE;
  }
  catch(boost::program_options::error& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    return ERROR_IN_COMMAND_LINE;
  }
  catch ( const std::exception& e )
  {
    std::cerr << e.what() << std::endl;
    return ERROR_IN_COMMAND_LINE;
}
  // should do the following without fear because everything is required to be present

  // --verbose option
  if (vm.count("verbose"))
  {
    verbose = true;
    EGGS::utils::PrintVariableMap(vm);
  }
  EGGS::utils::PrintVariableMap(vm);

  std::random_device rd;
  std::mt19937 gen(rd());
  if(verbose){std::cout << "Initialization" << std::endl;}
  SimulationContext s(vm, gen, verbose);
  if(verbose){std::cout << "Running ..." << std::endl;}
  try{
    s.run(gen);
  }catch(const std::domain_error &){
    return 1;
  }
  return 0;
}
