//
//  Copyright © 2020 Arnaud Becheler. All rights reserved.
//
#include "qegg1.h"

namespace
{
  const size_t ERROR_IN_COMMAND_LINE = 1;
  const size_t SUCCESS = 0;
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;

}

int main(int argc, char* argv[])
{
  bool verbose = false;
  auto vm = handle_options(argc, argv);

  // --help option
  if (vm.count("help"))
  {
    std::cout << "This is Quetzal-EGG 1 simulator" << std::endl;
    std::cout << "Generates Newick gene trees in heterogeneous landscape." << std::endl;
    std::cout << "Author: Arnaud Becheler, 2020." << std::endl;
    PrintVariableMap(vm);
    return SUCCESS;
  }

  // --verbose option
  if (vm.count("v"))
  {
    verbose = true;
    PrintVariableMap(vm);
  }

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
