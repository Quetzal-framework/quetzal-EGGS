//
//  Copyright © 2018 Arnaud Becheler. All rights reserved.
//
//
// g++ main.cpp -o main -I/usr/include/gdal -lgdal -I/home/becheler/dev -std=c++17 -I/home/becheler/dev -lboost_program_options -lsqlite3
//
//
//
#include "qegg1.h"

int main(int argc, char* argv[])
{
  bool verbose = false;
  auto vm = handle_options(argc, argv);
  if (vm.count("help")) {
      return 1;
  }
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
