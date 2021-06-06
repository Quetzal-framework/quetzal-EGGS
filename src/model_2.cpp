//
//  Copyright © 2018 Arnaud Becheler. All rights reserved.
//
//
// g++ main.cpp -o main -I/usr/include/gdal -lgdal -I/home/becheler/dev -std=c++17 -I/home/becheler/dev -lboost_program_options -lsqlite3
//
//
//
#include "model_2.h"

int main(int argc, char* argv[])
{

  std::cout << "Reading options" << std::endl;
  auto vm = handle_options(argc, argv);
  if (vm.count("help")) {
      return 1;
  }
  try{
    std::cout << "Simulation" << std::endl;
    auto result = SimulationContext::run(vm);
  }
  catch(const std::exception& e){
    std::cout << e.what() << std::endl;
  }
  return 0;
}
