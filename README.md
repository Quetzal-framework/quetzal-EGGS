# quetzal-EGGS <img align="right" width="145" src="https://github.com/Becheler/Becheler.github.io/blob/master/draw/logos/quetzal_eggs.png">

[![Becheler](https://circleci.com/gh/Becheler/quetzal-EGGS.svg?style=shield)](https://app.circleci.com/pipelines/github/Becheler)

Suit of coalescence-based simulation programs for spatial population genetics.

Website: https://becheler.github.io/pages/quetzal_eggs/home

## Description

A demographic history is simulated forward in time in a heterogeneous landscape
(the landscape patterns are given to the program options as a geographic raster layer).

The exact demogaphic model is different for each quetzal-EGG, and their parametrization
are different.

At sampling time, the demographic process is stopped and tip nodes
(equivalently: the gene copies sampled across the landscape) coordinates are
read by the program from a CSV file.

A backward in time process begins, that tracks genes ancestors back in the
demographic history. When the Most Recent Common Ancestor of the sample is found,
the simulation stops and records simulation parameters and simulated gene trees
in a SQLite database.

## Installation

### From source with Great Lakes

In the current state of the project, installing the project from source is the only
option. You will need to make sure that the dependencies are met on your system.

If you are using the [Great Lakes Slurm cluster](https://arc.umich.edu/greatlakes/)
from University of Michigan, you can load all the required dependencies by typing:

```
module load cmake/3.17.3 gcc/8.2.0 gdal/3.0.1 boost/1.75.0
```

Then you can build the project typing:

```
git clone --recurse-submodules https://github.com/Becheler/quetzal-EGGS
cd quetzal-EGGS
mkdir Release
cd Release
cmake .. -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc)
cmake --build . --config Release
```

### Dependencies

#### OS

The present project was tested with the following OS:

- Distributor ID:	Ubuntu
- Description:	Ubuntu 20.04.2 LTS
- Release:	20.04
- Codename:	focal

#### Git

The Git utility package is, by default, included in Ubuntu’s software repositories
that can be installed via APT. Enter ```sudo apt install git``` in a terminal
to download and install Git.

Check it has been properly installed with ```git --version```.

#### CMake

CMake is used as an open-source, cross-platform toolbox we use to build, test, and package the software.
Check if it is installed by typing ```cmake --version``` in a terminal.
If it is not installed, please visit the [installation pages](https://vitux.com/how-to-install-cmake-on-ubuntu/).

#### C++ compiler

gcc (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0

####  Boost

Install boost with: ```sudo apt-get install libboost-all-dev```

#### GDAL

The Geospatial Data Abstraction Library (GDAL) is essential to represent a spatially explicit landscape!
To install GDAL please visit: http://www.gdal.org/

### Build

When all dependencies are met, you can simply build the project by running the
following commands in a terminal:

```
git clone --recurse-submodules https://github.com/Becheler/quetzal-EGGS
cd quetzal-EGGS
mkdir Release
cd Release
cmake ..
cmake --build . --config Release
```

# Tests output
If tests are failing and you want to output the reason of the failure, use instead:

```
CTEST_OUTPUT_ON_FAILURE=TRUE cmake --build .
```
