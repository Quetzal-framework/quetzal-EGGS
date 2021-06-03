# quetzal-EGGS <img align="right" width="200" src="https://github.com/Becheler/Becheler.github.io/blob/master/draw/logos/quetzal_eggs.png">

[![Becheler](https://circleci.com/gh/Becheler/quetzal-EGGS.svg?style=shield)](https://app.circleci.com/pipelines/github/Becheler)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Quetzal-EGGS is a suit of spatial coalescence simulation programs for neutral populations.
In their present version, the programs intend to read spatial information about
model parameters in a map, simulate coalescence trees for independent genetic markers and
write outputs in a SQL database.

It compares to other available simulation resources ( like [SPLATCHE](http://splatche.com/), [simcoal2](http://cmpg.unibe.ch/software/simcoal2/), [egglib](http://mycor.nancy.inra.fr/egglib/index.html), [msprime](http://msprime.readthedocs.io/en/stable/index.html) or [necsim](https://pycoalescence.readthedocs.io/en/release/necsim/necsim_library.html)) by offering original demographic
(forward in time) models. Being built from [Quetzal-CoalTL](https://github.com/Becheler/quetzal), the
true strength of this project is how easily demographic assumptions (continental dispersal kernels, oceanic dispersal barriers,
  niche, density-dependence ...)  can be integrated in a new simulator.

- In the future, there will be (hopefully) several programs with different demographic models.
- But right now, it's work in progress, so there is only one model available!

<!-- Place this tag where you want the button to render. -->
<a class="github-button" href="https://github.com/ntkme/github-buttons" data-icon="octicon-star" aria-label="Star ntkme/github-buttons on GitHub">Star</a> Star us on GitHub — it keeps me motivated!

You can read more on the (to-be-updated) website: https://becheler.github.io/pages/quetzal_eggs/home

## Description <img align="right" width="400" src="https://github.com/Becheler/Becheler.github.io/blob/master/draw/quetzal-EGGS/readme_scheme.svg">

A demographic history is simulated forward in time in a heterogeneous landscape
(the landscape patterns are given to the program options as a geographic raster layer).

The exact demogaphic model is different for each quetzal-EGG, and their parametrization
are different.

At sampling time, the demographic process is stopped and tip nodes coordinates
(that is, the gene copies that have been sampled across the landscape) are
read by the program from a CSV file.

A backward in time process begins, that tracks genes ancestors back in the
demographic history. When the Most Recent Common Ancestor of the sample is found,
the simulation stops and records simulation parameters and simulated gene trees
in a SQLite database.

## Installation

In the current state of the project, installing the project from source is the only
option available. You will need to make sure that the dependencies are met on your system.

### From source with the Great Lakes Slurm Cluster

If you are a user of the [Great Lakes Slurm cluster](https://arc.umich.edu/greatlakes/)
from University of Michigan, you can load all the required dependencies by typing:

```
module load cmake/3.17.3 gcc/8.2.0 gdal/3.0.1 boost/1.75.0
```

Then build the project with:

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
