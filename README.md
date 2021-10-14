# quetzal-EGGS <img align="right" width="200" src="https://github.com/Becheler/Becheler.github.io/blob/master/draw/logos/quetzal_eggs.png">

[![Becheler](https://circleci.com/gh/Becheler/quetzal-EGGS.svg?style=shield)](https://app.circleci.com/pipelines/github/Becheler)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Website becheler.github.io](https://img.shields.io/website-up-down-green-red/https/becheler.github.io.svg)](https://becheler.github.io/pages/quetzal_eggs/home)

Quetzal-EGGS is a suit of geospatial programs for simulating neutral molecular diversity of populations.
In their present version, the programs intend to:
1. read parameters of in a georeferenced map and a configuration file
2. simulate a demography history (forward in time)
3. simulate coalescence trees for independent genetic markers (backward in time)
4. record outputs in a SQL database.

It compares to other available simulation resources ( like [SPLATCHE](http://splatche.com/), [simcoal2](http://cmpg.unibe.ch/software/simcoal2/), [egglib](http://mycor.nancy.inra.fr/egglib/index.html), [msprime](http://msprime.readthedocs.io/en/stable/index.html) or [necsim](https://pycoalescence.readthedocs.io/en/release/necsim/necsim_library.html)) by offering original demographic
(forward in time) models.

Being built from [Quetzal-CoalTL](https://github.com/Becheler/quetzal),
this project true strength lies in the ease with which new demographic assumptions
(continental dispersal kernels, oceanic dispersal barriers, niche, density-dependence ...)
can be integrated in a new simulator.

- :crystal_ball: looking forward, we expect to grow the list of existing programs with customized demographic models.
- :email: You are interested? Think of a new model? Want to give some feedback? Don't be shy, [contact me!](https://github.com/Becheler)
- :star: You think this is a cool project? Drop a star on GitHub :point_up:
- :bug: A bug? Oopsie daisy! I'll fix it asap if you [email me](https://github.com/Becheler) or open an issue :point_up:

| Program | Demographic process       |
| --------------| --------------------|
| **EGG1** Stochastic oceanic dispersal (rafting events) | <img src="https://github.com/Becheler/Becheler.github.io/blob/master/movies/animation_EGG1.gif" width="250" height="250"/> |
| **EGG2** Pulses in continental matrix connectivity |  <img src="https://github.com/Becheler/Becheler.github.io/blob/master/movies/animation_EGG2.gif" width="250" height="250"/> |

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

### Using Docker

If you want to try quetzal-EGGS, or just do some tests on local before to go wild on a cluster, Docker is for you.

[Docker](https://www.docker.com/) is a software that can package an application and its dependencies in a virtual
container that can run on any Linux, Windows, or macOS computer. It is more efficient
than virtual machines, but less than system-wide installation.

- [Install Docker](https://docs.docker.com/get-docker/)
- Then, download [our Docker image](https://hub.docker.com/r/arnaudbecheler/quetzal-eggs) by typing ```docker pull arnaudbecheler/quetzal-eggs``` in a terminal
- Enter a docker interactive session with ```docker run --name mycontainer -it arnaudbecheler/quetzal-eggs bash```
- Build and install the project to ```/home/EGGS``` directory with:
```
git clone --recurse-submodules https://github.com/Becheler/quetzal-EGGS
cd quetzal-EGGS
mkdir Release
cd Release
cmake .. -DCMAKE_INSTALL_PREFIX="/home/EGGS"
cmake --build . --config Release --target install
```
You can now leave the build directory to see what is available:
```
cd /home/EGGS/
ls
```
You should see a binary and some default configuration files:
```
EGG1  quetzal_EGG1.config  sample.csv  suitability.tif
```
- ```EGG1``` is the program binary implementing model 1
- ```quetzal_EGG1.config``` is the configuration file associated to the model 1
- ```sample.csv``` lists the coordinates of sampled nodes (tips)
- ```suitability.tif``` is the default landscape
You can test the binary with:
```
./EGG1 --version
```
and access its help menu with:
```
./EGG1 --help   
```
You can run a simulation with the default settings:
```
./EGG1 --config quetzal_EGG1.config --suitability suitability.tif --tips sample.csv
```
**In another terminal**, go to the destination folder of your choice and copy the simulation result with the following command:
```
docker cp mycontainer:/home/EGGS/test_pods.db .
```
You can finally [upload the .db here](https://inloop.github.io/sqlite-viewer/) to visualize the results.
You can exit ```mycontainer``` by typing ```exit```.

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
# Memory
```
git clone --recurse-submodules https://github.com/Becheler/quetzal-EGGS
cd quetzal-EGGS
mkdir RelDeb
cd RelDeb
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
make
valgrind --tool=massif --xtree-memory=full ./src/EGG1 --config ../test/data/quetzal_EGG1.config --tips ../test/data/sample.csv --suitability ../test/data/suitability.tif --duration 50
```
