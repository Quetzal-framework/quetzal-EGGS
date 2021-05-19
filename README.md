# quetzal-EGGS
suit of coalescence-based simulation programs

# Build
```
git clone --recurse-submodules https://github.com/Becheler/quetzal-EGGS
cd quetzal-EGGS
mkdir Release
cd Release
cmake ..
cmake --build . --config Release
```

# Tests output
```
CTEST_OUTPUT_ON_FAILURE=TRUE cmake --build .
```
