DE++ - Differential Evolution in C++14
======================================

# Introduction

DE++ is a modern header-only library written in C++14 containing recent
algorithms based on the Differential Evolution (DE) algorithm
([Storn and Price, 1995][TR_DE] and [Storn and Price, 1997][GO_DE]).

It provides an easy to use framework for the development of DE algorithms in
C++, by providing interfaces for DE algorithms and single objective fitness
functions. Furthermore, all Learning-Based Real-Parameter Single Objective
Numerical Optimization problems (CEC-2017), which are defined [here][BF_MO],
are defined and implemented in C++, allowing the user to rapidly test and
compare implementations.

DE++ can be used to:
* Rapidly develop and test new DE algorithms
* Use implemented DE algorithms to optimize a given fitness function

[TR_DE]: http://www1.icsi.berkeley.edu/ftp/pub/techreports/1995/tr-95-012.pdf
[GO_DE]: https://dx.doi.org/10.1023/A:1008202821328
[BF_MO]: http://www.ntu.edu.sg/home/EPNSugan/index_files/CEC2017/CEC2017.htm

# Dependencies
The library depends on:

* [Boost.Random][Boost_Random] (Required) for the generation of random numbers
* [ThreadPool][Thread_Pool] (Optional) for the 4th example
* [Google-test][Googletest_Main] (Optional) only for unit testing

Boost.Random was selected since it is generally faster than <random>.

[Boost_Random]: http://www.boost.org/doc/libs/release/doc/html/boost_random.html
[Thread_Pool]: https://github.com/progschj/ThreadPool
[Googletest_Main]: https://github.com/google/googletest

# Differential Evolution Algorithms

The following algorithms have been implemented so far:

| Algorithm | Citation        | Source code                 |
| --------- | --------------- | --------------------------- |
| DEGL      | [DOI][DEGL]    | include/algorithm/degl.hpp  |
| SHADE     | [DOI][SHADE]   | include/algorithm/shade.hpp |
| L-SHADE   | [DOI][L-SHADE] | include/algorithm/shade.hpp |

[DEGL]: http://dx.doi.org/10.1109/TEVC.2008.2009457
[SHADE]: http://dx.doi.org/10.1109/CEC.2013.6557555
[L-SHADE]: http://dx.doi.org/10.1109/CEC.2014.6900380

# Installation
On a new Linux installation the following must be run:

    sudo apt install build-essential git cmake libboost-random-dev

Then checkout the repository:

    git checkout --recursive https://github.com/tsakirin/deplusplus.git

If you only want to pull the source without its submodules (i.e. ThreadPool
and Google-test), don't use the --recursive flag. Please note that 4th example
and unit tests won't be able to build if you omit the --recursive flag.

If you want to use the CEC-2017 problems, you need to download the codes
folder, available [here][SuganthanCEC].

[SuganthanCEC]: http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC-2017%2fBound-Constrained

Then simply:

    unrar x codes.rar
    mv codes/C\ version/input_data ${deplusplus-root}/cec-2017

where deplusplus-root is the root folder of the project.

# Quick Start
## Building the Examples

Create a build directory and build the examples:

    mkdir build && cd !:1
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j 4 all

The executables are located under the bin folder.

Details for each example are to be found within the respective source code
files.

## Build options

Additional options that can be specified by -DOPTION=ON are the following:

OPTION          | Description
--------------- | -----------
-DBUILD_STATIC  | Builds static binaries
-DBUILD_TESTS   | Builds unit tests using gtest (requires lcov to be installed)
-DBUILD_DOC     | Builds the documentation using [doxygen][Doxygen]

[Doxygen]: http://www.stack.nl/~dimitri/doxygen/

## Tests

The current tests mainly assert that the results of the implementation of all
CEC-2017 functions match the official C version. The use lcov to create coverage
metrics.

In order to build it, first [compile googletest][Googletest_Doc] and then copy
libgtest.a and libgtest_main.a in a new lib folder, in the root directory.

Then install lcov:

    sudo apt install lcov

If you have done all of the above properly, you need to run:

    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON ..
    make -j 4 gtest

Which will automatically build and run all tests.

[Googletest_Doc]: https://github.com/google/googletest/tree/master/googletest

## Documentation

The documentation is written in Doxygen, following the Qt style. To build it,
you need to have doxygen installed.

    sudo apt install doxygen

Following the installation, build the documentation using:

    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DOC=ON ..
    make -C build documentation

