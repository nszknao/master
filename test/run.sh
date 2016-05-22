#!/bin/bash

export LD_LIBRARY_PATH=/usr/local/lib

if [ -e gatest.exe ]; then
        rm gatest.exe
fi
############ Compile-force #############################################
#g++ -g -O0 -lstdc++ ../test/testMain.cpp ../ga/src/Ga.cpp ../ga/src/GaCommon.cpp -o gatest.exe -lm -lgsl -lgslcblas -std=c++11
g++49 -g -O0 -lstdc++ -I /usr/include/boost148 ../ga/src/SchwefelEllipsoidCMA.cpp -o gatest.exe -lm -lgsl -lgslcblas -std=c++11

if [ -e gatest.exe ]; then
        ./gatest.exe
fi
