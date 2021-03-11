#!/usr/bin/bash
c++ -std=c++17 -shared -g0 -O3 -march=native -Wall -ffast-math -fPIC $(python3 -m pybind11 --includes) healpywrapper_test.cc -o healpywrapper_test$(python3-config --extension-suffix) -I .


