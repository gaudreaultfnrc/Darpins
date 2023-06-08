#!/bin/sh

DEBUG=0
FLAGS=
if [ $DEBUG = 1 ]; then
	FLAGS="-g -Wall"
else
	FLAGS="-O2"
fi
CFLAGS=$FLAGS
CXXFLAGS=$FLAGS
LDFLAGS=$FLAGS

PYTHON_INCLUDES=
PYTHON_LIBRARIES=

echo "compiling levenshtein"
# levenshtein
swig -c++ -python -shadow levenshtein.i
g++ -fPIC -c levenshtein_wrap.cxx -I${PYTHON_INCLUDES}
g++ -shared levenshtein_wrap.o -o _levenshtein.so -Xlinker -rpath -L${PYTHON_LIBRARIES} -lpython3.5m
cp _levenshtein.so levenshtein.py ..
