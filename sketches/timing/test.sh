#!/bin/bash

testuj() {
prefix=$(mktemp -d)
make -C ../.. distclean
make -C ../.. install FC="$1" FFLAGS="$2" prefix="$prefix"
cat input.par | /usr/bin/time -f %U\ %P -o A "$prefix/bin/dv-mag" -top 120 -n 70000 -precision 1e-7 
echo "+ $(cat A) $1 $2" | tee -a results.txt
}

echo "-------------------------" >> results.txt
date +%Y-%m-%d\ %H:%M:%S >> results.txt
python cpu.py | tee -a results.txt
echo "" >> results.txt

testuj "gfortran" "-O0"
testuj "gfortran" "-O1"
testuj "gfortran" "-O2"
testuj "gfortran" "-O2 -march=native"
testuj "gfortran" "-O2 -march=native -mtune=generic"
testuj "gfortran" "-O2 -mavx2"
testuj "gfortran" "-O3 -mavx2"
testuj "gfortran" "-O3"
testuj "gfortran" "-O3 -march=native"
testuj "gfortran" "-O3 -march=native -ffast-math"
testuj "gfortran" "-O3 -march=native -funsafe-math-optimizations"
testuj "gfortran" "-Ofast"
testuj "gfortran" "-Ofast -mavx2"
testuj "gfortran" "-Ofast -march=native"
testuj "gfortran" "-Os"

if which ifort; then
  testuj "ifort" "-O0"
  testuj "ifort" "-O2"
  testuj "ifort" "-O2 -xHost"
  testuj "ifort" "-O2 -mieee-fp"
  testuj "ifort" "-O2 -mieee-fp -xHost"
  testuj "ifort" "-O3"
  testuj "ifort" "-fast"
  testuj "ifort" "-O3 -xHost"
  testuj "ifort" "-O3 -mieee-fp"
  testuj "ifort" "-O3 -mieee-fp -xHost"
fi

echo "" >> results.txt
