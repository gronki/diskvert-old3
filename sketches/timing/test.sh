#!/bin/bash

testuj() {

  prefix=$(mktemp -d)

  make -C ../.. distclean
  make -C ../.. install FC="$1" FFLAGS="$2" prefix="$prefix"

  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$prefix/lib"

  cat input.par | /usr/bin/time -f %U -o T1 "$prefix/bin/dv-mag" -top 80 -n 10000 -precision 1e-7
  cat input.par | /usr/bin/time -f %U -o T2 "$prefix/bin/dv-mag" -top 80 -n 30000 -precision 1e-7
  cat input.par | /usr/bin/time -f %U -o T3 "$prefix/bin/dv-mag" -top 80 -n 90000 -precision 1e-7

  LC_NUMERIC=C printf "+ %4.1f %4.1f %5.1f + %s %s\n" $(cat T1) $(cat T2) $(cat T3) "$1" "$2" | tee -a results.txt
  rm -rfv T1 T2 T3 "$prefix"

}

echo "-------------------------" >> results.txt
date +%Y-%m-%d\ %H:%M:%S >> results.txt
python cpu.py | tee -a results.txt
echo "" >> results.txt

testuj "gfortran" "-O0"
testuj "gfortran" "-Os"

testuj "gfortran" "-O2"
# testuj "gfortran" "-O2 -flto"
testuj "gfortran" "-O2 -mavx2"
testuj "gfortran" "-O2 -march=native"
# testuj "gfortran" "-O2 -march=native -flto"
testuj "gfortran" "-O2 -march=native -funsafe-math-optimizations"
# testuj "gfortran" "-O2 -march=native -funsafe-math-optimizations -flto"

testuj "gfortran" "-O3"
# testuj "gfortran" "-O3 -flto"
testuj "gfortran" "-O3 -mavx2"
testuj "gfortran" "-O3 -march=native"
# testuj "gfortran" "-O3 -march=native -flto"
testuj "gfortran" "-O3 -march=native -funsafe-math-optimizations"
# testuj "gfortran" "-O3 -march=native -funsafe-math-optimizations -flto"
testuj "gfortran" "-O3 -march=native -ffast-math"

if which ifort; then
  testuj "ifort" "-O0"
  testuj "ifort" "-O2"
  testuj "ifort" "-O2 -fp-model precise"
  testuj "ifort" "-O2 -xHost"
  testuj "ifort" "-O2 -xHost -fp-model precise"
  testuj "ifort" "-O2 -xHost -fp-model precise -ipo"
  testuj "ifort" "-O2 -mieee-fp"
  testuj "ifort" "-O2 -mieee-fp -xHost"
  testuj "ifort" "-O3"
  testuj "ifort" "-fast"
  testuj "ifort" "-O3 -xHost"
  testuj "ifort" "-O3 -xHost -fp-model precise"
  testuj "ifort" "-O3 -xHost -fp-model precise -ipo"
  testuj "ifort" "-O3 -xHost -ipo"
  testuj "ifort" "-O3 -mieee-fp"
  testuj "ifort" "-O3 -mieee-fp -xHost"
fi

echo "" >> results.txt
