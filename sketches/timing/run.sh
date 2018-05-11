#!/bin/bash

set -e

trials=10
jobs=2

testuj() {

  prefix=$(mktemp -dt build_XXXXXXXXXXXX)
  FC=$1
  shift 1
  FFLAGS="$@"

  make -C ../../build clean
  make -C ../../build install CC="gcc" FC="$FC" CFLAGS="-O2" FFLAGS="$FFLAGS" prefix="$prefix"

  echo "prefix: $prefix"

  parallel --eta --jobs $jobs --delay 0.2s cat input.par \| /usr/bin/time -f %U -o $prefix/time-{} $prefix/bin/dv-mag-rx -n 1024 -corona -post-corona -top 120 -linear -o $prefix/disk-{} ::: $(seq $trials)

  LC_NUMERIC=C echo "$(python stat.py $prefix/time-*) + $FC $FFLAGS" | tee -a results.txt
  rm -rfv "$prefix"

}

printf "" > results.txt

testuj gfortran '-O0'
testuj gfortran '-O1'
testuj gfortran '-Og'
testuj gfortran '-Os'

for math in '' '-funsafe-math-optimizations'; do
  for lto in '' '-flto'; do
    for march in '' '-march=native'; do
      for opti in '-O2' '-O2 -ftree-vectorize' '-O3'; do
      testuj gfortran $opti $march $math $lto
      done
    done
  done
done

if which ifort; then
  for math in '-fp-model fast=1' '-fp-model precise'; do
    for lto in '' '-ipo'; do
      for march in '' '-xhost'; do
        for opti in '-O2' '-O3'; do
          testuj ifort $opti $march $math $lto
        done
      done
    done
  done
fi

echo '----------------------------------------------------------'
echo '----------------------------------------------------------'
echo
cat results.txt | sort -g
