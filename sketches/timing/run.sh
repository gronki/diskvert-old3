#!/bin/bash

set -e

jobs=3
samples=60
ngrid=1024

pardir=$(mktemp -d)
parallel python mkpar.py \| tee $pardir/{}.par ::: $(seq $samples)
# parallel shuf -i700-1800 -n1 \| tee $pardir/{}.n ::: $(seq $samples)

testuj() {

  prefix=$(mktemp -d)

  FC=$1
  shift 1
  FFLAGS="$@"

  make -C ../../build clean
  /usr/bin/time -f %U -o $prefix/buildtime make -C ../../build install \
    CC="gcc" FC="$FC" CFLAGS="-O2" FFLAGS="$FFLAGS" prefix="$prefix"

  echo "prefix: $prefix"

  parallel --eta --jobs $jobs cat $pardir/{}.par \| /usr/bin/time -f %U -o $prefix/time-{} \
    $prefix/bin/diskvert -n $ngrid -compton -post-corona -top 120 -linear -o $prefix/disk-{} \
    ::: $(seq $samples)

  LC_NUMERIC=C echo "$(python stat.py $prefix/time-*) + $(cat $prefix/buildtime) + $FC $FFLAGS" | tee -a results.txt
  rm -rfv "$prefix"

}

rm -f results.txt && touch results.txt

testuj gfortran -O0
testuj gfortran -Os
testuj gfortran -Ofast
testuj gfortran -O2
testuj gfortran -O3 -march=native
testuj gfortran -O3 -march=native -flto

# for math in '' '-funsafe-math-optimizations' '-funsafe-math-optimizations -fno-protect-parens'; do
  # for lto in '' '-flto'; do
  #   for march in '' '-march=native'; do
  #     for opti in '-O2' '-O3'; do
  #     testuj gfortran $opti $march $math $lto
  #     done
  #   done
  # done
# done

if which ifort; then
  testuj ifort -O0
  testuj ifort -Os
  testuj ifort -fast
  for lto in '' '-ipo'; do
    for march in '' '-xhost'; do
      for opti in '-O2' '-O3'; do
        testuj ifort $opti $march $lto
      done
    done
  done
fi

echo '----------------------------------------------------------'
echo '----------------------------------------------------------'
echo

fn=results.$(date +%y%m%d%H%M).txt
echo "$samples samples, $jobs jobs, n = $ngrid" > $fn
cat results.txt | sort -g | tee -a $fn
