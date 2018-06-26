#!/bin/bash

set -e

samples=48
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

  parallel --eta --jobs '25%' cat $pardir/{}.par \| /usr/bin/time -f %U -o $prefix/time-{} \
    $prefix/bin/diskvert -n $ngrid -compton -post-corona -top 120 -linear -o $prefix/disk-{} \
    ::: $(seq $samples)

  LC_NUMERIC=C echo "$(python stat.py $prefix/time-*) + $(cat $prefix/buildtime) + $FC $FFLAGS" | tee -a results.txt
  rm -rfv "$prefix"

}

rm -f results.txt && touch results.txt

testuj gfortran
testuj gfortran -O0
testuj gfortran -Os
testuj gfortran -Ofast
testuj gfortran -Ofast -flto
testuj gfortran $(rpm -E %optflags)

for lto in '' '-flto'; do
  for math in '' '-ffast-math'; do
    for opt in '-O2' '-O3'; do
      for arch in '' '-march=native'; do
        testuj gfortran $opt $arch $math $lto
      done
    done
  done
done

if which ifort; then
  testuj ifort
  testuj ifort -O0
  testuj ifort -Os
  testuj ifort -fast
  testuj ifort -xhost

  for lto in '' '-ipo'; do
    for opt in '-O2' '-O3'; do
      for arch in '' '-xhost'; do
        testuj ifort $opt $arch $lto
      done
    done
  done
fi

echo '----------------------------------------------------------'
echo '----------------------------------------------------------'
echo

fn=results.$(date +%y%m%d%H%M).txt
echo "$samples samples, $jobs jobs, n = $ngrid" > $fn
cat results.txt | sort -g | tee results.srt.txt | tee -a $fn
