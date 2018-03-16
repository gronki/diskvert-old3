#!/bin/bash

testuj() {

  prefix=$(mktemp -dt build_XXXXXXXXXXXX)

  make -C ../../build clean
  make -C ../../build -j6 install CC="gcc" FC="$1" FFLAGS="$2" prefix="$prefix"

  echo "prefix: $prefix"

  for k in $(seq 7); do
    cat input.par | /usr/bin/time -f %U -o "$prefix/time-$k" "$prefix/bin/dv-mag-rx" -n 600 -corona -post-corona -top 120 -linear -precision 1e-5 -o "$prefix/disk-$k"
    # cat input.par | $prefix/bin/dv-mag-rx -n 900 -corona -post-corona -top 120 -linear -precision 1e-5 -o "$prefix/disk-$k" | grep PERF | sed s/PERF// | tee $prefix/time-$k
  done

  LC_NUMERIC=C echo "$(python stat.py $prefix/time-*) + $1 $2" | tee -a results.txt
  rm -rfv "$prefix"

}

printf "" > results.txt

for math in '-mieee-fp'; do
  for flto in '' '-flto'; do
    for opti in '-O2' '-O2 -ftree-vectorize' '-O2 -finline-functions' '-O3'; do
      for march in '' '-march=native'; do
        testuj "gfortran" "$math $opti $flto $march"
      done
    done
  done
done

if which ifort; then
  for math in '-fp-model fast=1' '-fp-model precise'; do
    for flto in '' '-ipo'; do
      for opti in '-O2' '-O3'; do
        for march in '' '-xhost'; do
          testuj "ifort" "$math $opti $flto $march"
        done
      done
    done
  done
fi

cat results.txt | sort -g | tee results.srt.txt

