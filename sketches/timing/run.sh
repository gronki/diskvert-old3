#!/bin/bash

testuj() {

  prefix=$(mktemp -dt build_XXXXXXXXXXXX)

  make -C ../../build -B install FC="$1" FFLAGS="$2" prefix="$prefix"

  echo "prefix: $prefix"

  for k in $(seq 8); do
    cat input.par | LD_LIBRARY_PATH="$prefix/lib:$LD_LIBRARY_PATH" /usr/bin/time -f %U -o "$prefix/time-$k" "$prefix/bin/dv-mag" -top 90 -n 20000 -precision 1e-5 -o "$prefix/disk-$k" >/dev/null 2>&1
  done

  LC_NUMERIC=C echo "$(python stat.py $prefix/time-*) + $1 $2" | tee -a results.txt
  rm -rfv "$prefix"

}

# echo "-------------------------" >> results.txt
# date +%Y-%m-%d\ %H:%M:%S >> results.txt
# # python cpu.py | tee -a results.txt
# echo "" >> results.txt

for math in '' '-funsafe-math-optimizations' '-ffast-math'; do
  for flto in '' '-flto'; do
    for opti in '-Os' '-O2' '-O2 -ftree-vectorize' '-O3'; do
      for march in '' '-march=native' '-mavx2'; do
        for parens in '-fprotect-parens' '-fno-protect-parens'; do
          testuj "gfortran" "$opti $march $math $flto $parens"
        done
      done
    done
  done
done

if which ifort; then
  for math in '-mieee-fp' '-fp-model precise' '-fp-model fast=1'; do
    for flto in '' '-ipo'; do
      for opti in '-Os' '-O2' '-O3'; do
        for march in '' '-xHost' '-mavx2'; do
          for parens in '-fprotect-parens' '-fno-protect-parens'; do
            testuj "ifort" "$opti $march $math $flto $parens"
          done
        done
      done
    done
  done
fi

echo "" >> results.txt
