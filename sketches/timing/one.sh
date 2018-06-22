#!/bin/bash

set -e

jobs=12
samples=48
ngrid=1024

prefix=$(mktemp -d)
echo "prefix: $prefix"

parallel --eta --jobs $jobs python mkpar.py \| /usr/bin/time -f %U -o $prefix/time-{} \
  diskvert -n $ngrid -compton -post-corona -top 120 -linear -o $prefix/disk-{} \
  ::: $(seq $samples)

python stat.py $prefix/time-*
rm -rf "$prefix"
