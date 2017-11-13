#!/usr/bin/env bash

cat input.par | dv-mag-rx -corona -no-opacity-bf -o relx
diskvert-cooling2D relx.dat
display relx.cool2D.png &

cat input.par | dv-mag -o rkm -balance-multi -no-bf
python plot.py
display temp.png &
