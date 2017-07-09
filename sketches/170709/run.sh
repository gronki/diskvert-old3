#!/usr/bin/env bash

cat input.par | dv-mag -corona -no-bf -o MCFK &
cat input.par | dv-mag -compton -no-bf -o MWFK &
cat input.par | dv-mag -compton2 -no-bf -o MZFK &

wait

python plot.py &

for m in MCFK MWFK MZFK
do
	diskvert-cooling2D $m.dat -o $m.C.png &
	diskvert-new-plot $m.dat -o $m.Z.png &
	diskvert-new-plot --tau $m.dat -o $m.T.png &
	wait
done


