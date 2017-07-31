#!/usr/bin/env bash
# coding: utf-8

cat input.par | dv-mag-rx -no-bf -corona -all -o MCFX &
cat input.par | dv-mag-rx -no-bf -o MDFX &
cat input.par | dv-mag -no-bf -corona -o MCFK &

wait

diskvert-cooling2D --force MDFX.dat &
diskvert-cooling2D MCFK.dat &
diskvert-cooling2D --force MCFX.dat &

wait

eog M{DFX,CFK,CFX}.cool2D.png &

itmax=$(cat MCFX.txt | grep niter | awk '{ print $2 }')

parallel python plot.py MCFX ::: $(seq 0 $itmax)
parallel diskvert-cooling2D --force MCFX.dat[{}] ::: $(seq 0 $itmax)

