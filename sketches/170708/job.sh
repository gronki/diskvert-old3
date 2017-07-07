#!/usr/bin/env bash

cat D-${1}.par | dv-mag -compton -no-bf -o MWFK-${1}
diskvert-new-plot MWFK-${1}.dat
diskvert-cooling2D MWFK-${1}.dat

cat D-${1}.par | dv-mag -corona -no-bf -o MCFK-${1}
diskvert-new-plot MCFK-${1}.dat
diskvert-cooling2D MCFK-${1}.dat

cat D-${1}.par | dv-mag -compton -o MWFBK-${1}
diskvert-new-plot MWFBK-${1}.dat
diskvert-cooling2D MWFBK-${1}.dat

cat D-${1}.par | dv-mag -corona -o MCFBK-${1}
diskvert-new-plot MCFBK-${1}.dat
diskvert-cooling2D MCFBK-${1}.dat

echo ${1}
