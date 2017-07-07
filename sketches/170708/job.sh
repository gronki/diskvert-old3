#!/usr/bin/env bash

rysuj() {
diskvert-new-plot --dpi 144 ${1}.dat -o ${1}.Z.png
diskvert-new-plot --dpi 144 --tau ${1}.dat -o ${1}.T.png
diskvert-cooling2D --dpi 144 ${1}.dat -o ${1}.C.png
}

cat D.${1}.par | dv-mag -compton -no-bf -o MWFK.${1}
rysuj MWFK.${1}

cat D.${1}.par | dv-mag -compton -o MWFBK.${1}
rysuj MWFBK.${1}

cat D.${1}.par | dv-mag -corona -no-bf -o MCFK.${1}
rysuj MCFK.${1}


echo ${1} | tee -a finished.log
