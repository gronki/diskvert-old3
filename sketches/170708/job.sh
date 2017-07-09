#!/usr/bin/env bash

rysuj() {
diskvert-new-plot ${1}.dat -o ${1}.Z.png
diskvert-new-plot --tau ${1}.dat -o ${1}.T.png
diskvert-cooling2D ${1}.dat -o ${1}.C.png
}

cat D.${1}.par | dv-alpha -no-bf -o ADFK.${1}
tar cfv ADFK.${1}.tar ADFK.$1.{dat,col,txt} && rm -f ADFK.$1.{dat,col,txt}
gzip ADFK.${1}.tar

cat D.${1}.par | dv-mag -no-bf -o MDFK.${1}
rysuj MDFK.${1}
tar cfv MDFK.${1}.tar MDFK.$1.{dat,col,txt} && rm -f MDFK.$1.{dat,col,txt}
gzip MDFK.${1}.tar

cat D.${1}.par | dv-mag -compton2 -no-bf -o MWFK.${1}
rysuj MWFK.${1}
tar cfv MWFK.${1}.tar MWFK.$1.{dat,col,txt} && rm -f MWFK.$1.{dat,col,txt}
gzip MWFK.${1}.tar

# cat D.${1}.par | dv-mag -compton2 -o MWFBK.${1}
# rysuj MWFBK.${1}
# tar cfv MWFBK.${1}.tar MWFBK.$1.{dat,col,txt} && rm -f MWFBK.$1.{dat,col,txt}
# gzip MWFBK.${1}.tar

cat D.${1}.par | dv-mag -corona -no-bf -o MCFK.${1}
rysuj MCFK.${1}
tar cfv MCFK.${1}.tar MCFK.$1.{dat,col,txt} && rm -f MCFK.$1.{dat,col,txt}
gzip MCFK.${1}.tar

echo ${1} | tee -a finished.log
