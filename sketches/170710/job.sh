#!/usr/bin/env bash

set -e

F=MCFK
cat D.$1.par | dv-mag -corona -no-bf -o $F.$1
tar czf  $F.$1.tar.gz $F.$1.{dat,col,txt,log} D.$1.par
rm -f $F.$1.{dat,col,txt,log}

F=MZFK
cat D.$1.par | dv-mag -compton2 -n 12000 -no-bf -o $F.$1
tar czf  $F.$1.tar.gz $F.$1.{dat,col,txt,log} D.$1.par
rm -f $F.$1.{dat,col,txt,log}

rm -f D.$1.par

echo $1 | tee -a finished.log
