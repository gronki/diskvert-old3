#!/usr/bin/env bash

#F=MCFK
#cat $1.par | tee $F.$1.par | dv-mag -corona -no-bf -o $F.$1
#tar czf $F.$1.tar.gz $F.$1.{dat,col,txt,log,par}
#rm -f $F.$1.{dat,col,txt,log,par}

F=MWFK
cat $1.par | tee $F.$1.par | dv-mag -compton2 -n 12000 -no-bf -o $F.$1
tar czf $F.$1.tar.gz $F.$1.{dat,col,txt,log,par}
rm -f $F.$1.{dat,col,txt,log,par}

F=MWFBK
cat $1.par | tee $F.$1.par | dv-mag -compton2 -n 12000 -o $F.$1
tar czf $F.$1.tar.gz $F.$1.{dat,col,txt,log,par}
rm -f $F.$1.{dat,col,txt,log,par}

F=MCFX
cat $1.par | tee $F.$1.par | dv-mag-rx -n 1200 -no-bf -o $F.$1
tar czf $F.$1.tar.gz $F.$1.{dat,col,txt,par}
rm -f $F.$1.{dat,col,txt,par}

F=MCFBX
cat $1.par | tee $F.$1.par | dv-mag-rx -n 1200 -o $F.$1
tar czf $F.$1.tar.gz $F.$1.{dat,col,txt,par}
rm -f $F.$1.{dat,col,txt,par}

rm -f $1.par

echo $1 | tee -a finished.log
