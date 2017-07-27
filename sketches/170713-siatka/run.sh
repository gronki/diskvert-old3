#!/usr/bin/env bash

rm -f *.{png,col,dat,txt,log,tgz}


cat input.par | dv-alpha-rx  -no-bf -o ADFX &
cat input.par | dv-alpha-rx   -o ADFBX &
cat input.par | dv-mag-rx  -no-bf -o MDFX &
cat input.par | dv-mag-rx   -o MDFBX &
wait
cat input.par | dv-mag-rx -corona -no-bf -o MCFX &
cat input.par | dv-mag-rx -corona  -o MCFBX &
wait
cat input.par | dv-alpha  -no-bf -o ADFK &
cat input.par | dv-alpha   -o ADFBK &
cat input.par | dv-mag  -no-bf -o MDFK &
cat input.par | dv-mag   -o MDFBK &
wait
cat input.par | dv-mag -corona -no-bf -o MCFK &
cat input.par | dv-mag -corona  -o MCFBK &
wait
python plot.py
tar czf MDFK.tgz MDFK.{dat,col,txt}
diskvert-cooling2D MDFK.tgz
tar czf MDFBK.tgz MDFBK.{dat,col,txt}
diskvert-cooling2D MDFBK.tgz
tar czf MCFK.tgz MCFK.{dat,col,txt}
diskvert-cooling2D MCFK.tgz
tar czf MCFBK.tgz MCFBK.{dat,col,txt}
diskvert-cooling2D MCFBK.tgz
