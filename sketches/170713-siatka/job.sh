J=${1%.par}

DVFLAGS="-no-bf"

cat $J.par | dv-mag-rx $DVFLAGS          -o "${J}D"
cat $J.par | dv-mag-rx $DVFLAGS -compton -o "${J}W"
cat $J.par | dv-mag-rx $DVFLAGS -corona  -o "${J}C"

python plot.py $J

