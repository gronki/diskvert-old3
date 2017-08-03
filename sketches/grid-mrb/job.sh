DVFLAGS="-no-bf -top 150"

cat $1.par | dv-mag-rx $DVFLAGS          -o "${1}D"
cat $1.par | dv-mag-rx $DVFLAGS -compton -o "${1}W"
cat $1.par | dv-mag-rx $DVFLAGS -corona  -o "${1}C"

