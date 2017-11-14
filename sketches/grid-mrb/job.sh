DVFLAGS="-no-bf -top 150"

cat "par/${1}" | dv-mag-rx $DVFLAGS          -o "data/${1}D"
cat "par/${1}" | dv-mag-rx $DVFLAGS -compton -o "data/${1}W"
cat "par/${1}" | dv-mag-rx $DVFLAGS -corona  -o "data/${1}C"

python plot.py "${1}"
diskvert-cooling2D "data/${1}C.dat" -o "cool/${1}C.png"
