DVFLAGS="-no-bf -top 140"

mkdir -p data
cat "${1}.par" | dv-mag-rx $DVFLAGS          -o "data/${1}D"
cat "${1}.par" | dv-mag-rx $DVFLAGS -compton -o "data/${1}W"
cat "${1}.par" | dv-mag-rx $DVFLAGS -corona  -o "data/${1}C"

mv "${1}.par" data/

mkdir -p img cool
diskvert-cooling2D "data/${1}C.dat" -o "cool/${1}C.png"
