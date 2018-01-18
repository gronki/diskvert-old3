DVFLAGS="-no-bf"

cat "${1}.par" | dv-mag-rx $DVFLAGS -compton -o "${1}W"
cat "${1}.par" | dv-mag-rx $DVFLAGS -corona  -o "${1}C"
