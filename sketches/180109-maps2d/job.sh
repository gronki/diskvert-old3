DVFLAGS="-no-bf -top 120 -n 840"
PARFILE="par/${1}.par"

mkdir -p datac dataw
cat "$PARFILE" | dv-mag-rx $DVFLAGS -compton -o "dataw/${1}"
cat "$PARFILE" | dv-mag-rx $DVFLAGS -corona  -o "datac/${1}"
