JOBID="${1%.par}"

cat "$1" | dv-mag-rx -compton $DVFLAGS -o "${JOBID}W"
cat "$1" | dv-mag-rx -corona  $DVFLAGS -o "${JOBID}C"
