N=1024
DVFLAGS= -no-bf

cat $1.par | dv-mag-rx -n $N $DVFLAGS          -o ${1}D
cat $1.par | dv-mag-rx -n $N $DVFLAGS -compton -o ${1}W
cat $1.par | dv-mag-rx -n $N $DVFLAGS -corona  -o ${1}C
