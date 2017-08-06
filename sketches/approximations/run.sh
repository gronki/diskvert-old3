dvflags='-no-bf'

#cat - | parallel <<EOF
cat input.par | dv-mag-rx ${dvflags} -corona    -o MCFX
cat input.par | dv-mag-rx ${dvflags} -compton   -o MWFX
cat input.par | dv-mag-rx ${dvflags}            -o MDFX
cat input.par | dv-mag    ${dvflags} -compton   -o MWFK
cat input.par | dv-mag    ${dvflags} -compton2  -o MQFK
cat input.par | dv-mag    ${dvflags} -corona    -o MCFK
#EOF

python plot.py &

for d in *.dat
do
diskvert-cooling2D --dpi 144 $d
done

eog *.cool2D.png &
