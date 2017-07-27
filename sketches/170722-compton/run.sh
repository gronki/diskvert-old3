opacities= -no-bf

cat input.par | dv-mag-rx ${opacities} -corona    -o MCFX &
cat input.par | dv-mag-rx ${opacities} -compton   -o MWFX &
cat input.par | dv-mag-rx ${opacities}            -o MDFX &

cat input.par | dv-mag    ${opacities} -compton   -o MWFK &
cat input.par | dv-mag    ${opacities} -compton2  -o MQFK &
cat input.par | dv-mag    ${opacities} -corona    -o MCFK &

wait

python plot.py &

for d in *.dat
do
diskvert-cooling2D --dpi 144 $d
done

eog *.cool2D.png &
