rm -f *.{col,txt,dat,log,png}

cat input.par | dv-mag-rx -no-bf -corona -all -o MCFX &
cat input.par | dv-mag-rx -no-bf -o MDFX &
cat input.par | dv-mag -no-bf -o MDFK &
cat input.par | dv-mag -no-bf -compton -o MWFK &
cat input.par | dv-mag -no-bf -corona -o MCFK &

wait

diskvert-cooling2D MDFK.dat &
diskvert-cooling2D MDFX.dat &
diskvert-cooling2D MWFK.dat &
diskvert-cooling2D MCFK.dat &

wait

itmax=$(cat MCFX.txt | grep niter | awk '{ print $2 }')
bash plot.sh $(seq 0 $itmax)
