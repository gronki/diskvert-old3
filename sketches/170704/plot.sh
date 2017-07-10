

#rm MCFX.* && cat input.par | dv-mag-rx -top 37 -n 500 -all -corona -no-bf -o MCFX

frames=$@

rm -f MCFX.*.png

pla() {
    python plot_temp.py MCFX $frames
    eog MCFX.T???.png &
}

plb() {
    for i in $frames
    do
        diskvert-cooling2D --force MCFX.dat[$i]  &
	if test $(( (i+1) % $(nproc) )) == 0; then wait; fi
    done

    wait && eog MCFX.???.cool2D.png &
}


pla &
plb &

wait
