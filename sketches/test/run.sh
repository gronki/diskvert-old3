dvflags="-linear"

tee input.par <<EOF
mbh    10
mdot   0.05
radius 6
alpha  0.1
eta   0.5
nu     0.0
EOF

parallel <<EOF
cat input.par | dv-mag-rx $dvflags -o diskd
cat input.par | dv-mag-rx $dvflags -o diskds  -quench 
cat input.par | dv-mag-rx $dvflags -o diskw   -compton
cat input.par | dv-mag-rx $dvflags -o diskdp  -post-corona
cat input.par | dv-mag-rx $dvflags -o diskwp  -compton -post-corona
cat input.par | dv-mag-rx $dvflags -o diskwps -compton -post-corona -quench
cat input.par | dv-mag-rx $dvflags -o diskc   -corona 
cat input.par | dv-mag-rx $dvflags -o diskcr  -corona -relcompt
cat input.par | dv-mag-rx $dvflags -o diskcs  -corona -quench 
EOF

parallel python plotone.py ::: disk*.dat
# parallel diskvert-cooling2D ::: disk*.dat

