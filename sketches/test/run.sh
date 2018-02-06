dvflags=""

cat input.par | dv-mag-rx $dvflags -o diskd
cat input.par | dv-mag-rx $dvflags -o diskw -compton
cat input.par | dv-mag-rx $dvflags -o diskp -compton -post-corona
cat input.par | dv-mag-rx $dvflags -o diskq -post-corona
cat input.par | dv-mag-rx $dvflags -o diskc -corona 

parallel python plotone.py ::: disk*.dat
parallel diskvert-cooling2D ::: disk*.dat

