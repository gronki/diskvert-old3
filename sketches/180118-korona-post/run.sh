dvf=""

cat input.par | dv-mag-rx $dvf          -o diskd
cat input.par | dv-mag-rx $dvf -compton -o diskw
cat input.par | dv-mag-rx $dvf          -post-corona -o diskd2
cat input.par | dv-mag-rx $dvf -compton -post-corona -o diskw2
cat input.par | dv-mag-rx $dvf -corona  -o diskc
