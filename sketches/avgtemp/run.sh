dvflags='-no-bf'

cat input.par | dv-mag-rx ${dvflags} -corona    -o MCFX
cat input.par | dv-mag-rx ${dvflags} -compton   -o MWFX
cat input.par | dv-mag-rx ${dvflags}            -o MDFX

