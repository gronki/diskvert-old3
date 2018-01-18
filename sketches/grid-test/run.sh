mkdir -p data

cat input.par | dv-mag-rx -compton -top  60 -linear -o data/lin060
cat input.par | dv-mag-rx -compton -top 120 -linear -o data/lin120
cat input.par | dv-mag-rx -compton -top 240 -linear -o data/lin240
cat input.par | dv-mag-rx -compton -top 480 -linear -o data/lin480

H=240
cat input.par | dv-mag-rx -compton -top $H -linear -o data/lin
cat input.par | dv-mag-rx -compton -top $H -log    -o data/log
cat input.par | dv-mag-rx -compton -top $H -asinh  -o data/ash
cat input.par | dv-mag-rx -compton -top $H -pow2   -o data/gsq

cat input.par | dv-mag-rx -compton -top  60 -log -o data/log060
cat input.par | dv-mag-rx -compton -top 120 -log -o data/log120
cat input.par | dv-mag-rx -compton -top 240 -log -o data/log240
cat input.par | dv-mag-rx -compton -top 480 -log -o data/log480
