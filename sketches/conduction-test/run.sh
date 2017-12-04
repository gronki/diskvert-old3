dvflags="-no-bf -corona -top 80 -all"

cat input.par | dv-mag-rx $dvflags -cond -o cond
dv-plot-rx cond
diskvert-cooling2D cond.dat

cat input.par | dv-mag-rx $dvflags -o nocond
dv-plot-rx nocond
diskvert-cooling2D nocond.dat

python plot.py
