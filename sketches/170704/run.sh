
f() {
	i=$1
	dv-plot-rx MCFX $i
	diskvert-cooling2D --force MCFX.dat[$i]
}
set -e

rm MCFX.*
cat input.par | dv-mag-rx -top 37 -n 600 -all -corona -no-bf -o MCFX

for i in $(seq 44 120)
do
f $i &
if test $(( (i+1) % 12 )) == 0; then
	wait
fi
done

wait

animate -delay 9 MCFX.???.png &
animate -delay 9 MCFX.???.cool2D.png &

