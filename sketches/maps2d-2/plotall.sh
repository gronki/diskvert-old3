
plot() {
  python plot2.py $@ &
  python plot3.py $@ &
  python plot4.py $@ &
  python plot6.py $@ &
  wait
}

for a in $(seq 3)
do
  for b in $(seq 7); do plot $a $b X  X ; done
  for b in $(seq 7); do plot $a $b X  X ; done
  for b in $(seq 5); do plot $a X  $b X ; done
  for b in $(seq 4); do plot $a X  X  $b; done
done
