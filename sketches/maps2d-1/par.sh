for a in $(seq 3)
do
  for b in $(seq 7); do python par.py $a $b X  X ; done
  for b in $(seq 7); do python par.py $a $b X  X ; done
  for b in $(seq 5); do python par.py $a X  $b X ; done
  for b in $(seq 4); do python par.py $a X  X  $b; done
done
