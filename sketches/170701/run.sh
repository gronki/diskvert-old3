rm -f disk.???.png disk.???.dat diskmulti*.png
cat input.par | dv-mag-rx -o disk -n 500 -top 24
#dvpl-alpha-rx-combine $(seq 0 100)
#eog -n diskmulti*.png &
dvpl-alpha-rx $(seq 0 100)
eog -n disk.???.png &
