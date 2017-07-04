rm -f disk.???.png disk.???.dat diskmulti*.png
cat ../input.par | dv-alpha-rx -o disk -n 210
dvpl-alpha-rx-combine $(seq 0 44)
eog -n diskmulti*.png &
dvpl-alpha-rx $(seq 0 44)
eog -n disk.???.png &
