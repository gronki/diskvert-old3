
## How to?

First, create the data directory

```bash
mkdir data
```

Then you can do the computation:

```bash
# generates files
python par.py
# do the computation
cd data
parallel --eta cat {} \| dv-mag-rx -compton -post-corona -o {.} ::: *.par
cd ..
# make plots
python maps.py
python modelpar.py
python cooling.py
python instabil.py
```

To clear the data directory, invoke:
```sh
find data -type f -delete
```
