
mkdir -p data.{0,1,2}

find data.0 -type f -delete
find data.1 -type f -delete
find data.2 -type f -delete

python par.py

find data.0 -name \*.par | parallel --eta cat {} \| dv-mag-rx -post-corona -o {.}
find data.1 -name \*.par | parallel --eta cat {} \| dv-mag-rx -compton -post-corona -o {.}
find data.2 -name \*.par | parallel --eta cat {} \| dv-mag-rx -corona -o {.}

bash plot.sh
