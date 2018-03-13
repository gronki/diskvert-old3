
find data -type f -delete
python par.py

cd data
parallel --eta cat {} \| dv-mag-rx -compton -post-corona -o {.} ::: {01,02,03,05}-*.par
parallel --eta cat {} \| dv-mag-rx          -post-corona -o {.} ::: 04-*.par
cd ..

parallel python ::: cooling.py maps.py plotpar.py instabil.py
