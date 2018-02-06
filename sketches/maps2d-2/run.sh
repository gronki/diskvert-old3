find . -path './data/*' -delete
python par.py
cd data
cat jobs.lst | parallel --eta bash ../job.sh
cd ..
python plot1.py
