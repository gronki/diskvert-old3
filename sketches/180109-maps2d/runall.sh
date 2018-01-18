for d in data.*
do
    cd "$d"
    cat jobs.lst | parallel -j6 --bar bash ../job.sh
    cd ..
    echo "$d" >> completed.txt
done
