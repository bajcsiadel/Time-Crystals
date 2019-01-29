#!/bin/bash

mkdir -p ./results/outs ./results/movies ./results/stats
cd ./source
make
cd ../

for param in `ls -a ./parameters | egrep '^[0-9]+ero.json$'`
do
    NAME=`echo $param | cut -d'.' -f1`
    DATE=`date +%Y%m%d_%H%M%S`
    build/softsim parameters/$param > results/outs/out_$NAME\_$DATE.txt &

    nr_runnings=`ps -e | grep "softsim" | wc -l`
    while [ $nr_runnings -ge 30 ]
    do
        sleep 3
        nr_runnings=`ps -e | grep "softsim" | wc -l`
    done
done
