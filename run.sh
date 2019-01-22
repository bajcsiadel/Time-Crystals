#!/bin/bash

mkdir -p ./results/outs ./results/movies ./results/stats
cd ./source
make
cd ../

for param in `ls -a ./parameters | egrep '.json$'`
do
    NUMBER=`echo $param | tr -dc '0-9'`
    nr=$(($nr+1))
    build/softsim parameters/$param > results/outs/out_$NUMBER.txt &

    nr_runnings=`ps -e | grep "softsim" | wc -l`
    while [ $nr_runnings -ge 30 ]
    do
        sleep 3
        nr_runnings=`ps -e | grep "softsim" | wc -l`
    done
done
