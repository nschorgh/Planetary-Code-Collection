#!/bin/bash

# merges and sorts output files from parallel runs of programs shadows and fieldofviews

/bin/rm -f horizons.dat 
touch horizons.dat
#/bin/rm -f fieldofviews.dat
#touch fieldofviews.dat


for i in `ls horizon.* | sort -V`;
do
    cat $i >> horizons.dat
done
#for i in `ls fieldofview.* | sort -V`;
#do
#    cat $i >> fieldofviews.dat
#done


# sort -V sorts numerically






