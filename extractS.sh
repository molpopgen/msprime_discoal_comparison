#!/bin/bash

Ne=10000
L=50000
POS=`echo "$L/2" | bc -l`

grep segsites msprime_neutral.ms | cut -d":" -f 2 > msprime_neutral.segsites
for i in 1 100 500 1000
do
    grep segsites msprime_$i.ms | cut -d":" -f 2 > msprime_$i.segsites
    grep segsites discoal_$i.ms | cut -d":" -f 2 > discoal_$i.segsites
done


