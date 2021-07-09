#!/bin/bash

Ne=10000
L=50000
POS=`echo "$L/2" | bc -l`
SWEEPSITE=0`echo "$POS / $L" | bc -l` 
NREPS=100

# echo "python3 run_msprime.py neutral --nsam 50 --nreps $NREPS --theta 1000 --rho 1000 --length $L --seed 1234512 --prefix msprime_neutral"
for i in 1 100 500 1000
do
    freq=0`echo "$i / 2 / $Ne" | bc -l`
    # echo "python3 run_msprime.py soft --nsam 50 --nreps $NREPS --theta 1000 --rho 1000 --length $L --Ne $Ne --alpha 1000 --position $POS --freq $freq --seed 1236 --prefix msprime_$i"
    # echo "./discoal 50 $NREPS $L -ws 0.0 -t 1000.0 -r 1000.0 -x $SWEEPSITE -a 1000 -f $freq -N $Ne > discoal_$i.ms"
    echo "PYTHONPATH=../fwdpy11 python3 fp11sweeps.py $i $NREPS $RANDOM > fp11_$i.segsites"
done

