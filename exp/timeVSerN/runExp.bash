#!/bin/bash

n=1024;
nLayer=10
#p=0.0008
p=`echo "scale=8;5/$n" | bc`
echo p=$p

mdir='mult/'
ov=0.1
k=100;
alg='-K'
x=20;

for nr in `seq 1 10`;
do
    mkdir $mdir
    echo n=$n
    synth $n $p $nLayer $ov $mdir
    mim -M $mdir $alg -k $k -x $x
    rm -rf $mdir
    n=`echo "scale=8;2*$n" | bc`
    p=`echo "scale=8;5/$n" | bc`
done

