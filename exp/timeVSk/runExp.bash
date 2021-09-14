#!/bin/bash

n=100000;
nLayer=5
#p=0.0008
p=`echo "scale=8;5/$n" | bc`
echo p=$p

mdir='mult/'
ov=0.1
k=100;
alg='-K'
x=20;

for nLayer in `seq 5 1 20`;
do
    mkdir $mdir
    synth $n $p $nLayer $ov $mdir
    mim -M $mdir $alg -k $k -x $x
    rm -rf $mdir
done

