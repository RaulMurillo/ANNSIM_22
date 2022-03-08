#!/bin/bash

mkdir -p P8

sed -i "s/num_format = .*/num_format = bit_8/g" cifar10_posit.py

python cifar10_posit.py

mv train.hist P8
mv test.hist P8

#######
mkdir -p P10

sed -i "s/num_format = .*/num_format = bit_10/g" cifar10_posit.py

python cifar10_posit.py

mv train.hist P10
mv test.hist P10

#######
mkdir -p FP16

sed -i "s/num_format = .*/num_format = IEEE_Half/g" cifar10_posit.py

python cifar10_posit.py

mv train.hist FP16
mv test.hist FP16

#######
mkdir -p bfloat16

sed -i "s/num_format = .*/num_format = bfloat16/g" cifar10_posit.py

python cifar10_posit.py

mv train.hist bfloat16
mv test.hist bfloat16

