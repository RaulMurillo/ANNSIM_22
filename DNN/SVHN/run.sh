#!/bin/bash

mkdir -p P8

sed -i "s/num_format = .*/num_format = bit_8/g" svhn_posit.py

python svhn_posit.py

mv train.hist P8
mv test.hist P8

#######
mkdir -p P10

sed -i "s/num_format = .*/num_format = bit_10/g" svhn_posit.py

python svhn_posit.py

mv train.hist P10
mv test.hist P10

#######
mkdir -p P12

sed -i "s/num_format = .*/num_format = bit_12/g" svhn_posit.py

python svhn_posit.py

mv train.hist P12
mv test.hist P12

#######
mkdir -p P14

sed -i "s/num_format = .*/num_format = bit_14/g" svhn_posit.py

python svhn_posit.py

mv train.hist P14
mv test.hist P14

#######
mkdir -p P16

sed -i "s/num_format = .*/num_format = bit_16/g" svhn_posit.py

python svhn_posit.py

mv train.hist P16
mv test.hist P16

#######
mkdir -p FP16

sed -i "s/num_format = .*/num_format = IEEE_Half/g" svhn_posit.py

python svhn_posit.py

mv train.hist FP16
mv test.hist FP16

#######
mkdir -p bfloat16

sed -i "s/num_format = .*/num_format = bfloat16/g" svhn_posit.py

python svhn_posit.py

mv train.hist bfloat16
mv test.hist bfloat16

