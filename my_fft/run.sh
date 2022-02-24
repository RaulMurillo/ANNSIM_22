#!/bin/bash

cd build

for type in double float
do
    sed -i "s/using Real.*/using Real = ${type};/g" ../src/test.cpp
    make
    ./my_fft > ../results/${type}.log
done

for es in 8 5
do
    sed -i "s/using Real.*/using Real = cfloat<16,${es}>;/g" ../src/test.cpp
    make
    ./my_fft > ../results/float_16_${es}.log
done

for n in 32 28 24 20 16 14 12 10 8
do
    sed -i "s/using Real.*/using Real = posit<${n},2>;/g" ../src/test.cpp
    make
    ./my_fft > ../results/posit_${n}.log
done

cd ..