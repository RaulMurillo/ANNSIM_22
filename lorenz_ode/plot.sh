#!/bin/bash

PICDIR=pics2
mkdir -p ${PICDIR}

cd build

for format in double float float16 bfloat posit_32_2 posit_28_2 posit_24_2 posit_20_2 posit_18_2 posit_16_2 posit_14_2 posit_12_2 posit_10_2 posit_8_2
do
    echo ${format}
    # gnuplot < lorenz_ode_commands_${format}.txt
    # mv xyz_time.png ../pics/${format}_time.png
    # mv xyz_3d.png ../pics/${format}_3d.png
    python ../plot.py "lorenz_ode_data_${format}.txt" "${format}"
    mv *.pdf ../${PICDIR}/
done