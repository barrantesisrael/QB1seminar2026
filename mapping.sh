#!/bin/bash

# Cancer samples
for rep in 1 2 3; do
    salmon quant -i salmon_index -l A \
        -1 reads/UHR_${rep}_R1.fq -2 reads/UHR_${rep}_R2.fq \
        -o salmon_quant/UHR_rep${rep} -q
done

# Normal samples
for rep in 1 2 3; do
    salmon quant -i salmon_index -l A \
        -1 reads/HBR_${rep}_R1.fq -2 reads/HBR_${rep}_R2.fq \
        -o salmon_quant/HBR_rep${rep} -q
done