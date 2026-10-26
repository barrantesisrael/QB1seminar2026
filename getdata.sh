#!/bin/bash

# RNAseq data
wget -q https://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xf HBR_UHR_ERCC_ds_5pc.tar

# Chromosome 22 sequence incl. ERCC spike-ins
wget -q https://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

# Chromosome 22 annotations incl. ERCC spike-ins
wget -q https://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf

# adapt filenames for mapping
for f in *.read[12].fastq.gz; do
    new=$(echo "$f" | sed -E 's/^(HBR|UHR)_Rep([0-9]+)_.*\.read([12])\.fastq\.gz$/\1_\2_R\3.fq/')
    gunzip -c "$f" > "$new"
done
