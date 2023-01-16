#!/bin/bash

# First download a transcriptome reference from Essemble called CF_002263795.2_ARS-UCD1.3_rna.fna

#Index with salmon
salmon index -t GCF_002263795.2_ARS-UCD1.3_rna.fna -i salmon_index

# run the for loop to get counts


for fn in reads/SRR131070{18..23};
do
samp=${fn}
echo "Processing sample ${samp}"
salmon quant -i salmon_index -l A \
         -1 ${samp}_1.fastq \
         -2 ${samp}_2.fastq \
         -p 8 --validateMappings -o quants/${samp}_quant
done 
