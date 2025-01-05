#!/bin/bash

# Comando de paralelización
parallel -j 7 "gatk HaplotypeCaller -R hg19.fa -I BAM/{1} -O VARIANT_CALLING/{1}.vcf" ::: ls BAM/*_unique.bam
