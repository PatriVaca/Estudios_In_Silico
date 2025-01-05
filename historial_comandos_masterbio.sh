prefetch -p ERR3013393 ERR3013392 ERR3013391 ERR3013389 ERR3013388 ERR3013387 ERR3013386 # -p para mostrar el progreso de descarga
fasterq-dump -t 12 -p -O FASTQ/ SRA/sra/ERR3013393.sra SRA/sra/ERR3013392.sra SRA/sra/ERR3013391.sra SRA/sra/ERR3013389.sra SRA/sra/ERR3013388.sra SRA/sra/ERR3013387.sra SRA/sra/ERR3013386.sr
#for file in *.fastq; do fastqc "$file" -o ~/EINSILICO/FASTQC; done # Ejecutándolo desde el interior de la carpeta FASTQ con todos los .fastq
# Otra opción para hacer los FASTQC (más rápida y automatizada)
fastqc -t 12 ~/EINSILICO/FASTQ/*.fastq -O ~/EINSILICO/FASTQC

ls FASTQ/*.fastq | parallel -j 12 'cutadapt -q 30 --minimum-length 50 --maximum-length 300 -o FASTQ_Filtered/{/.}_trimmed30.fastq {}' # Filtrado por calidad > 30 de los FASTQ y longitud mínima de las reads finales de 50 y máxima de 300 bases (estándar para lecturas largas)

## FASTQC reports de nuevo post-filtrado de las reads
fastqc -t 12 ~/EINSILICO/FASTQ_Filtered/*.fastq -O ~/EINSILICO/FASTQC_Filtered

## MAPEO DE LAS READS SOBRE EL GENOMA DE REFERENCIA hg19 DE UCSC(correspondiente a GRCH37) ##

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

bwa index hg19.fa # 1er paso: indexamos el genoma de referencia

ls FASTQ_Filtered/*.fastq | parallel -j 7 "bwa mem -t 2 hg19.fa {} > SAM/{/.}.sam" # Mapeo de las reads frente al genoma de referencia (ya indexado) con 7 procesos, 2 hilos cada uno (7x2=14 hilos lógicos en total de 16 que tiene el sistema)

ls SAM/*.sam | parallel -j 7 "samtools view -h {} | awk '\$0 ~ /^@/ || \$0 ~ /NM:i:[0-2]/' > SAM/{/.}_filtered.sam" # Una vez mapeadas las lecturas, nos quedamos solo con aquellas lecturas del .sam que tengan como máximo 2 mismatches con respecto al genoma de referencia, como indican en el artículo

ls SAM/*_filtered.sam | parallel -j 7 "samtools view -S -b {} > BAM/{/.}.bam"

ls BAM/*.bam | parallel -j 7 "samtools sort {} -o BAM/{/.}_sorted.bam"

ls BAM/*_sorted.bam | parallel -j 7 "samtools index {}" # Hay que indexar los BAM cada vez que se procesan por una nueva herramienta (puesto que cambian sus índices),
# como con Picard MarkDuplicates, que cambia la estructura del BAM o con samtools addreplacerg, que cambia los metadatos asociados a cada read
