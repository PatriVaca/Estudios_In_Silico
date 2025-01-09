#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##########################################
#         PIPELINE MASTERBIO             #
# Autora: Patricia del Carmen Vaca Rubio #
##########################################

################################################################################
# La parte inicial de este trabajo (desde la descarga de los ficheros FASTQ hasta
# la obtención de las lecturas mapeadas en formato BAM) se realizó en masterbio
# (computador para los alumnos del máster en Bioinformática de la Universitat de
# València), dado que cuenta con más memoria RAM y cores que mi computador,
# facilitando así la descarga y el procesamiento de datos de tales dimensiones.
################################################################################

### CONFIGURACIÓN PREVIA ###
set -euo pipefail # detiene el script si hay errores o comandos no definidos

# Ejecutando todo desde el directorio raíz: /home/vapadel/EINSILICO/ (MASTERBIO), en este caso ##
## Crear directorios de salida
mkdir -p FASTQ
mkdir -p SRA
mkdir -p FASTQC
mkdir -p FASTQ_Filtered
mkdir -p FASTQC_Filtered
mkdir -p SAM
mkdir -p BAM

## Rutas base
BASE_DIR=$(pwd) # Directorio donde se ejecuta el script
FASTQ_PATH="$BASE_DIR/FASTQ"
SRA_PATH="$BASE_DIR/SRA"
FASTQC_PATH="$BASE_DIR/FASTQC"
FASTQ_FILTERED_PATH="$BASE_DIR/FASTQ_Filtered"
FASTQC_FILTERED_PATH="$BASE_DIR/FASTQC_Filtered"
SAM_PATH="$BASE_DIR/SAM"
BAM_PATH="$BASE_DIR/BAM"
N=12 # Número de procesos que se ejecutarán en paralelo en los procesos
# paralelizables (en función del computador). Cuando j -7, un proceso para cada
# una de las 7 muestras del estudio

## Lista de muestras a descargar (puedes agregar más u otros IDs aquí)
SRA_IDS=(
  "ERR3013393"
  "ERR3013392"
  "ERR3013391"
  "ERR3013389"
  "ERR3013388"
  "ERR3013387"
  "ERR3013386"
)

### DESCARGA DE LAS MUESTRAS ###
## Descarga de muestras
echo "Descargando muestras de SRA..."
for SRA_ID in "${SRA_IDS[@]}"; do
  prefetch -p "$SRA_ID" -O "$SRA_PATH/" --resume no # -p muestra progreso
done

## Conversión de SRA a FASTQ
echo "Convirtiendo archivos SRA a FASTQ..."
for SRA_ID in "${SRA_IDS[@]}"; do
  fasterq-dump -t $N -p -O "$FASTQ_PATH/" "$SRA_PATH/$SRA_ID/" --force
done

### INFORME DE CALIDAD PREVIO AL PROCESADO DE LAS LECTURAS ###
## Generación de informes FASTQC
echo "Generando informes FASTQC..."
fastqc -t $N "$FASTQ_PATH"/*.fastq -o "$FASTQC_PATH/"

### FILTRADO DE CALIDAD DE LAS LECTURAS ###
echo "Filtrando lecturas por calidad con cutadapt..."
ls "$FASTQ_PATH"/*.fastq | parallel -j $N "cutadapt -q 30 --minimum-length 50 --maximum-length 300 -o $FASTQ_FILTERED_PATH/{/.}_trimmed30.fastq {}" # Filtrado por calidad > 30 de los FASTQ y longitud mínima de las reads finales de 50 y máxima de 300 bases (estándar para lecturas largas)

### Generación de informes FASTQC post-filtrado ###
echo "Generando informes FASTQC después del filtrado..."
fastqc -t $N "$FASTQ_FILTERED_PATH"/*.fastq -o "$FASTQC_FILTERED_PATH/"

### MAPEO DE LAS LECTURAS SOBRE EL GENOMA DE REFERENCIA hg19 DE UCSC (CORRESPONDIENTE A GRCH37) ###
echo "Descargando genoma de referencia hg19..."
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O "$BASE_DIR/hg19.fa.gz"
gunzip "$BASE_DIR/hg19.fa.gz"

## Indexado del genoma de referencia
echo "Indexando el genoma de referencia..."
bwa index "$BASE_DIR/hg19.fa"

## Mapeo de lecturas con bwa mem
echo "Mapeando lecturas contra el genoma de referencia..."
ls "$FASTQ_FILTERED_PATH"/*.fastq | parallel -j 7 "bwa mem -t 2 $BASE_DIR/hg19.fa {} > $SAM_PATH/{/.}.sam"
# Mapeo de las reads frente al genoma de referencia (ya indexado) con 7 procesos, 2 hilos cada uno (7x2=14 hilos lógicos en total de 16 que tiene el sistema)

## Filtrado de mapeos con samtools
echo "Filtrando mapeos con samtools..."
ls "$SAM_PATH"/*.sam | parallel -j 7 "samtools view -h {} | awk '\$0 ~ /^@/ || \$0 ~ /NM:i:[0-2]/' > $SAM_PATH/{/.}_filtered.sam"
# Una vez mapeadas las lecturas, nos quedamos solo con aquellas lecturas del .sam que tengan como máximo 2 mismatches con respecto al genoma de referencia, como indican en el artículo

## Conversión de SAM a BAM
echo "Convirtiendo archivos SAM a BAM..."
ls "$SAM_PATH"/*_filtered.sam | parallel -j 7 "samtools view -S -b {} > $BAM_PATH/{/.}.bam"

## Ordenando los archivos BAM
echo "Ordenando archivos BAM..."
ls "$BAM_PATH"/*.bam | parallel -j 7 "samtools sort {} -o $BAM_PATH/{/.}_sorted.bam"

## Indexado de los archivos BAM ordenados
echo "Indexando archivos BAM ordenados..."
ls "$BAM_PATH"/*_sorted.bam | parallel -j 7 "samtools index {}"

### Mensaje final ###
echo "Pipeline completado. Revisa los resultados en:"
echo "  - FASTQ: $FASTQ_PATH"
echo "  - SRA: $SRA_PATH"
echo "  - FASTQC: $FASTQC_PATH"
echo "  - FASTQ Filtered: $FASTQ_FILTERED_PATH"
echo "  - SAM: $SAM_PATH"
echo "  - BAM: $BAM_PATH"

# Una vez que hemos obtenido los ficheros BAM (y sus ficheros derivados),
# muevo dichos ficheros hasta mi computador para seguir con el análisis.
