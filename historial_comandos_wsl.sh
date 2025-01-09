#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##########################################
#         PIPELINE COMPUTADOR            #
# Autora: Patricia del Carmen Vaca Rubio #
##########################################

################################################################################
# La parte final de este trabajo (a partir de los ficheros BAM generados en el 
# supercomputador de masterbio) se realizó en mi propio computador para poder
# utilizar programas de forma local para cuya instalación se requerían permisos
# de superusuario, como gatk o annovar.
################################################################################

# PASO PREVIO A ESTA PIPELINE: CREAR UN DIRECTORIO BAM QUE YA CONTENGA LOS FICHEROS
# GENERADOS POR EL SCRIPT ANTERIOR: 'historial_comandos_masterbio.sh'

### CONFIGURACIÓN PREVIA ###
set -euo pipefail # detiene el script si hay errores o comandos no definidos

## Ejecutando todo desde el directorio raíz: ~/Estudios_in_silico/, en este caso ###
# mkdir BAM # YA CONTIENE LOS FICHEROS BAM TRAÍDOS DESDE MASTERBIO
mkdir -p RECAL
mkdir -p VARIANT_CALLING
mkdir -p PROGRAMS # Contiene aquellos programas que usamos en el análisis y que no
# están disponibles en los repositorios estándar de Ubuntu o Debian
# (no instalables por sudo apt-get install), como gatk o ANNOVAR
mkdir -p ANNOVAR_ANNOTATED

## Rutas base
BASE_DIR=$(pwd) # Directorio donde se ejecuta el script
BAM_PATH="$BASE_DIR/BAM"
RECAL_PATH="$BASE_DIR/RECAL"
VARIANT_CALLING_PATH="$BASE_DIR/VARIANT_CALLING"
pGATK="$BASE_DIR/PROGRAMS/gatk-4.6.1.0"
pANNOVAR_PATH="$BASE_DIR/PROGRAMS/annovar"
ANNOVAR_ANNOTATED_PATH="$BASE_DIR/ANNOVAR_ANNOTATED"
HUMANDB_PATH="$pANNOVAR_PATH/humandb"

### Eliminación de lecturas mapeadas duplicadas con Picard tools ###
echo "*******************************************************"
echo "Eliminando duplicados de los archivos BAM..."
echo "*******************************************************"
ls "$BAM_PATH"/*_sorted.bam | parallel -j 7 "picard MarkDuplicates I={} O=$BAM_PATH/{/.}_unique.bam M=$BAM_PATH/{/.}_metrics.txt REMOVE_DUPLICATES=true"

for file in "$BAM_PATH"/*_trimmed30_filtered_sorted_unique.bam; do
	mv "$file" "${file/_trimmed30_filtered_sorted_unique.bam/_unique.bam}"
done # Para simplificar los nombres de los ficheros BAM

### Añadiendo metadatos y reindexando BAM ###
# Para ejecutar HaplotypeCaller de GATK, primero DEBEMOS añadir información
# sobre los grupos a los que pertenece cada read (@RG, SM, LB y PL) en los
# metadatos de cada BAM y volvemos a indexar el fichero BAM
for file in "$BAM_PATH"/*_unique.bam; do
    basename=$(basename $file | cut -d'_' -f1)
    samtools addreplacerg -r "@RG\tID:${basename}\tSM:${basename}\tLB:${basename}\tPL:Ion_Torrent" "$file" -o "$BAM_PATH/${basename}_unique_RG.bam"
    samtools index "$BAM_PATH/${basename}_unique_RG.bam" # Hay que indexar los BAM cada vez que se procesan por una nueva herramienta (puesto que cambian sus índices),
# como con Picard MarkDuplicates, que cambia la estructura del BAM o con samtools addreplacerg, que cambia los metadatos asociados a cada read
done

### Recalibración de calidad de bases con GATK ###
echo "*******************************************************"
echo "Descargando genoma de referencia hg19..."
echo "*******************************************************"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O "$BASE_DIR/hg19.fa.gz"
gunzip "$BASE_DIR/hg19.fa.gz"

## Indexado del genoma de referencia
echo "*******************************************************"
echo "Indexando el genoma de referencia..."
echo "*******************************************************"
bwa index "$BASE_DIR/hg19.fa"

## GATK Recalibración de la calidad de las bases
# (paso de machine learning previo al llamado de variantes recomendado)
echo "*******************************************************"
echo "Descargando y creando archivos necesarios para la recalibración..."
echo "*******************************************************"
"$pGATK/./gatk" CreateSequenceDictionary -R hg19.fa
samtools faidx hg19.fa
wget -P "$BASE_DIR" http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget -P "$BASE_DIR" http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz

# Para solventar una serie de problemas con los zip descargados
gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
bgzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz


gunzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
bgzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf
tabix -p vcf 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz

echo "*******************************************************"
echo "Iniciando recalibración de bases con BaseRecalibrator..."
echo "*******************************************************"
ls "$BAM_PATH"/*_unique_RG.bam | parallel -j 7 "$pGATK/./gatk BaseRecalibrator \
                                                    -R $BASE_DIR/hg19.fa \
                                                    -I {} \
                                                    --known-sites $BASE_DIR/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
                                                    --known-sites $BASE_DIR/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
                                                    -O $RECAL_PATH/{/.}_recal_data.table"

for file in "$RECAL_PATH"/*_recal_data.table; do
    basename=$(basename $file | cut -d'_' -f1)
    "$pGATK/./gatk" ApplyBQSR -R "$BASE_DIR/hg19.fa" -I "$BAM_PATH/${basename}_unique_RG.bam" --bqsr-recal-file ${file} -O "$RECAL_PATH/${basename}_finalrecal.bam"
done

### Llamado de variantes con HaplotypeCaller ###
echo "*******************************************************"
echo "Ejecutando HaplotypeCaller para llamado de variantes..."
echo "*******************************************************"

# GATK HaplotypeCaller con parallel (7 procesos a la vez)
ls "$RECAL_PATH"/*_finalrecal.bam | parallel -j 7 "$pGATK/./gatk HaplotypeCaller -R $BASE_DIR/hg19.fa -I {} -O $VARIANT_CALLING_PATH/{/.}.vcf"

### PROCESO DE ANOTACIÓN DE LAS VARIANTES GENÉTICAS CON ANNOVAR ###
## Conversión de los VCF a ficheros tipo input de ANNOVAR:
ls "$VARIANT_CALLING_PATH"/*.vcf | parallel -j 7 "$pANNOVAR_PATH/./convert2annovar.pl -format vcf4 {} > $VARIANT_CALLING_PATH/{/.}.avinput"

## Descarga de bases de datos necesarias para ANNOVAR
# Anotación de las variantes con ANNOVAR (protocolos refGene para anotaciones
# genómicas, dbnsfp30a para predicciones funcionales, avsnp150 para variantes en
# dbSNP), operaciones g=genómica y f=frecuencia)
echo "*******************************************************"
echo "Descargando bases de datos para ANNOVAR..."
echo "*******************************************************"
"$pANNOVAR_PATH/./annotate_variation.pl" -buildver hg19 -downdb -webfrom annovar refGene "$HUMANDB_PATH/" # Para la anotación de genes y efectos funcionales
"$pANNOVAR_PATH/./annotate_variation.pl" -buildver hg19 -downdb -webfrom annovar dbnsfp30a "$HUMANDB_PATH/" # Para determinar el efecto de las variantes en proteínas (SIFT, PolyPhen-2, MutationTaster, LRT, FATHMM)
"$pANNOVAR_PATH/./annotate_variation.pl" -buildver hg19 -downdb -webfrom annovar avsnp150 "$HUMANDB_PATH/" # Para variantes ya reportadas (SNPs e INDELs)
"$pANNOVAR_PATH/./annotate_variation.pl" -buildver hg19 -downdb -webfrom annovar 1000g2012apr "$HUMANDB_PATH/"
# variantes genéticas de alta calidad del Proyecto 1000 genomas (fase 1, no la más reciente que es la 3 de 2015 porque los VCF se obtuvieron a partir de ficheros de SNPs e INDEL de la fase 1 (los unicos disponibles en internet))

#----------MUESTRAS POR SEPARADO----------
## Anotación de variantes por muestras
echo "*******************************************************"
echo "Anotando variantes con ANNOVAR por muestra..."
echo "*******************************************************"
ls "$VARIANT_CALLING_PATH"/*.avinput | parallel -j 7 "$pANNOVAR_PATH/./table_annovar.pl {} $HUMANDB_PATH/ \
                                                            -buildver hg19 \
                                                            -out $ANNOVAR_ANNOTATED_PATH/{/.} \
                                                            -remove \
                                                            -protocol refGene,dbnsfp30a,avsnp150,ALL.sites.2012_04 \
                                                            -operation g,f,f,f \
                                                            -nastring . \
                                                            -csvout"

## Filtrado de variantes exónicas por muestra
echo "*******************************************************"
echo "Filtrando variantes exónicas por muestra..."
echo "*******************************************************"
ls "$ANNOVAR_ANNOTATED_PATH"/*.csv | parallel -j 7 "grep -E 'exonic' {} > $ANNOVAR_ANNOTATED_PATH/{/.}_exonic.csv"

#----------MUESTRAS JUNTAS----------
## Anotación y análisis conjunto de todas las muestras
# Un único VCF con formato de input de annovar concatenado
echo "*******************************************************"
echo "Concatenando y anotando todas las muestras juntas..."
echo "*******************************************************"
cat "$VARIANT_CALLING_PATH"/*finalrecal.avinput > "$VARIANT_CALLING_PATH/allsamples.avinput"

"$pANNOVAR_PATH/./table_annovar.pl" "$VARIANT_CALLING_PATH/allsamples.avinput" "$HUMANDB_PATH/" \
    -buildver hg19 \
    -out "$ANNOVAR_ANNOTATED_PATH/allsamples" \
    -remove \
    -protocol refGene,dbnsfp30a,avsnp150,ALL.sites.2012_04 \
    -operation g,f,f,f \
    -nastring . \
    -csvout

echo "*******************************************************"
echo "Filtrando variantes exónicas de todas las muestras juntas..."
echo "*******************************************************"
grep -E "exonic" "$ANNOVAR_ANNOTATED_PATH/allsamples.hg19_multianno.csv" > "$ANNOVAR_ANNOTATED_PATH/allsamples.hg19_multianno_exonic.csv"

## Filtrado del resultado de SNPs e INDELs anotado con ANNOVAR sobre allsamples (todas las muestras juntas) para solo mostrar los INDELS (nuevo fichero)
## Filtrado por INDELs
echo "*******************************************************"
echo "Filtrando INDELs de las variantes anotadas..."
echo "*******************************************************"
awk -F ',' '$9 ~ /frameshift/ { print $0 }' "$ANNOVAR_ANNOTATED_PATH/allsamples.hg19_multianno_exonic.csv" > "$ANNOVAR_ANNOTATED_PATH/allsamples.hg19_multianno_exonic_INDELS.csv" # Supplementary Table 1
# Nos quedamos con aquellas líneas del csv original cuya columna 9 contenga la palabra "frameshift", que serán las deleciones (frameshift deletion) y las inserciones (frameshift insertion)

## Generación de estadísticas
# Para sacar las estadísticas de la Figura 2 del artículo
echo "*******************************************************"
echo "Generando estadísticas para la Figura 2..."
echo "*******************************************************"
awk -F ',' -f "$BASE_DIR/numSNPSanalyzed_GATK.awk" "$ANNOVAR_ANNOTATED_PATH/allsamples.hg19_multianno.csv" > Figure2.txt

# Para sacar las estadísticas de la Figura 3 del artículo
awk -F ',' -f "$BASE_DIR/distribnumSNPs_exonicfunction.awk" "$ANNOVAR_ANNOTATED_PATH/allsamples.hg19_multianno_exonic.csv" > Figure3.txt

# Para sacar las estadísticas de la Tabla 3 del artículo
echo "*******************************************************"
echo "Generando estadísticas para la Table 3..."
echo "*******************************************************"
ls "$ANNOVAR_ANNOTATED_PATH"/ERR*_multianno_exonic.csv | parallel -j 7 "awk -F ',' -f $BASE_DIR/numSNPsINDELS_persample.awk {} > $ANNOVAR_ANNOTATED_PATH/{/.}_TABLE3.txt"

## Filtrado para SNPs no sinónimos y análisis adicional con gnomAD
# Non-synonymous SNPs, SIFTscore<0.05 y Polyphen2score entre 0.453–0.956  (No especifican si es PolyPhen2_HDIV_score o PolyPhen2_HVAR_score)
awk -F "," -f table4.awk ANNOVAR_ANNOTATED/allsamples.hg19_multianno_exonic.csv > ANNOVAR_ANNOTATED/allsamples.hg19_exonicSIFT_Polyphen2.csv

# Una vez obtenidos los nsSNPs, comparamos frente a la base de datos de gnomAD y con filtrado MAF.
# Descargamos la base de datos gnomad211_exome para hg19 desde las bds de annovar
# (base de datos de gnomAD versión 2.1.1 del exoma humano hg19 que contiene las frecuencias
# de variantes genéticas en regiones CODIFICANTES)
echo "*******************************************************"
echo "Descargando base de datos gnomAD..."
echo "*******************************************************"
"$pANNOVAR_PATH/./annotate_variation.pl" -buildver hg19 -downdb -webfrom annovar gnomad211_exome "$HUMANDB_PATH/"

# Una vez descargada la base de datos de anotación, para anotar las variantes genéticas, recordamos
# que primero debemos transformar el csv en un fichero de entrada de ANNOVAR (avinput)
echo "*******************************************************"
echo "Preparando entrada para anotación de gnomAD..."
echo "*******************************************************"
awk -F"," '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' "$ANNOVAR_ANNOTATED_PATH/allsamples.hg19_exonicSIFT_Polyphen2.csv" > "$ANNOVAR_ANNOTATED_PATH/allsamples_gnomAD.avinput"

echo "*******************************************************"
echo "Anotando variantes con gnomAD..."
echo "*******************************************************"
"$pANNOVAR_PATH/./table_annovar.pl" "$ANNOVAR_ANNOTATED_PATH/allsamples_gnomAD.avinput" "$HUMANDB_PATH/" \
    -buildver hg19 \
    -out "$ANNOVAR_ANNOTATED_PATH/allsamples_gnomAD" \
    -remove \
    -protocol refGene,gnomad211_exome \
    -operation g,f \
    -nastring . \
    -csvout

# table_annovar.pl devuelve ficheros con .hg19_multianno
echo "*******************************************************"
echo "Filtrando variantes basadas en MAF..."
echo "*******************************************************"
head -n 1 "$ANNOVAR_ANNOTATED_PATH/allsamples_gnomAD.hg19_multianno.csv" > "$ANNOVAR_ANNOTATED_PATH/novel_afr_variants.csv"
awk -F',' '$16 == "." || $16 < 0.0001 {print $0}' "$ANNOVAR_ANNOTATED_PATH/allsamples_gnomAD.hg19_multianno.csv" >> "$ANNOVAR_ANNOTATED_PATH/novel_afr_variants.csv"
# variantes donde AF_afr no está disponible (representada como .) o tiene valores extremadamente bajos

awk -F',' 'BEGIN{OFS=","} NR==1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' "$ANNOVAR_ANNOTATED_PATH/allsamples_gnomAD.hg19_multianno.csv" > "$ANNOVAR_ANNOTATED_PATH/MAF_afr_filtered.csv"
awk -F',' 'BEGIN{OFS=","} $16 != "." && $16 < 0.06 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16}' "$ANNOVAR_ANNOTATED_PATH/allsamples_gnomAD.hg19_multianno.csv" >> "$ANNOVAR_ANNOTATED_PATH/MAF_afr_filtered.csv"

echo "### PROCESO ANNOVAR COMPLETADO ###"
# variantes donde AF_afr tiene un MAF<0.05 (las que se listan en la Tabla 4 del artículo)
### CHECKEAR SI SUPPLEMENTARY_TABLE2 Y SUPPLEMENTARY_TABLE3 se obtienen por filtrado de SIFTScore y Polyphen a cada muestra por separado o todas juntas.
