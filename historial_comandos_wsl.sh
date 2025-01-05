## Eliminación de lecturas mapeadas duplicadas con Picard tools
ls BAM/*_sorted.bam | parallel -j 7 "picard MarkDuplicates I={} O=BAM/{/.}_unique.bam M=BAM/{/.}_metrics.txt REMOVE_DUPLICATES=true"

for file in BAM/*_trimmed30_filtered_sorted_unique.bam; do
	  mv "$file" "${file/_trimmed30_filtered_sorted_unique.bam/_unique.bam}"
  done # Para simplificar los nombres de los ficheros BAM

# Para ejecutar HaplotypeCaller de GATK, primero DEBEMOS añadir información
# sobre los grupos a los que pertenece cada read (@RG, SM, LB y PL) en los metadatos de cada BAM
# y volvemos a indexar el fichero BAM
for file in BAM/*_unique.bam; do
    basename=$(basename $file | cut -d'_' -f1)
    samtools addreplacerg -r "@RG\tID:${basename}\tSM:${basename}\tLB:${basename}\tPL:Ion_Torrent" $file -o BAM/${basename}_unique_RG.bam
    samtools index BAM/${basename}_unique_RG.bam
done

# GATK Recalibración de la calidad de las bases (paso de machine learning recomendado)
./gatk CreateSequenceDictionary -R ../../hg19.fa
./gatk samtools faidx ../../hg19.fa

wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
#wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz
wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
#wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz

mv 1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz.tbi
mv Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz.tbi

# Desde ~/Estudios_in_silico/programs/gatk-4.6.1.0
ls ../../BAM/*_unique_RG.bam | parallel -j 7 "./gatk BaseRecalibrator \
                                                    -R ../../hg19.fa \
                                                    -I {} \
                                                    --known-sites ../../Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
                                                    --known-sites ../../1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
                                                    -O ../../RECAL/{/.}_recal_data.table"

for file in ../../RECAL/*_recal_data.table; do
    basename=$(basename $file | cut -d'_' -f1)
    ./gatk ApplyBQSR -R ../../hg19.fa -I ../../BAM/${basename}_unique_RG.bam --bqsr-recal-file ${file} -O ../../RECAL/${basename}_finalrecal.bam
done

# GATK HaplotypeCaller pero con parallel (7 procesos a la vez SIN GVCF)
ls ../../RECAL/*_finalrecal.bam | parallel -j 7 "./gatk HaplotypeCaller -R ../../hg19.fa -I {} -O ../../VARIANT_CALLING/{/.}.vcf"

# GATK HaplotypeCaller pero con parallel (7 procesos a la vez con GVCF)
# ls ../../RECAL/*_finalrecal.bam | parallel -j 7 "./gatk HaplotypeCaller -R ../../hg19.fa -I {} -O ../../VARIANT_CALLING_GVCF/{/.}G.vcf -ERC GVCF"


#### SI LA ORDENACIÓN DE LAS COORDENADAS NO COINCIDE: VOLVER A ORDENAR LOS BAM CON PICARD SORTSAM
# EN LUGAR DE CON SAMTOOLS SORT :)

## ANNOVAR
# Conversión de los VCF a ficheros tipo input de ANNOVAR:
ls ../../VARIANT_CALLING/*.vcf | parallel -j 7 "./convert2annovar.pl -format vcf4 {} > ../../VARIANT_CALLING/{/.}.avinput"

# Anotación de las variantes con ANNOVAR (protocolos refGene para anotaciones genómicas,
# dbnsfp30a para predicciones funcionales, avsnp150 para variantes en dbSNP)
# operaciones g=genómica y f=frecuencia:
# Primero descargamos las bds en humandb/:
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/ # Para la anotación de genes y efectos funcionales
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/ # Para determinar el efecto de las variantes en proteínas (SIFT, PolyPhen-2, MutationTaster, LRT, FATHMM)
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/ # Para variantes ya reportadas (SNPs e INDELs)
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2012apr humandb/ # variantes genéticas de alta calidad del Proyecto 1000 genomas (fase 1, no la más reciente que es la 3 de 2015 porque los VCF
# se obtuvieron a partir de ficheros de SNPs e INDEL de la fase 1 (los unicos disponibles en internet))

#### MUESTRAS POR SEPARADO #####
ls ../../VARIANT_CALLING/*.avinput | parallel -j 7 "./table_annovar.pl {} humandb/ \
                                                            -buildver hg19 \
                                                            -out ../../ANNOVAR_ANNOTATED/{/.} \
                                                            -remove \
                                                            -protocol refGene,dbnsfp30a,avsnp150,ALL.sites.2012_04 \
                                                            -operation g,f,f,f \
                                                            -nastring . \
                                                            -csvout"

# Filtrado por exónicos (dentro de la carpeta de ANNOVAR_ANNOTATED)
ls | parallel -j 7 "grep -E "exonic" {} > {/.}_exonic.csv"

#### MUESTRAS JUNTAS (UN ÚNICO VCF VERSIÓN INPUT DE ANNOVAR CONCATENADO) ####
~/Estudios_in_silico/VARIANT_CALLING$ cat *finalrecal.avinput > allsamples.avinput

./table_annovar.pl ../../VARIANT_CALLING/allsamples.avinput humandb/ \
        -buildver hg19 \
        -out ../../ANNOVAR_ANNOTATED/allsamples \
        -remove \
        -protocol refGene,dbnsfp30a,avsnp150,ALL.sites.2012_04 \
        -operation g,f,f,f \
        -nastring . \
        -csvout

grep -E "exonic" allsamples.hg19_multianno.csv > allsamples.hg19_multianno_exonic.csv


## Filtrado del resultado de SNPs e INDELs anotado con ANNOVAR sobre allsamples (todas las muestras juntas) para solo mostrar los INDELS (nuevo fichero)
# En directorio main:

(base) patriciacvr01@LAPTOP-HMHETCR2:~/Estudios_in_silico$ awk -F ',' '$9 ~ /frameshift/ { print $0 }' ANNOVAR_ANNOTATED/allsamples.hg19_multianno_exonic.csv > ANNOVAR_ANNOTATED/allsamples.hg19_multianno_exonic_INDELS.csv # Supplementary Table 1
# Nos quedamos con aquellas líneas del csv original cuya columna 9 contenga la palabra "frameshift", que serán las deleciones (frameshift deletion) y las inserciones (frameshift insertion)

# Para sacar las estadísticas de la Figura 2 del artículo
(base) patriciacvr01@LAPTOP-HMHETCR2:~/Estudios_in_silico$ awk -F ',' -f numSNPSanalyzed_GATK.awk ANNOVAR_ANNOTATED/allsamples.hg19_multianno.csv > Figure2.txt

# Para sacar las estadísticas de la Figura 3 del artículo
(base) patriciacvr01@LAPTOP-HMHETCR2:~/Estudios_in_silico$ awk -F ',' -f distribnumSNPs_exonicfunction.awk ANNOVAR_ANNOTATED/allsamples.hg19_multianno_exonic.csv > Figure3.txt

# Para sacar las estadísticas de la Table 3 del artículo
ls ANNOVAR_ANNOTATED/ERR*_multianno_exonic.csv | parallel -j 7 "awk -F "," -f numSNPsINDELS_persample.awk {} > ANNOVAR_ANNOTATED/{/.}_TABLE3.txt"

# Non-synonymous SNPs, SIFTscore<0.05 y Polyphen2score entre 0.453–0.956  (No especifican si es PolyPhen2_HDIV_score o PolyPhen2_HVAR_score)
awk -F "," -f table4.awk ANNOVAR_ANNOTATED/allsamples.hg19_multianno_exonic.csv > ANNOVAR_ANNOTATED/allsamples.hg19_exonicSIFT_Polyphen2.csv


#### AHORA TOCA FILTRAR LAS MUTACIONES EN FUNCIÓN DEL SIFT SCORE Y POLYPHEN (Y NO SINÓNIMAS) Y QUEDARNOS CON ESAS
# CON ESAS MUTACIONES SON LAS QUE ENFRETAMOS A LA BASE DE DATOS gnomAD (+MAF) --> y sacar reports y gráficas.
# Una vez obtenidos los nsSNPs, comparamos frente a la base de datos de gnomAD y con filtrado MAF:
# Descargamos la base de datos gnomad211_exome para hg19 desde las bds de annovar
# (base de datos de gnomAD versión 2.1.1 del exoma humano hg19 que contiene las frecuencias
# de variantes genéticas en regiones CODIFICANTES)
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome humandb/
# Falta lo de MAF y la anotación de esos SNPs filtrados frente a tal base de datos gnomAD

# Una vez descargada la base de datos de anotación, para anotar las variantes genéticas, recordamos
# que primero debemos transformar el csv en un fichero de entrada de ANNOVAR (avinput)
awk -F"," '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' ANNOVAR_ANNOTATED/allsamples.hg19_exonicSIFT_Polyphen2.csv > ANNOVAR_ANNOTATED/allsamples_gnomAD.avinput

# Ahora realizamos la anotación como tal (desde el directorio con el programa annovar)
./table_annovar.pl ../../ANNOVAR_ANNOTATED/allsamples_gnomAD.avinput humandb/ \
        -buildver hg19 \
        -out ../../ANNOVAR_ANNOTATED/allsamples_gnomAD \
        -remove \
        -protocol refGene,gnomad211_exome \
        -operation g,f \
        -nastring . \
        -csvout
# table_annovar.pl devuelve ficheros con .hg19_multianno
head -n 1 ANNOVAR_ANNOTATED/allsamples_gnomAD.hg19_multianno.csv > ANNOVAR_ANNOTATED/novel_afr_variants.csv
awk -F',' '$16 == "." || $16 < 0.0001 {print $0}' ANNOVAR_ANNOTATED/allsamples_gnomAD.hg19_multianno.csv >> ANNOVAR_ANNOTATED/novel_afr_variants.csv
# variantes donde AF_afr no está disponible (representada como .) o tiene valores extremadamente bajos:

awk -F',' 'BEGIN{OFS=","} NR==1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' ANNOVAR_ANNOTATED/allsamples_gnomAD.hg19_multianno.csv > ANNOVAR_ANNOTATED/MAF_afr_filtered.csv
awk -F',' 'BEGIN{OFS=","} $16 != "." && $16 < 0.05 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16}' ANNOVAR_ANNOTATED/allsamples_gnomAD.hg19_multianno.csv >> ANNOVAR_ANNOTATED/MAF_afr_filtered.csv
# variantes donde AF_afr tiene un MAF<0.05 (las que se listan en la Table 4 del artículo)
### CHECKEAR SI SUPPLEMENTARY_TABLE2 Y SUPPLEMENTARY_TABLE3 se obtienen por filtrado de SIFTScore y Polyphen a cada muestra por separado o todas juntas.
