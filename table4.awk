## Filtrado del csv de exones de mutaciones Non-synonymous, SIFT-score<0.05 y Polyphen2-score entre 0.453-0.956 ##
# Sobre el fichero allsamples.hg19_multianno_exonic.csv

BEGIN {
    OFS="," # Define la coma como separador de campos en la salida (en lugar de espacios por defecto con print $0)
}
{
    gsub(/"/, "", $9); # Elimina las comillas en la columna 9
}
    ($9=="nonsynonymous SNV") && ($11 < 0.05) && ($15 >= 0.453 && $15 <= 0.956) { print $0 }
