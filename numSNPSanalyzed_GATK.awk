## Obtención de las estadísticas de la Figura 2 del artículo ##

BEGIN {
    downstream=0;
    exonic=0;
    exonic_splicing=0;
    intergenic=0;
    intronic=0;
    ncRNA_exonic=0;
    ncRNA_intronic=0;
    ncRNA_splicing=0;
    splicing=0;
    upstream=0;
    upstream_downstream=0;
    utr3=0;
    utr5=0;
}

{
    gsub(/"/, "", $6); # Elimina las comillas en la columna 6
}

$6 == "downstream" { downstream += 1 }
$6 == "exonic" { exonic += 1 }
$6 == "exonic_splicing" { exonic_splicing += 1 }
$6 == "intergenic" { intergenic += 1 }
$6 == "intronic" { intronic += 1 }
$6 == "ncRNA_exonic" { ncRNA_exonic += 1 }
$6 == "ncRNA_intronic" { ncRNA_intronic += 1 }
$6 == "ncRNA_splicing" { ncRNA_splicing += 1 }
$6 == "splicing" { splicing += 1 }
$6 == "upstream" { upstream += 1 }
$6 == "upstream:downstream" { upstream_downstream += 1 }
$6 == "UTR3" { utr3 += 1 }
$6 == "UTR5" { utr5 += 1 }

END {
    print "Downstream: " downstream;
    print "Exonic: " exonic;
    print "Exonic Splicing: " exonic_splicing;
    print "Intergenic: " intergenic;
    print "Intronic: " intronic;
    print "ncRNA Exonic: " ncRNA_exonic;
    print "ncRNA Intronic: " ncRNA_intronic;
    print "ncRNA Splicing: " ncRNA_splicing;
    print "Splicing: " splicing;
    print "Upstream: " upstream;
    print "Upstream:Downstream: " upstream_downstream;
    print "UTR3: " utr3;
    print "UTR5: " utr5;
}


