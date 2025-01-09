## Para obtener la Table 3 del artículo ##

BEGIN {
		frameshift_deletion=0;
		frameshift_insertion=0;
		nonframeshift_deletion=0;
		nonframeshift_insertion=0;
		non_synonymous=0;
		synonymous=0;
		total_lines=0;
	}
{	
	total_lines+=1
    gsub(/"/, "", $9); # Elimina las comillas en la columna 9
}	

	$9=="frameshift deletion" { frameshift_deletion += 1 }
	$9=="frameshift insertion" { frameshift_insertion += 1 }
	$9=="nonframeshift deletion" { nonframeshift_deletion += 1 }
	$9=="nonframeshift insertion" { nonframeshift_insertion += 1 }
	$9=="nonsynonymous SNV" { non_synonymous += 1 }
	$9=="synonymous SNV" { synonymous += 1 }

END {
		print "frameshift deletion: " frameshift_deletion;
		print "frameshift insertion: " frameshift_insertion;
		print "nonframeshift deletion: " nonframeshift_deletion;
		print "nonframeshift insertion: " nonframeshift_insertion;
		print "nonsynonymous SNP: " non_synonymous;
		print "synonymous SNP: " synonymous;
		print "total SNPs and INDELs: " total_lines;
	} 