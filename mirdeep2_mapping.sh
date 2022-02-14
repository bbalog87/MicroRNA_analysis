#!/bin/bash
## mirna detection pipeline with miRDeep2 package
grep -v ">" mature.miRNA.catalog.fa | sed 's/U/T/g' | sort | uniq | awk
'BEGIN{i=1}; {print ">seq_"i++; print}' > mirnas.fa

fastaparse.pl mirnas.fa -a 18 > mirnas_no_short.fa 2>mirnas_too_short

bowtie -p 1 -f -n 0 -e 80 -l 18 -a -m 5 \
       --best --strata zander.genome.for.miRNA.dectection \
	   --al mature.fa_mapped --un mature.fa_not_mapped \
	     mirnas_no_short.fa mappings.bwt
		 
# reads processed: 30600
# reads with at least one reported alignment: 2048 (6.69%)
# reads that failed to align: 28429 (92.91%)
# reads with alignments suppressed due to -m: 123 (0.40%)
#Reported 3352 alignments

convert_bowtie_output.pl mappings.bwt > mappings.arf


cat mappings.arf | sed 's/seq_//' | awk 'BEGIN{OFS="\t"}
{$1="001_"$1"_x1"; print}' > mappings_for_mirdeep2.arf
cat mature.fa_mapped | sed 's/>seq//' | awk '/_/{print ">001"$0"_x1";
next} {print}' > mature_for_mirdeep2.fa

###
miRDeep2.pl mature_for_mirdeep2.fa zander.genome.for.miRNA.dectection.fa \
            mappings_for_mirdeep2.arf none none none

# ended: 15:11:53
# total:0h:1m:11s
#
# fasta and bed files have been created in subfolder
#mirna_results_01_04_2019_t_15_02_15
#
#
# miRDeep runtime:
#
# started: 15:2:15
# ended: 15:11:53
# total:0h:9m:38s
#
#
#
# miRDeep runtime:
#
# started: 15:2:15
# ended: 15:11:53
# total:0h:9m:38s
#
