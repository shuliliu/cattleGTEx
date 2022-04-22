cd ~/ASE/Adipose

    vcf_dir="~/Adipose_new/permutation"
################################################################
### at least 10 individuals with ASE data and a minimum of 8 reads per individual ( don't know how to filter; maybe ignore-----default parameter: --min-cov 8 in phaser_expr_matrix.py).
### aFC means var_het_afc? ( assume so)

###All the nominal-sig sQTLs
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05) print $1,$2,$8}' OFS="\t" Adipose_results.txt | awk 'FNR==NR{a[$1"\t"$2]=$5; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Adipose.nominals.2rd.sig.txt - > Adipose_slope.ASE.correlation.txt


####Only top significant sQTLs (FDR-storey < 0.05) per gene
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05 ) print $1,$2,$8}' OFS="\t" Adipose_results_intron.txt  | awk -F ":|\t" 'FNR==NR{a[$1"_"$2"_"$3"_"$4"_"$9]=$12; next}{k=$1"_"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Adipose_new1_perind2.permutation.storey.txt - > Adipose_slope.ASE_intron.correlation_top_eQTL.txt

#############################################
#####Correlation between ASE_afc and cis-sQTL slope.
cd ~/ASE/Muscle

vcf_dir="~/Muscle_new/permutation"
################################################################
### at least 10 individuals with ASE data and a minimum of 8 reads per individual ( don't know how to filter; maybe ignore-----default parameter: --min-cov 8 in phaser_expr_matrix.py).
### aFC means var_het_afc? ( assume so)

###All the nominal-sig sQTLs
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05) print $1,$2,$8}' OFS="\t" Muscle_results.txt | awk 'FNR==NR{a[$1"\t"$2]=$5; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Muscle.nominals.sig.txt - > Muscle_slope.ASE.correlation.txt


####Only top significant sQTLs (FDR-storey < 0.05) per gene
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05 ) print $1,$2,$8}' OFS="\t" Muscle_results_intron.txt  | awk -F ":|\t" 'FNR==NR{a[$1"_"$2"_"$3"_"$4"_"$9]=$12; next}{k=$1"_"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Muscle_new1_perind2.permutation.storey.txt - > Muscle_slope.ASE_intron.correlation_top_sQTL.txt


#############################################
#####Correlation between ASE_afc and cis-eQTL slope.
cd /lustre/project/uvm_mckay/WGBS_others/RNA-seq/Revise_round_1/ASE/Liver

vcf_dir="/lustre/project/bull_scr/RNA-seq/sQTL/sqtl_detect/3_fastqtl/Liver_new/permutation"
################################################################
### at least 10 individuals with ASE data and a minimum of 8 reads per individual ( don't know how to filter; maybe ignore-----default parameter: --min-cov 8 in phaser_expr_matrix.py).
### aFC means var_het_afc? ( assume so)

###All the nominal-sig sQTLs
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05) print $1,$2,$8}' OFS="\t" Liver_results.txt | awk 'FNR==NR{a[$1"\t"$2]=$5; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Liver.nominals.2rd.sig.txt - > Liver_slope.ASE.correlation.txt


####Only top significant sQTLs (FDR-storey < 0.05) per gene
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05 ) print $1,$2,$8}' OFS="\t" Liver_results_intron.txt  | awk -F ":|\t" 'FNR==NR{a[$1"_"$2"_"$3"_"$4"_"$9]=$12; next}{k=$1"_"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Liver_new1_perind2.permutation.storey.txt - > Liver_slope.ASE_intron.correlation_top_sQTL.txt
