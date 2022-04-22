cd ~/ASE/Adipose

vcf_dir="~/vcf_tissues/Adipose"
################################################################
### at least 10 individuals with ASE data and a minimum of 8 reads per individual ( don't know how to filter; maybe ignore-----default parameter: --min-cov 8 in phaser_expr_matrix.py).
### aFC means var_het_afc? ( assume so)

###All the nominal-sig eQTLs
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05) print $1,$2,$8}' OFS="\t" Adipose_results.txt | awk 'FNR==NR{a[$1"\t"$2]=$5; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Adipose.nominals.sig.txt - > Adipose_slope.ASE.correlation.txt


####Only top significant eQTLs (FDR-storey < 0.05) per gene
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05 ) print $1,$2,$8}' OFS="\t" Adipose_results.txt  | awk 'FNR==NR{a[$1"\t"$6]=$9; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Adipose.permutations.storey.txt - > Adipose_slope.ASE.correlation_top_eQTL.txt

#############################################
#####Correlation between ASE_afc and cis-eQTL slope.


cd ~/ASE/Muscle

vcf_dir="/lustre/project/bull_scr/RNA-seq/SNP_calling/Final_vcf_file/filter_MAF_0.05_DR2_0.8/eQTL/vcf_tissues/Muscle"
################################################################
### at least 10 individuals with ASE data and a minimum of 8 reads per individual ( don't know how to filter; maybe ignore-----default parameter: --min-cov 8 in phaser_expr_matrix.py).
### aFC means var_het_afc? ( assume so)

###All the nominal-sig eQTLs
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05) print $1,$2,$8}' OFS="\t" Muscle_results.txt | awk 'FNR==NR{a[$1"\t"$2]=$5; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Muscle.nominals.sig.txt - > Muscle_slope.ASE.correlation.txt


####Only top significant eQTLs (FDR-storey < 0.05) per gene
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05 ) print $1,$2,$8}' OFS="\t" Muscle_results.txt  | awk 'FNR==NR{a[$1"\t"$6]=$9; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Muscle.permutations.storey.txt - > Muscle_slope.ASE.correlation_top_eQTL.txt



#############################################
#####Correlation between ASE_afc and cis-eQTL slope.


cd ~/ASE/Liver

vcf_dir="~/vcf_tissues/Liver"
################################################################
### at least 10 individuals with ASE data and a minimum of 8 reads per individual ( don't know how to filter; maybe ignore-----default parameter: --min-cov 8 in phaser_expr_matrix.py).
### aFC means var_het_afc? ( assume so)

###All the nominal-sig eQTLs
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05) print $1,$2,$8}' OFS="\t" Liver_results.txt | awk 'FNR==NR{a[$1"\t"$2]=$5; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Liver.nominals.sig.txt - > Liver_slope.ASE.correlation.txt


####Only top significant eQTLs (FDR-storey < 0.05) per gene
awk -F "\t" '{if ($5>=10 &&  $11 < 0.05 ) print $1,$2,$8}' OFS="\t" Liver_results.txt  | awk 'FNR==NR{a[$1"\t"$6]=$9; next}{k=$1"\t"$2; if ( k in a) print $0"\t"a[k]}' ${vcf_dir}/Liver.permutations.storey.txt - > Liver_slope.ASE.correlation_top_eQTL.txt



