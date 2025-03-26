mkdir /media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/
for file in $(find /media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP/T-N_snp/* -type f -name "*_T-N.vcf") ; 
do
  vcf_filename_with_extn=$(basename $file) ;
	vcf_filename_without_extn="${vcf_filename_with_extn%_T-N.vcf}" ;
       
    sed -e '1,2d' < $file | cut -f 1,2,4,5,10 > /media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/$vcf_filename_without_extn"_GATK_SNP.vcf";

done 

for file in $(find /media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel/T-N_indels/* -type f -name "*_T-N.vcf") ; 
do
  vcf_filename_with_extn=$(basename $file) ;
	vcf_filename_without_extn="${vcf_filename_with_extn%_T-N.vcf}" ;
       
    sed -e '1,1d' < $file | cut -f 1,2,4,5,10 > /media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/$vcf_filename_without_extn"_GATK_INDEL.vcf";

done 

for file in $(find /media/LUN1/NIlesh_analysis/Exome/MuTECT/* -type f -name "*T_fxd_sorted_DupRm_realn_recal.mutect.txt") ; 
do
  vcf_filename_with_extn=$(basename $file) ;
	vcf_filename_without_extn="${vcf_filename_with_extn%T_fxd_sorted_DupRm_realn_recal.mutect.txt}" ;
       
    sed -e '1,2d' < $file | cut -f 1,2,3,4,20,21 > /media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/$vcf_filename_without_extn"_mutect.vcf";

done 

mkdir /media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/merge_same_tumor/
for file in $(find /media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/* -type f -name "*.vcf") ; 
do
  vcf_filename_with_extn=$(basename $file) ;
	vcf_filename_without_extn="${vcf_filename_with_extn%_*}" ;
       
    cat /media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/$vcf_filename_without_extn* > /media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/merge_same_tumor/$vcf_filename_without_extn"_combine.vcf";

done 






