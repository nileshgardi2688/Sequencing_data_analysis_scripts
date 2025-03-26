##########################5X filter for SNP###################################################################
#!/bin/bash
#mkdir /media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP/T-N_snp/coverage_5x_filter_snp
for file in $(find /media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP/T-N_snp/* -type f -name "2_T-N.vcf") ; 
do
vcf_filename_with_extn=$(basename $file) ;
#vcf_filename=$(basename "$1")
vcf_base_name="${vcf_filename_with_extn%_T-N*}"
alt_base_depth=5

	cat $file | while read LINE 
	do 
		if echo $LINE | grep -Eq '^#'
		then
				continue;
        #echo $LINE >> $vcf_base_name"_"$alt_base_depth"x_filtered.vcf"
		else
				
				echo $LINE | awk -v alt_base_count="$alt_base_depth" -v line="$LINE" '{split($10,a,":|,");if(a[3]>=alt_base_count) print line;}' >> /media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP/T-N_snp/coverage_5x_filter_snp/$vcf_base_name"_T-N_"$alt_base_depth"x_filtered.vcf"

		fi
	done
	echo "Completed!!! Output File generated as $vcf_base_name"_T-N_"$alt_base_depth"x_filtered.vcf""

done;


##########################5X filter for Indel###################################################################
#!/bin/bash
#mkdir /media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel/T-N_indels/coverage_5x_filter_Indel
for file in $(find /media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel/T-N_indels/* -type f -name "2_T-N.vcf") ; 
do
vcf_filename_with_extn=$(basename $file) ;
vcf_filename=$(basename "$1")
vcf_base_name="${vcf_filename_with_extn%.*}"
alt_base_depth=5

  cat $file | while read LINE 
	do 
		if echo $LINE | grep -Eq '^#'
		then
				continue;
        #echo $LINE >> $vcf_base_name"_"$alt_base_depth"x_filtered.vcf"
		else
				
				echo $LINE | awk -v alt_base_count="$alt_base_depth" -v line="$LINE" '{split($10,a,":|,");if(a[3]>=alt_base_count) print line;}' >> /media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel/T-N_indels/coverage_5x_filter_Indel/$vcf_base_name"_"$alt_base_depth"x_filtered.vcf"

	fi
	done
	echo "Completed!!! Output File generated as $vcf_base_name"_"$alt_base_depth"x_filtered.vcf""

done;

echo "Analysis of filtering completed!!!!!"


#################Extracting 5 fields from filtered file (chr_start_end_ref_alt) for SNP######################################

#mkdir /media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP/T-N_snp/coverage_5x_filter_snp/chr_start_end_ref_alt_snp/
for file in $(find /media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP/T-N_snp/coverage_5x_filter_snp/* -type f -name "2_T-N_5x_filtered.vcf") ; 
do
vcf_filename_with_extn=$(basename $file) ;
vcf_base_name="${vcf_filename_with_extn%.*}"
awk 'BEGIN { print "chr","\t","start","\t","end","\t","ref_allele","\t","alt_allele"} { print $1,"\t",$2,"\t",$3=length($4)+$2-1,"\t",$4,"\t",$5 }' $file | sed 's/ //g' >> /media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP/T-N_snp/coverage_5x_filter_snp/chr_start_end_ref_alt_snp/$vcf_base_name".vcf" 
done;

#################Extracting 5 fields from filtered file (chr_start_end_ref_alt) for INDEL#########################################
#mkdir /media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel/T-N_indels/coverage_5x_filter_Indel/chr_start_end_ref_alt_indel/
for file in $(find /media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel/T-N_indels/coverage_5x_filter_Indel/* -type f -name "2_T-N_5x_filtered.vcf") ; 
do
vcf_filename_with_extn=$(basename $file) ;
vcf_base_name="${vcf_filename_with_extn%.*}"
awk 'BEGIN { print "chr","\t","start","\t","end","\t","ref_allele","\t","alt_allele"} { print $1,"\t",$2,"\t",$3=length($4)+$2-1,"\t",$4,"\t",$5 }' $file | sed 's/ //g' >> /media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel/T-N_indels/coverage_5x_filter_Indel/chr_start_end_ref_alt_indel/$vcf_base_name".vcf" 
done;