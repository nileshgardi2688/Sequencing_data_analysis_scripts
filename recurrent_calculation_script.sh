mkdir /media/LUN1/NIlesh_analysis/Exome/Recurrence_analysis_exome/recal_mutect_ori/

#################Exome data#########################################
for f in $(find /media/LUN1/NIlesh_analysis/Exome/GATK_SNP_INDEL_MUTECT_merge/recal_mutect_ori/* -name "*.txt" -type f ); 
do 
    sample=$(basename $f);
	sample_name="${sample%.txt}";
	#vcf_filename_without_extn="${vcf_filename_with_extn%.*}" ;
	sed 's/\t/_/g' $f >> /media/LUN1/NIlesh_analysis/Exome/Recurrence_analysis_exome/recal_mutect_ori/$sample_name"exo.txt"
done;


############Remove header line frm each file#############################################
sed -i '/chr_/d' /media/LUN1/NIlesh_analysis/Exome/Recurrence_analysis_exome/recal_mutect_ori/*.*
sed -i 's/chr//g' /media/LUN1/NIlesh_analysis/Exome/Recurrence_analysis_exome/recal_mutect_ori/*.*

############################Add sample name to end of each line in each file
cd /media/LUN1/NIlesh_analysis/Exome/Recurrence_analysis_exome/recal_mutect_ori/
for f in *.*; do sed -i "s/$/\t$f/" $f; done
sed -i 's/\.txt//g' /media/LUN1/NIlesh_analysis/Exome/Recurrence_analysis_exome/recal_mutect_ori/*.txt 

######################Run perl program for recurrence analysis
perl /media/LUN1/NIlesh_analysis/Exome/recurrence_calculation.pl

##################Annotate with COSMIC, DBSNP, and mylab######################
perl /media/LUN1/NIlesh_analysis/Exome/Map_cos_snp_myLab.pl


####################Merge all files into single#######################
#cat /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/*.txt | sort | uniq >> /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/Recurrence_output.out

########################Recurrence calculation###########################
#while read line;  
#do recurrent_pair=$(grep $line /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/*.txt| cut -f 6 -d '/'| cut -f 1 -d '.'| tr '\n' '.'| sed 's/.$//g'); echo -e $line'\t'$recurrent_pair >> /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/output_with_samples_exome.txt ;
#done < "/media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/Recurrence_output.out"

#####################Making oncotator input file#######################################################
#sed -i 's/_/\t/g' /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/output_with_samples_exome.txt
#sed -i '1ichr\tstart\tend\tref_allele\talt_allele\tsource' /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/output_with_samples_exome.txt

################################Oncotator#################################
#oncotator -v --input_format MAFLITE --db-dir /home/adlab/ngs_data/databases/oncotator_ds_tmp --output_format=TCGAMAF /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/Recurrent_Map_cos_snp_mylab.tab /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/output_with_samples.maf hg19