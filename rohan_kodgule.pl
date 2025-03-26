use warnings;
###########Open the Parental folder in opendir command. It will contain all the subfolders name with sample id. In these subfolder you place your R1 and R2 fastq files. Name of the fastq files would start with your sample ID followed by underscore R1.fastq or undersc R2.fastq. Pleaes make it generic, so trhat program run smoothely.
opendir DIR, "/home/darthvader/NGS_Analysis/";
@only_directories = grep(!/\./,readdir(DIR));


map{
$sample_id=$1;

$echo1="echo "Sai generation started for R1 sample $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command1="bwa aln -t 12 -f /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R1.sai" /home/darthvader/Downloads/references/hg19_all.fasta /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R1.fastq.gz""
system($echo1);
system($command1);

$echo2="echo "Sai generation started for R2 sample $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command2= "bwa aln -t 12 -f /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R2.sai" /home/darthvader/Downloads/references/hg19_all.fasta /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R2.fastq.gz""
system($echo2);
system($command2); 

$echo3="echo "Sam generation started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command3= "bwa sampe -f /home/darthvader/NGS_Analysis/$sample_id/$sample_id".sam" -r "@RG\tID:TSCC-12T\tPL:ILLUMINA\tLB:LIB-EXO\tSM:UNKNOWN\tPI:200" /home/darthvader/Downloads/references/hg19_all.fasta /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R1.sai" /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R2.sai" /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R1.fastq.gz" /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R2.fastq.gz""
system($echo3);
system($command3);

$echo4="echo "fxd mate generation started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command4= "java -Xmx20G -jar /home/darthvader/Downloads/Picard/picard.jar FixMateInformation I= /home/darthvader/NGS_Analysis/$sample_id/$sample_id".sam" o= /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_fxd.sam" VALIDATION_STRINGENCY=SILENT TMP_DIR=tmp/" 
system($echo4);
system($command4);

$echo5="echo "bam generation started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command5= "/home/darthvader/Downloads/samtools-1.3.1/samtools view -bT /home/darthvader/Downloads/references/hg19_all.fa /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_fxd.sam" > /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_fxd.bam""
system($echo5);
system($command5);

$echo6="echo "Sorting bam started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command6= "/home/darthvader/Downloads/samtools-1.3.1/samtools sort /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_fxd.bam" > /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_fxd_sorted.bam"" 
system($echo6);
system($command6);

$echo7="echo "Indexing bam started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command7= "/home/darthvader/Downloads/samtools-1.3.1/samtools index /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_fxd_sorted.bam" > /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_fxd_sorted_index.bai"" 
system($echo7);
system($command7);

$echo8="echo "mpileup started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command8= "/home/darthvader/Downloads/samtools-1.3.1/samtools mpileup -f /home/darthvader/Downloads/references/hg19_all.fa /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_fxd_sorted.bam" > /home/darthvader/NGS_Analysis/$sample_id/$sample_id".mpileup""
system($echo8);
system($command8);

$echo9="echo "mpileup2snp started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command9= "java -jar /home/darthvader/Downloads/VarScan.v2.3.9.jar mpileup2snp /home/darthvader/NGS_Analysis/$sample_id/$sample_id".mpileup" --min- coverage 8 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.03 --p-value 1e-4 --output-vcf 1 > /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_varscan_snp.vcf"" 
system($echo9);
system($command9);

$echo10="echo "mpileup2indel started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command10= "java -jar /home/darthvader/Downloads/VarScan.v2.3.9.jar mpileup2indel /home/darthvader/NGS_Analysis/$sample_id/$sample_id".mpileup" --min- coverage 8 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.03 --p-value 1e-4 --output-vcf 1 > /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_varscan_IND.vcf"" 
system($echo10);
system($command10);

$echo11="echo "AnnovarSNP started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command11= "perl /home/darthvader/Downloads/annovar/table_annovar.pl /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_varscan_snp.vcf" /home/darthvader/Downloads/annovar/humandb/ -buildver hg19 -out /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_snp" -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp147,cosmic70,exac03,kaviar,icgc21,clinvar_20131105,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -otherinfo"
system($echo11);
system($command11);

$echo12="echo "AnnovarINDEL started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command12=  "perl /home/darthvader/Downloads/annovar/table_annovar.pl /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_varscan_IND.vcf" /home/darthvader/Downloads/annovar/humandb/ -buildver hg19 -out /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_ind" -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp147,cosmic70,exac03,kaviar,icgc21,clinvar_20131105,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -otherinfo"
system($echo12);
system($command12); 


$echo13="echo "PINDEL started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command13= "/home/darthvader/Downloads/samtools-1.3.1/./samtools view /home/darthvader/NGS_Analysis/$sample_id/$sample_id".sam" | /home/darthvader/pindel/./sam2pindel - /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_output.pindel" 500 tumour 0 Illumina-PairEnd"
system($echo13);
system($command13); 

$echo14="echo "PINDEL_analysis started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command14= "/home/darthvader/pindel/./pindel -f /home/darthvader/Downloads/references/hg19_all.fa  -p /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_output.pindel" -c ALL -M 30 -a 5 -u 0.3 -o /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel""
system($echo14);
system($command14);

$echo15="echo "PINDEL2vcf started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command15= "/home/darthvader/pindel/./pindel2vcf -p /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel_SI" -r /home/darthvader/Downloads/references/hg19_all.fasta -R hg19 -d 20111123 -v /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel_SI.vcf"" 
system($echo15);
system($command15);


$echo16="echo "PINDEL2vcf started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command16= "/home/darthvader/pindel/./pindel2vcf -p /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel_D" -r /home/darthvader/Downloads/references/hg19_all.fasta -R hg19 -d 20111123 -v /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel_D.vcf""
system($echo16);
system($command16);

$echo17="echo "annovar started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command17= "perl /home/darthvader/Downloads/annovar/table_annovar.pl /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel_SI.vcf" /home/darthvader/Downloads/annovar/humandb/ -buildver hg19 -out /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel_SI" -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp147,cosmic70,exac03,kaviar,icgc21,clinvar_20131105,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -otherinfo "
system($echo17);
system($command17);

$echo18="echo "Annovar started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command18= "perl /home/darthvader/Downloads/annovar/table_annovar.pl /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel_D.vcf" /home/darthvader/Downloads/annovar/humandb/ -buildver hg19 -out /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_Pindel_D" -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp147,cosmic70,exac03,kaviar,icgc21,clinvar_20131105,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -otherinfo"
system($echo18);
system($command18);


$echo19="echo "BWA mem started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command19= "bwa mem -R '@RG\tID:TSCC\tSM:AMLsample' -M -t 12 /home/darthvader/Downloads/references/hg19_all.fasta /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R1.fastq.gz" /home/darthvader/NGS_Analysis/$sample_id/$sample_id"_R2.fastq.gz" | /home/darthvader/Downloads/samtools-1.3.1/samtools view -bS - | /home/darthvader/Downloads/samtools-1.3.1/samtools sort > /home/darthvader/NGS_Analysis/$sample_id/$sample_id".bam""
system($echo19);
system($command19);

$echo20="echo "samtoold index started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command20= "/home/darthvader/Downloads/samtools-1.3.1/samtools index /home/darthvader/NGS_Analysis/$sample_id/$sample_id".bam" > /home/darthvader/NGS_Analysis/$sample_id/$sample_id".bam.bai"" 
system($echo20);
system($command20);

$echo21="echo "itdseek started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command21= "bash /home/darthvader/Downloads/itdseek-master/itdseek.sh /home/darthvader/NGS_Analysis/$sample_id/$sample_id".bam" /home/darthvader/Downloads/references/hg19_all.fasta /home/darthvader/Downloads/samtools-1.3.1/samtools > /home/darthvader/NGS_Analysis/$sample_id/$sample_id".itdseek.vcf"" 
system($echo21);
system($command21);

$echo22="echo "annovar_last started for $sample_id at $(date)">>/home/darthvader/NGS_Analysis/$sample_id/$sample_id"_variantcall_echo.txt""
$command22= "perl /home/darthvader/Downloads/annovar/table_annovar.pl /home/darthvader/NGS_Analysis/$sample_id/$sample_id".itdseek.vcf" /home/darthvader/Downloads/annovar/humandb/ -buildver hg19 -out /home/darthvader/NGS_Analysis/ITDSEEK/$sample_id"_itdseek" -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp147,cosmic70,exac03,kaviar,icgc21,clinvar_20131105,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -otherinfo"
system($echo22);
system($command22);


}@only_directories;
