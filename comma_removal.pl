use warnings;
use strict;
opendir DIR1,"/media/LUN1/NIlesh_analysis/Exome/GATK/SNP/";
mkdir "/media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP";
my @files=grep{s/\.vcf$//}readdir(DIR1);



map{
my $name=$_;
open FH1,"/media/LUN1/NIlesh_analysis/Exome/GATK/SNP/$name.vcf";
open FH2,">/media/LUN1/NIlesh_analysis/Exome/GATK/SNP/corrected_SNP/$name.vcf";
my @array=<FH1>; 
foreach my $line(@array){
	chomp($line);
	if($line=~/\#/g){
		print FH2 $line,"\n";
	}
	else
	{
		my @arr=split("\t",$line);
		if($arr[4]=~/,/g){
			my @alt=split(",",$arr[4]);
			for(my $i=0;$i<$#alt+1;$i++){
				my @cov=split(":",$arr[9]);
				my @alt_depth=split(",",$cov[1]);
				my $j=$i+1;
				my $str=$cov[0]."\:".$alt_depth[0].",".$alt_depth[$j]."\:".$cov[2]."\:".$cov[3]."\:".$cov[4];
				#print "$alt[$i]\t$str\n";
				print FH2 $arr[0]."\t".$arr[1]."\t".$arr[2]."\t".$arr[3]."\t".$alt[$i]."\t".$arr[5]."\t".$arr[6]."\t".$arr[7]."\t".$arr[8]."\t".$str,"\n";
			}
		}else{
			print FH2 $line,"\n";
		}
	}
	}
	
}@files;



opendir DIR2,"/media/LUN1/NIlesh_analysis/Exome/GATK/Indel/";
mkdir "/media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel";
my @files1=grep{s/\.vcf$//}readdir(DIR2);

map{
my $name1=$_;
open FH3,"/media/LUN1/NIlesh_analysis/Exome/GATK/Indel/$name1.vcf";
open FH4,">/media/LUN1/NIlesh_analysis/Exome/GATK/Indel/corrected_Indel/$name1.vcf";
my @arr1ay1=<FH3>; 
foreach my $line1(@arr1ay1){
	chomp($line1);
	if($line1=~/\#/g){
		print FH4 $line1,"\n";
	}
	else
	{
		my @arr1=split("\t",$line1);
		if($arr1[4]=~/,/g){
			my @alt1=split(",",$arr1[4]);
			for(my $k=0;$k<$#alt1+1;$k++){
				my @cov1=split(":",$arr1[9]);
				my @alt_depth1=split(",",$cov1[1]);
				my $m=$k+1;
				my $str1=$cov1[0]."\:".$alt_depth1[0].",".$alt_depth1[$m]."\:".$cov1[2]."\:".$cov1[3]."\:".$cov1[4];
				#print "$alt[$i]\t$str\n";
				print FH4 $arr1[0]."\t".$arr1[1]."\t".$arr1[2]."\t".$arr1[3]."\t".$alt1[$k]."\t".$arr1[5]."\t".$arr1[6]."\t".$arr1[7]."\t".$arr1[8]."\t".$str1,"\n";
			}
		}else{
			print FH4 $line1,"\n";
		}
	}
	}
	
}@files1;
