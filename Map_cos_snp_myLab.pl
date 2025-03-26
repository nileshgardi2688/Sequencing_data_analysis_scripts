use strict;
use warnings;
use Data::Dumper;
#my $inDir="chunks/";
#open(FW,">>"."OncotatorOutput.txt")||die "Can not open FW.\n";
open(F1,"/media/LUN1/NIlesh_analysis/COSMIC_DBSNP/Cosmic_key_20150318.txt")||die "Can not open F1.\n";
open(F2,"/media/LUN1/NIlesh_analysis/COSMIC_DBSNP/DBSNP_key_20150318.txt")||die "Can not open F2.\n";
open(F3,"/media/LUN1/NIlesh_analysis/COSMIC_DBSNP/mylab_key_20150318.txt")||die "Can not open F3.\n";
open(F4,"/media/LUN1/NIlesh_analysis/Exome/Recurrence_analysis_exome/recal_mutect_ori/recurrent_data_with_sample_names.out")||die "Can not open F4.\n";
open(FW,">"."/media/LUN1/NIlesh_analysis/Exome/Recurrence_analysis_exome/recal_mutect_ori/recurrent_Map_cos_snp_mylab.txt")||die "Can not open FW.\n";
my (%cosmic, %snp, %mylab)=();
%cosmic=&file(<F1>);
%snp=&file(<F2>);
%mylab=&file(<F3>);
#print Dumper(%cosmic);
print FW "chr\tstart\tend\tref_allele\talt_allele\tSource\tCosmic\tdbSNP\tMyLabDB\n";

while(<F4>){
chomp($_);
            my @file=split("\t",$_);
            my @key=split("_",$file[0]);
            my $k=$key[0]."_".$key[1]."_".$key[3]."_".$key[4];
            print FW "chr$key[0]\t$key[1]\t$key[2]\t$key[3]\t$key[4]\t$file[1]\t";
            local $"=",";
            if(defined $cosmic{$k}){print FW "@{$cosmic{$k}}"."\t"}else{print FW "-\t"}
            if(defined $snp{$k}){print FW "@{$snp{$k}}\t"}else{print FW "-\t"}
            if(defined $mylab{$k}){print FW "1\n"}else{print FW "-\n"}
}

####subroutine : To make hash of file
sub file{
    my @arr=@_;
    my %file1=();
    foreach(@arr){
                  chomp($_);
                  my @file2=split("\t",$_);
                  #$file1{$file2[0]}=$file2[3];
                push(@{$file1{$file2[0]}},$file2[1]); 
    
    }
    return %file1;
}

