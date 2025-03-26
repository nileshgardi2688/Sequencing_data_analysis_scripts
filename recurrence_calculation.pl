use strict;
opendir DIR,"/media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/";
open F3,">/media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/recurrent_data_with_sample_names.out";
#mkdir("/media/LUN1/NIlesh_analysis/Exome/MuTECT/filtered");
my @files = grep(/\.txt/,readdir(DIR));
my %hash1=();
map{
my $code=$_;
open F1,"/media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/$code";
my @arr=<F1>; close(F1);
foreach my $line(@arr){
        chomp($line);
        my @array=split("\t",$line);
        push(@{$hash1{$array[0]}},$array[1]);
          }
}@files;

foreach(keys(%hash1)){
                      local $"=",";
      print F3 "$_\t@{$hash1{$_}}\n";
}