#use warnings;
opendir DIR,"/media/LUN1/NIlesh_analysis/Exome/MuTECT/";
mkdir("/media/LUN1/NIlesh_analysis/Exome/MuTECT/filtered");
@files = grep(/\.txt/,readdir(DIR));
#@files=grep(-f,readdir(DIR));
#shift(@files);
#shift(@files);

map{
$code=$_;
$code=~s/\.txt//g;
open F1,"/media/LUN1/NIlesh_analysis/Exome/MuTECT/$code.txt";
open F3,">/media/LUN1/NIlesh_analysis/Exome/MuTECT/filtered/$code"."_5_fields".".txt";
#undef $/;
@arr1=<F1>;
	
	if($arr1[0] eq "## muTector v1.0.27200\n")
	{
		foreach $line(@arr1)
		{
			#shift(@arr1);
			#shift(@arr1);
			@array=split("\t",$line);
			
      			if($array[28]==0 && (($array[19]+$array[20])>=5))
			    {
				        if($array[19]==0)
				          {
				            print F3 $array[0]."\t".$array[1]."\t".$array[1]."\t".$array[2]."\t".$array[3],"\n";
				          }

				        else
				          {
					              if(($array[20]/$array[19])>=0.2)
					              {
					              #print F3 $array[3],"\n";
					              print F3 $array[0]."\t".$array[1]."\t".$array[1]."\t".$array[2]."\t".$array[3],"\n";
					              }
				          }
			    }
					
					#print F3 $array[3],"\n";
			#		print F3 $array[0]."\t".$array[1]."\t".$array[1]."\t".$array[2]."\t".$array[3],"\n";
					
				
		 
		}
	
}
}@files;


=head

opendir DIR1,"/media/LUN1/NIlesh_analysis/Exome/MuTECT/filtered";
mkdir("/media/LUN1/NIlesh_analysis/Exome/MuTECT/filtered/oncotator_input");
$copy_command="cp oncotator.sh /home/adlab/Nilesh_data/Journey_from_MuTECT_to_oncotator/MuTECT/filtered/oncotator_input/";
system($copy_command);
mkdir("MuTECT/filtered/oncotator_input/oncotator_output");
@vcf_files = grep(//,readdir(DIR1));
shift(@vcf_files);
shift(@vcf_files);

map{
$vcf_code=$_;
$vcf_code=~s/\.txt//g;
#$vcf_code=~s/\.MuTECT.*//g;
open F4,"MuTECT/filtered/$vcf_code.txt";
open F5,"12T-N-N_Oncotator_input.txt";
$vcf_code=~s/\.MuTECT.*//g;
open F6,">MuTECT/filtered/oncotator_input/$vcf_code"."_oncotator".".txt";
@arr2=<F4>;
@arr3=<F5>;
$unknown="Unknown\t" x 34;
#$unknown=~s/\n//g;
#chomp($unknown);

print F6 $arr3[0];
	foreach $line1(@arr2)
	{
	chomp($line1);
	print F6 "hg19\t".$line1."\t".$unknown,"\n";
	}

}@vcf_files;


system("sh /home/adlab/Nilesh_data/Journey_from_MuTECT_to_oncotator/MuTECT/filtered/oncotator_input/oncotator.sh");





