use warnings;
opendir DIR,"C:/Users/Nilesh/Desktop/excelincommand/";
#@files = grep(//,readdir(DIR));
@folders=grep(-d,readdir(DIR));
shift(@folders);
shift(@folders);

map{
$folder_name=$_;
opendir DIR1,"C:/Users/Nilesh/Desktop/excelincommand/$folder_name";
@files=grep(/\.txt/,readdir(DIR1));
open COMPILE,">C:/Users/Nilesh/Desktop/excelincommand/$folder_name/compile.txt";

map{
$file_name=$_;


##########ITDSEEK
	if($file_name=~m/itdseek/g)
	{
	open itdseek,"C:/Users/Nilesh/Desktop/excelincommand/$folder_name/$file_name";
	@arr5=<itdseek>;
		foreach $line5(@arr5)
		{
		print COMPILE "DP\tAD\tFREQ\t",$line5;			
		}	
	
	}

##########VARSCAN_SNP
	if($file_name=~m/varscan_snp/g)
	{
	open varscan_snp,"C:/Users/Nilesh/Desktop/excelincommand/$folder_name/$file_name";
	@arr1=<varscan_snp>;
		foreach $line1(@arr1)
		{
		@array1=split("\t",$line1);
		@info1=split("\:",$array1[173]);
		#print COMPILE $info5[3],"\t",$info5[5],"\t",$info5[6],"\t",$line5;
			if(($array1[5]=~ /exonic/g) && ($array1[8]=~ /nonsynonymous/g))
			{
				if($array1[12]<1 && $array1[13]<1 && $array1[19]<1)
				{
					print COMPILE $info1[3],"\t",$info1[5],"\t",$info1[6],"\t",$line1; 	
				}	
							
			}
		}	
	
	}


############VARSCAN_INDEL
	if($file_name=~m/varscan_IND/g)
	{
	open varscan_ind,"C:/Users/Nilesh/Desktop/excelincommand/$folder_name/$file_name";
		@arr2=<varscan_ind>;
			foreach $line2(@arr2)
			{
			@array2=split("\t",$line2);
			@info2=split("\:",$array2[173]);
				if($array2[5]=~ /exonic/)
				{
					if($array2[12]<1 && $array2[13]<1 && $array2[19]<1)
					{
						print COMPILE $info2[3],"\t",$info2[5],"\t",$info2[6],"\t",$line2; 	
					}			
								
				}
		}
	}
	
##########PINDEL_SI
	if($file_name=~m/Pindel_SI/g)
	{
	open pindel_si,"C:/Users/Nilesh/Desktop/excelincommand/$folder_name/$file_name";
	@arr3=<pindel_si>;
		foreach $line3(@arr3)
		{
		@array3=split("\t",$line3);
			if($array3[5]=~ /exonic/g)
			{
				if($array3[12]<1 && $array3[13]<1 && $array3[19]<1)
				{
					print COMPILE "\t\t\t",$line3; 	
				}	
							
			}
		}	
	
	}
	
##########PINDEL_D
	if($file_name=~m/Pindel_D/g)
	{
	open pindel_d,"C:/Users/Nilesh/Desktop/excelincommand/$folder_name/$file_name";
	@arr4=<pindel_d>;
		foreach $line4(@arr4)
		{
		@array4=split("\t",$line4);
			if($array4[5]=~ /exonic/g)
			{
				if($array4[12]<1 && $array4[13]<1 && $array4[19]<1)
				{
					print COMPILE "\t\t\t",$line4; 	
				}	
							
			}
		}	
	
	}	


	
}@files;
}@folders;


map{
my %seen;
open COMPILE,"C:/Users/Nilesh/Desktop/excelincommand/$folder_name/compile.txt";
while ( <COMPILE> ) { 
    print if $seen{$_}++;
}
}@folder;





	
	
	

