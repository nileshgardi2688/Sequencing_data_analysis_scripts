use warnings;
open FH1, "</media/LUN1/NIlesh_analysis/Exome/a.txt";
open FH2, ">>/media/LUN1/NIlesh_analysis/Exome/coverage.info";

@arr=<FH1>;

foreach $line(@arr)
{
chop($line);
 @content=split('\t',$line);
 @data=split('\,',$content[0]);
 @samples=split('\,',$content[1]);
 
 if($#samples > 0)
 {
 
 for($i=0;$i<=$#samples;$i++)
 {
 open file1,"/media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/merge_same_tumor/".$samples[$i]."_combine.vcf";
 
      foreach $line1(<file1>)
      {
       @array=split('\t',$line1);
 
                           if($data[0] eq $array[0] && $data[1] eq $array[1] && $data[2] eq $array[2] && $data[3] eq $array[3])
                           {
                           print FH2 $samples[$i]."\t".$line1;
                           }
      }
 }
 }
 
 else
 {
 open file2,"/media/LUN1/NIlesh_analysis/Exome/accessing_coverage_info_of_last_mutations/merge_same_tumor/".$content[1]."_combine.vcf";
      
      foreach $line2(<file2>)
      {
       @array1=split('\t',$line2);
 
                           if($data[0] eq $array1[0] && $data[1] eq $array1[1] && $data[2] eq $array1[2] && $data[3] eq $array1[3])
                           {
                           print FH2 $content[1]."\t".$line2;
                           }
      }
 
 }
 


}