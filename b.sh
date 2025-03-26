#sed -i 's/ /_/g' "/media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/corrected_Oncotator_output_Recurrent_Map_cos_snp_mylab.maf"

# Extracting individual tumor data from recurrent file
mkdir /media/LUN1/NIlesh_analysis/Individual_tumor/
sample_list=($(cat /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/corrected_Oncotator_output_Recurrent_Map_cos_snp_mylab.maf | cut -f 170 | egrep -v '\,' |sort |uniq))
for i in ${sample_list[0]}
do
 while IFS=$'\t' read -r -a element
 do 
    chr=${element[169]}
    #echo $i"\t"$chr"\n";
    
    if [ "$chr" == "$i" ] || [[ "$chr" =~ .*\,"$i"\,.* ]] || [[ "$chr" =~ ^"$i"\,.* ]] || [[ "$chr" =~ .*\,"$i"$ ]]; then
    #if [[ "$chr" =~ .*\."$i".* ]]; then
    #echo ${element[@]} >> xyz/$i".txt"
    ( IFS=$'\t'; echo ${element[@]} >> /media/LUN1/NIlesh_analysis/Individual_tumor/$i".txt")
    fi
 done < /media/LUN1/NIlesh_analysis/Recurrence_analysis_exome_transcriptome/corrected_Oncotator_output_Recurrent_Map_cos_snp_mylab.maf
done

sed -i 's/ /\t/g' /media/LUN1/NIlesh_analysis/Individual_tumor/*.txt

# Counting cosmic, dbsnp, cosmic+dbsnp, novel entries for individual tumor

for f in $(find /media/LUN1/NIlesh_analysis/Individual_tumor/* -name "*.txt" -type f ); 
do
  sample=$(basename $f);
	sample_name="${sample%.*}"
awk -v file="$sample_name" -F'\t' 'BEGIN{novel=0; cosmic=0; dbsnp=0; mylab=0; cosdb=0; cosmy=0; dbmy=0; cosdbmy=0} 
    {if($121 ~ /-/ && $257 ~ /-/ && $168 ~ /-/) novel++;
    else if($121 !~ /-/ && $257 ~ /-/ && $168 ~ /-/) cosmic++; 
    else if($121 ~ /-/ && $257 !~ /-/ && $168 ~ /-/) dbsnp++;
    else if($121 ~ /-/ && $257 ~ /-/ && $168 !~ /-/) mylab++;
    else if($121 !~ /-/ && $257 !~ /-/ && $168 ~ /-/) cosdb++;
    else if($121 !~ /-/ && $257 ~ /-/ && $168 !~ /-/) cosmy++;
    else if($121 ~ /-/ && $257 !~ /-/ && $168 !~ /-/) dbmy++;
    else if($121 !~ /-/ && $257 !~ /-/ && $168 !~ /-/) cosdbmy++;}
    
    END{print "Sample Name=""\t"file
    print "Total number of variants=""\t"novel+cosmic+dbsnp+mylab+cosdb+cosmy+dbmy+cosdbmy;
    print "Novel entries=""\t"novel;
    print "Exclusive Cosmic entries=""\t"cosmic;
    print "Exclusive DBsnp entries=""\t"dbsnp;
    print "Exclusive MyLAB entries=""\t"mylab;
    print "Cosmic+DBsnp common entries=""\t"cosdb;
    print "Cosmic+MyLAB common entries=""\t"cosmy;
    print "DBsnp+MyLAB common entries=""\t"dbmy;
    print "Cosmic+DBsnp+MyLAB common entries=""\t"cosdbmy;
    print "Total Cosmic entries=""\t"cosmic+cosdb+cosmy+cosdbmy"\n";
    #print "Total DBsnp entries=""\t"dbsnp+cosdb"\n";
    }' < $f >> /media/LUN1/NIlesh_analysis/Individual_tumor/COSMIC_DBSNP_INFO.info
done



###################################### all types of mutations calculations
for f in $(find /media/LUN1/NIlesh_analysis/Individual_tumor/* -name "*.txt" -type f ); 
do
  sample=$(basename $f);
	sample_name="${sample%.*}"
awk -v file="$sample_name" -F'\t' 'BEGIN{missense=0; nonsense=0; nonstop=0; splice=0; threeutr=0; fiveflank=0; fiveutr=0; igr=0; intron=0; silent=0; DenovoStartInFrame=0; DenovoStartOutOfFrame=0; FrameShiftDel=0; FrameShiftIns=0; InFrameDel=0; InFrameIns=0; lincRNA=0; RNA=0; StartCodonDel=0; StartCodonIns=0; StartCodonSNP=0; StopCodonDel=0; StopCodonIns=0;} 
    {if($9 ~ "3" && $9 ~ "UTR") threeutr++;
    else if($9 ~ "5" && $9 ~ "Flank") fiveflank++;
    else if($9 ~ "5" && $9 ~ "UTR") fiveutr++;
    else if($9 == "IGR") igr++;
    else if($9 == "Intron") intron++;
    else if($9 == "De_novo_Start_InFrame") DenovoStartInFrame++;
    else if($9 == "De_novo_Start_OutOfFrame") DenovoStartOutOfFrame++;
    #else if($9 == "Non-coding_Transcript") noncoading++;
    else if($9 == "Frame_Shift_Del") FrameShiftDel++;
    else if($9 == "Frame_Shift_Ins") FrameShiftIns++;
    else if($9 == "In_Frame_Del") InFrameDel++;
    else if($9 == "In_Frame_Ins") InFrameIns++;
    else if($9 == "lincRNA") lincRNA++;
    else if($9 == "RNA") RNA++;
    else if($9 == "Start_Codon_Del") StartCodonDel++;
    else if($9 == "Start_Codon_Ins") StartCodonIns++;
    else if($9 == "Start_Codon_SNP") StartCodonSNP++;
    else if($9 == "Stop_Codon_Del") StopCodonDel++;
    else if($9 == "Stop_Codon_Ins") StopCodonIns++;
    else if($9 == "Missense_Mutation") missense++;
    else if($9 == "Nonsense_Mutation") nonsense++;
    else if($9 == "Nonstop_Mutation") nonstop++;
    else if($9 == "Silent") silent++;
    else if($9 == "Splice_Site") splice++;}
    END{print "Sample Name=""\t"file
    print "Total number of variants=""\t"missense+nonsense+nonstop+splice+threeutr+fiveflank+fiveutr+igr+intron+noncoading+silent+DenovoStartInFrame+DenovoStartOutOfFrame+FrameShiftDel+FrameShiftIns+InFrameDel+InFrameIns+lincRNA+RNA+StartCodonDel+StartCodonIns+StartCodonSNP+StopCodonDel+StopCodonIns;
    print "3UTR entries=""\t"threeutr;
    print "5Flank entries=""\t"fiveflank;
    print "5UTR entries=""\t"fiveutr;
    print "IGR entries=""\t"igr;
    print "Intron entries=""\t"intron;
    print "De_novo_Start_InF000rame entries=""\t"DenovoStartInFrame;
    print "De_novo_Start_OutOfFrame entries=""\t"DenovoStartOutOfFrame;
    #print "Non-coading transcript entries=""\t"noncoading;
    print "Frame_Shift_Del entries=""\t"FrameShiftDel;
    print "Frame_Shift_Ins entries=""\t"FrameShiftIns;
    print "In_Frame_Del entries=""\t"InFrameDel;
    print "In_Frame_Ins entries=""\t"InFrameIns;
    print "lincRNA entries=""\t"lincRNA;
    print "RNA entries=""\t"RNA;
    print "Start_Codon_Del entries=""\t"StartCodonDel;
    print "Start_Codon_Ins entries=""\t"StartCodonIns;
    print "Start_Codon_SNP entries=""\t"StartCodonSNP;
    print "Stop_Codon_Del entries=""\t"StopCodonDel;
    print "Stop_Codon_Ins entries=""\t"StopCodonIns;
    print "Missense mutation entries=""\t"missense;
    print "Nonsense mutation entries=""\t"nonsense;
    print "Nonstop entries=""\t"nonstop;
    print "Silent entries=""\t"silent;
    print "Splice entries=""\t"splice"\n"
    }' < $f >> /media/LUN1/NIlesh_analysis/Individual_tumor/all_mutations.info
done

## Extracting Cosmic and novel information for individual tumor
mkdir /media/LUN1/NIlesh_analysis/Individual_tumor/Novel_cosmic/
for f in $(find /media/LUN1/NIlesh_analysis/Individual_tumor/* -name "*.txt" -type f ); 
do
  sample=$(basename $f);
	sample_name="${sample%.*}"
awk -F'\t'  '{if($121 ~ /-/ && $257 ~ /-/ && $168 ~ /-/) print $0; 
              else if($121 !~ /-/ && $257 ~ /-/ && $168 ~ /-/) print $0;
              else if($121 !~ /-/ && $257 !~ /-/ && $168 ~ /-/) print $0;
              else if($121 !~ /-/ && $257 ~ /-/ && $168 !~ /-/) print $0;
              else if($121 !~ /-/ && $257 !~ /-/ && $168 !~ /-/) print $0;
             }' < $f >> /media/LUN1/NIlesh_analysis/Individual_tumor/Novel_cosmic/$sample_name"_novel_cosmic.txt"
done


###################################### Into cosmic and novel we make table of all types of mutations calculations
for f in $(find /media/LUN1/NIlesh_analysis/Individual_tumor/Novel_cosmic/* -name "*.txt" -type f ); 
do
  sample=$(basename $f);
	sample_name="${sample%.*}"
awk -v file="$sample_name" -F'\t' 'BEGIN{missense=0; nonsense=0; nonstop=0; splice=0; threeutr=0; fiveflank=0; fiveutr=0; igr=0; intron=0; silent=0; DenovoStartInFrame=0; DenovoStartOutOfFrame=0; FrameShiftDel=0; FrameShiftIns=0; InFrameDel=0; InFrameIns=0; lincRNA=0; RNA=0; StartCodonDel=0; StartCodonIns=0; StartCodonSNP=0; StopCodonDel=0; StopCodonIns=0;} 
    {if($9 ~ "3" && $9 ~ "UTR") threeutr++;
    else if($9 ~ "5" && $9 ~ "Flank") fiveflank++;
    else if($9 ~ "5" && $9 ~ "UTR") fiveutr++;
    else if($9 == "IGR") igr++;
    else if($9 == "Intron") intron++;
    else if($9 == "De_novo_Start_InFrame") DenovoStartInFrame++;
    else if($9 == "De_novo_Start_OutOfFrame") DenovoStartOutOfFrame++;
    #else if($9 == "Non-coding_Transcript") noncoading++;
    else if($9 == "Frame_Shift_Del") FrameShiftDel++;
    else if($9 == "Frame_Shift_Ins") FrameShiftIns++;
    else if($9 == "In_Frame_Del") InFrameDel++;
    else if($9 == "In_Frame_Ins") InFrameIns++;
    else if($9 == "lincRNA") lincRNA++;
    else if($9 == "RNA") RNA++;
    else if($9 == "Start_Codon_Del") StartCodonDel++;
    else if($9 == "Start_Codon_Ins") StartCodonIns++;
    else if($9 == "Start_Codon_SNP") StartCodonSNP++;
    else if($9 == "Stop_Codon_Del") StopCodonDel++;
    else if($9 == "Stop_Codon_Ins") StopCodonIns++;
    else if($9 == "Missense_Mutation") missense++;
    else if($9 == "Nonsense_Mutation") nonsense++;
    else if($9 == "Nonstop_Mutation") nonstop++;
    else if($9 == "Silent") silent++;
    else if($9 == "Splice_Site") splice++;}
    END{print "Sample Name=""\t"file
    print "Total number of variants=""\t"missense+nonsense+nonstop+splice+threeutr+fiveflank+fiveutr+igr+intron+noncoading+silent+DenovoStartInFrame+DenovoStartOutOfFrame+FrameShiftDel+FrameShiftIns+InFrameDel+InFrameIns+lincRNA+RNA+StartCodonDel+StartCodonIns+StartCodonSNP+StopCodonDel+StopCodonIns;
    print "3UTR entries=""\t"threeutr;
    print "5Flank entries=""\t"fiveflank;
    print "5UTR entries=""\t"fiveutr;
    print "IGR entries=""\t"igr;
    print "Intron entries=""\t"intron;
    print "De_novo_Start_InF000rame entries=""\t"DenovoStartInFrame;
    print "De_novo_Start_OutOfFrame entries=""\t"DenovoStartOutOfFrame;
    #print "Non-coading transcript entries=""\t"noncoading;
    print "Frame_Shift_Del entries=""\t"FrameShiftDel;
    print "Frame_Shift_Ins entries=""\t"FrameShiftIns;
    print "In_Frame_Del entries=""\t"InFrameDel;
    print "In_Frame_Ins entries=""\t"InFrameIns;
    print "lincRNA entries=""\t"lincRNA;
    print "RNA entries=""\t"RNA;
    print "Start_Codon_Del entries=""\t"StartCodonDel;
    print "Start_Codon_Ins entries=""\t"StartCodonIns;
    print "Start_Codon_SNP entries=""\t"StartCodonSNP;
    print "Stop_Codon_Del entries=""\t"StopCodonDel;
    print "Stop_Codon_Ins entries=""\t"StopCodonIns;
    print "Missense mutation entries=""\t"missense;
    print "Nonsense mutation entries=""\t"nonsense;
    print "Nonstop entries=""\t"nonstop;
    print "Silent entries=""\t"silent;
    print "Splice entries=""\t"splice"\n"
    }' < $f >> /media/LUN1/NIlesh_analysis/Individual_tumor/Novel_cosmic/Cosmic_novel_all_mutations.info
done