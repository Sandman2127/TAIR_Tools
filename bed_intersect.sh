#!/bin/bash

export PATH=:/mnt/ws/home/dsanders/bedtools-2.17.0/bin/:/mnt/ws/progs/bin/python2.7:/mnt/ws/home/dsanders/Data/bsmap-2.90/:/mnt/ws/home/dsanders/sratoolkit.2.8.2/bin/:/usr/bin/samtools/:$PATH


val=$(echo $1 " total is:");
val2=$(wc -l $1 | awk '{print $1}');
val3=$(echo $2 " total is:");
val4=$(wc -l $2 | awk '{print $1}');

/mnt/ws/home/dsanders/bedtools-2.17.0/bin/bedtools intersect -a $1 -b $2 > ${1%.DMR}.olp.${2%.DMR}.DMR ;

echo " overlap is:" ;
wc -l ${1%.DMR}.olp.${2%.DMR}.DMR ;

val5=$(echo "overlap of" $1 $2 "is")
val6=$(wc -l ${1%.DMR}.olp.${2%.DMR}.DMR | awk '{print $1}')  

/mnt/ws/home/dsanders/bedtools-2.17.0/bin/bedtools intersect -v -a $1 -b $2 > ${1%.DMR}.not.${2%.DMR}.DMR ;

echo  $1 " only is:" ;
wc -l ${1%.DMR}.not.${2%.DMR}.DMR ;

/mnt/ws/home/dsanders/bedtools-2.17.0/bin/bedtools intersect -v -a $2 -b $1 > ${2%.DMR}.not.${1%.DMR}.DMR ;

echo $2 " only is:" ;
wc -l ${2%.DMR}.not.${1%.DMR}.DMR ;

echo $val $val2 >> overlap_summary
echo $val3 $val4 >> overlap_summary
echo $val5 $val6 >> overlap_summary


# for changing 1s to chr1s

#awk '$1==1{print "chr1""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $1 > ${1}.1.cat
#awk '$1==2{print "chr2""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $1 > ${1}.2.cat
#awk '$1==3{print "chr3""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $1 > ${1}.3.cat
#awk '$1==4{print "chr4""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $1 > ${1}.4.cat
#awk '$1==5{print "chr5""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $1 > ${1}.5.cat

#cat ${1}.1.cat ${1}.2.cat ${1}.3.cat ${1}.4.cat ${1}.5.cat > ${1}.met.DMR

#rm ${1}.1.cat ${1}.2.cat ${1}.3.cat ${1}.4.cat ${1}.5.cat



#for i in *.bed; 

#do

#/mnt/ws/progs/bin/python2.7 /mnt/ws/home/dsanders/Peak_Overlap_tool/Gene_map_overlap.py /mnt/ws/home/dsanders/Peak_Overlap_tool/TAIR10_genes_only.bed $i -ma 500 > ${i%.bed}.500.map.txt ; 

#/mnt/ws/progs/bin/python2.7 /mnt/ws/home/dsanders/Peak_Overlap_tool/Gene_map_overlap.py /mnt/ws/home/dsanders/Peak_Overlap_tool/TAIR10_genes_only.bed $i -ma 1000 > ${i%.bed}.1000.map.txt ;

#done
 


