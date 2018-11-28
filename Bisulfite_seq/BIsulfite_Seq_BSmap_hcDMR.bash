#!/bin/bash

HOME="/mnt/ws/home/dsanders"

PLOTTER="/mnt/ws/home/dsanders/Full_genome_methylation_binner"

PATH="/mnt/ws/home/dsanders/Data/bsmap-2.90/":$PATH

export PATH PLOTTER HOME

#fastqc

#$HOME/Data/FastQC/fastqc $1

#rm *.zip

#mv *.html ./fastqc_html

#/mnt/ws/home/dsanders/sratoolkit.2.8.2/bin/fastq-dump.2.8.2 -A $3 $4 $5 $6 $7 $8 

#gzip ${3}.fastq

#mv ./${3}.fastq.gz ./${1}


/mnt/ws/home/dsanders/Data/bsmap-2.90/bsmap -f 5 -q 33 -n 0 -a $1 -d $2 -o ${1%.fq.gz}.bam -p 2 -v 0.08 ;

# v = 0.12 or 12 % of the 100 bp read can be mismatched, 0.08 = default 8% incorrect

# -f filter 5 n or less per read
# -L map first 70bp
# -q 33 trim anything not 33 quality or better
# -p  is the number of processors/threads you would like to utilize
# -v 0.08 
# -n  [0,1]   set mapping strand information. default: 0
#                   -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+), 
#                   for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.
#                   -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, -- 


python /mnt/ws/home/dsanders/Data/bsmap-2.90/methratio.py -z -r -s /mnt/ws/home/dsanders/Data/bsmap-2.90/samtools/ -o ${1%.fq.gz}.outy -d $2 ${1%.fq.gz}.bam ;
# calc total methylation levels

cat <(echo ${1%.fq.gz}.CHH.mC_level) <(awk '$4=="CHH"{MC+=$7;TO+=$8;}END{print MC/TO}'  ${1%.fq.gz}.outy) > ${1%.fq.gz}.CHH.GMC.count ;

cat <(echo ${1%.fq.gz}.CHG.mC_level) <(awk '$4=="CHG"{MC+=$7;TO+=$8;}END{print MC/TO}' ${1%.fq.gz}.outy) > ${1%.fq.gz}.CHG.GMC.count ;

cat <(echo ${1%.fq.gz}.CG.mC_level) <(awk '$4=="CG"{MC+=$7;TO+=$8;}END{print MC/TO}' ${1%.fq.gz}.outy) > ${1%.fq.gz}.CG.GMC.count ;

cat ${1%.fq.gz}.CG.GMC.count ${1%.fq.gz}.CHG.GMC.count ${1%.fq.gz}.CHH.GMC.count > ${1%.fq.gz}.GMC.count ;

rm ${1%.fq.gz}.CG.GMC.count ${1%.fq.gz}.CHG.GMC.count ${1%.fq.gz}.CHH.GMC.count ;

# separating contexts and doing WG methylation plots

for i in "CG" "CHG" "CHH" ;

do

grep $i ${1%.fq.gz}.outy > ${1%.fq.gz}.${i}.tmp ;

/mnt/ws/home/dsanders/anaconda2/bin/python2.7 $PLOTTER/WG_methylation_plotter.py $PLOTTER/genome_1Mb_binned.txt ${1%.fq.gz}.${i}.tmp > ${1%.fq.gz}.${i}.WG_plot ;

rm ${1%.fq.gz}.${i}.tmp ;

done


# methdiff DMR calling

#python /mnt/ws/home/dsanders/Data/bsmap-2.90/methdiff.py -o ${1%.fq.gz}.CG -l col,${1%.fq.gz} -x CG -d $2 -r 0.4 -b 100 -p 0.01 ${CONTROL_LIB%.fq.gz}.outy ${1%.fq.gz}.outy  ;

#python /mnt/ws/home/dsanders/Data/bsmap-2.90/methdiff.py -o ${1%.fq.gz}.CHG  -l col,${1%.fq.gz} -x CHG -d $2 -r 0.2 -b 100 -p 0.01 ${CONTROL_LIB%.fq.gz}.outy ${1%.fq.gz}.outy ;

#python /mnt/ws/home/dsanders/Data/bsmap-2.90/methdiff.py -o ${1%.fq.gz}.CHH  -l col,${1%.fq.gz} -x CHH -d $2 -r 0.1 -b 100 -p 0.01 ${CONTROL_LIB%.fq.gz}.outy ${1%.fq.gz}.outy ;


# CHH context analysis:

#$HOME/anaconda2/bin/python2.7 $HOME/Data/bsmap-2.5/methratio.py -z -u -s /mnt/ws/home/dsanders/Data/bsmap-2.90/samtools/ -o ${1%.fq.gz}.u.2.5.outy -d $2 ${1%.fq.gz}.bam ;

#$HOME/anaconda2/bin/python2.7 $HOME/Data/CHH_context_analysis.py -in ${1%.fq.gz}.u.2.5.outy >> ${1%.fq.gz}.CHH.context.txt ;


# transform for jbrowse tracks

#non stranded bGph separation

#awk '$4=="CHH"{print $1"\t"$2-1"\t"$2"\t"$5}' $1 | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${1%.fq.gz}.bGph

#stranded bGph: note the $3$5 in the first awk statement, use non-stranded bGph separation for feeding into deeptools

awk '$4=="CHH"{print $1"\t"$2-1"\t"$2"\t"$3$5}' ${1%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${1%.fq.gz}.CHH.bGph ;
awk '$4=="CHG"{print $1"\t"$2-1"\t"$2"\t"$3$5}' ${1%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${1%.fq.gz}.CHG.bGph ;
awk '$4=="CG"{print $1"\t"$2-1"\t"$2"\t"$3$5}' ${1%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${1%.fq.gz}.CG.bGph ;


$HOME/bedGraphToBigWig/bedGraphToBigWig ${1%.fq.gz}.CHH.bGph $HOME/bedGraphToBigWig/tairchrs ${1%.fq.gz}.CHH.bw ;
$HOME/bedGraphToBigWig/bedGraphToBigWig ${1%.fq.gz}.CHG.bGph $HOME/bedGraphToBigWig/tairchrs ${1%.fq.gz}.CHG.bw ;
$HOME/bedGraphToBigWig/bedGraphToBigWig ${1%.fq.gz}.CG.bGph $HOME/bedGraphToBigWig/tairchrs ${1%.fq.gz}.CG.bw ;

rm ${1%.fq.gz}.CG.bGph ;
rm ${1%.fq.gz}.CHG.bGph ;
rm ${1%.fq.gz}.CHH.bGph ;


## Run hcDMR caller

### Step 1 - generate methratio file from aligned bam file:

##### Example usage:

python /mnt/ws/home/dsanders/Data/hcDMR_caller/Main/methratio_alt.py --Steve --ref=./$2 --out=${1%.fq.gz}.out -u -z -r ${1%.fq.gz}.bam ;

#python /mnt/ws/home/dsanders/Data/hcDMR_caller/Main/methratio_alt.py --Steve --ref=./$2 --out=${1%.fq.gz}.noU.out -z -r ${1%.fq.gz}.bam ;

##### This step will output a methratio file: input_file_name.gz

### Step 2 - count the C and CT count at every position in the genome

##### Required files:
#* input_file_name.gz *output from methratio_alt.py script*
#* TAIR10_v2.cytosine.gz *can be found in the folder /Reference*

##### Input scripts:
#* BSmap_to_cytosine.pl

##### Example usage:

perl /mnt/ws/home/dsanders/Data/hcDMR_caller/Main/BSmap_to_cytosine.pl --input_file ${1%.fq.gz}.out.gz --reference_cytosine /mnt/ws/home/dsanders/Data/hcDMR_caller/Ref_data/TAIR10_v2.cytosine.gz ;

# reference has to be remade for --reference cytosine, we cannot just use TAIR10 the C's will be in the wrong places: 
#chr1
#1	+	z
#2	+	z


##### This step will output a C and CT count file: input_file_name.cytosine.gz

### Step 3 - Bin the genome into 100bp bins.

perl /mnt/ws/home/dsanders/Data/hcDMR_caller/Main/Cytosine_to_100bp.pl ${1%.fq.gz}.out.cytosine.gz ;

##### This step will output three files (types of CHH, CHG, CG methylation) containing methylation level along the genome in 100bp bin:
#* input_file_name.CHH.100.gz
#* input_file_name.CHG.100.gz input_file_name.CG.100.gz

### Step 4 - Call hcDMRs against 54 WT dataset:

perl /mnt/ws/home/dsanders/Data/hcDMR_caller/Main/hcDMR_caller_DSmod2.pl -ref /mnt/ws/home/dsanders/Data/hcDMR_caller/Ref_data/CHH.100.54WT.Ref.txt.gz -input ${1%.fq.gz}.out.CHH.100.gz -dif 0.1 -n 33 ;

mv ${1%.fq.gz}.DMR ${1%.fq.gz}.CHH.DMR ;

awk '$7>=0.10 && $2!="begin" && $5<$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${1%.fq.gz}.CHH.DMR | grep "hypo" > ${1%.fq.gz}.CHH.hypo.DMR ;
awk '$7>=0.10 && $2!="begin" && $5>$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${1%.fq.gz}.CHH.DMR | grep "hyper"  > ${1%.fq.gz}.CHH.hyper.DMR ;

perl /mnt/ws/home/dsanders/Data/hcDMR_caller/Main/hcDMR_caller_DSmod2.pl -ref /mnt/ws/home/dsanders/Data/hcDMR_caller/Ref_data/CHG.100.54WT.Ref.txt.gz -input ${1%.fq.gz}.out.CHG.100.gz -dif 0.2 -n 33 ;

mv ${1%.fq.gz}.DMR ${1%.fq.gz}.CHG.DMR  ;

awk '$7>=0.2 && $2!="begin" && $5<$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${1%.fq.gz}.CHG.DMR | grep "hypo" > ${1%.fq.gz}.CHG.hypo.DMR  ;
awk '$7>=0.2 && $2!="begin" && $5>$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${1%.fq.gz}.CHG.DMR | grep "hyper" > ${1%.fq.gz}.CHG.hyper.DMR  ;

perl /mnt/ws/home/dsanders/Data/hcDMR_caller/Main/hcDMR_caller_DSmod2.pl -ref /mnt/ws/home/dsanders/Data/hcDMR_caller/Ref_data/CG.100.54WT.Ref.txt.gz -input ${1%.fq.gz}.out.CG.100.gz -dif 0.4 -n 33 ;

mv ${1%.fq.gz}.DMR ${1%.fq.gz}.CG.DMR ;

awk '$7>=0.4 && $2!="begin" && $5<$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${1%.fq.gz}.CG.DMR | grep "hypo" > ${1%.fq.gz}.CG.hypo.DMR  ;
awk '$7>=0.4 && $2!="begin" && $5>$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${1%.fq.gz}.CG.DMR | grep "hyper" > ${1%.fq.gz}.CG.hyper.DMR  ;


#making bw tracks of DMR files:

for i in ${1%.fq.gz}.CHH.hyper.DMR ${1%.fq.gz}.CHH.hypo.DMR ${1%.fq.gz}.CHG.hyper.DMR ${1%.fq.gz}.CHG.hypo.DMR ${1%.fq.gz}.CG.hyper.DMR ${1%.fq.gz}.CG.hypo.DMR ; 

do

awk '{print $1"\t"$2"\t"$3"\t"$6$5}' $i | sort -k1,1 -k2,2n > ${i%.DMR}.bGph ;

$HOME/bedGraphToBigWig/bedGraphToBigWig ${i%.DMR}.bGph $HOME/bedGraphToBigWig/tairchrs ${i%.DMR}.bw ;

# bedGraphToBigWig is a unix executable that can be downloaded from: http://hgdownload.cse.ucsc.edu/admin/exe/

rm ${i%.DMR}.bGph ;

done


