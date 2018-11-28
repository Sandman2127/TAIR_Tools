#!/bin/bash
PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/mnt/ws/progs/bowtie2-2.1.0/:/mnt/ws/home/dsanders/Data/SRA_toolkit/bin:/mnt/ws/home/dsanders/Data/IGVTools:
export PATH
# remember you need a genome fasta file in the same directory with same index name such as TAIR10.Refseq.sm.fa

/mnt/ws/progs/bowtie2-2.1.0/bowtie2 -p 8 --fast --no-1mm-upfront --qc-filter -x TAIR10.Refseq.sm -q ./"$1" -S ${1%.fastq}.sam  ;

samtools view -bS ${1%.fastq}.sam > ${1%.fastq}_NOSORT.bam ;  # convert sam to unsorted bam

samtools sort ${1%.fastq}_NOSORT.bam ${1%.fastq} ;  # sort bam

samtools index ${1%.fastq}.bam ;   # index bam

samtools flagstat ${1%.fastq}.bam > ${1%.fastq}.flagstat.txt ;  # counts the # of aligned reads 

rm ${1%.fastq}_NOSORT.bam ;  # remove unsorted bam

rm ${1%.fastq}.sam ;  # remove massive sam file

/mnt/ws/home/dsanders/bedtools-2.17.0/bin/bamToBed -i ${1%.fastq}.bam > ${1%.fastq}.bed ;  # convert aligned reads in bam to bed format for handling by other programs

#for i in *.bam ; do /mnt/ws/home/dsanders/Data/IGVTools/igvtools count -z 7 -w 25 -e 300 ./$i ./${i%.bam}.tdf ./Tair10.chrom.sizes ;  done

# for peak calling I used SICER and MACS2, simple 1 liners using the recommended program parameters can be used to call peaks for your data