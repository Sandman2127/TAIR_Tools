#!/bin/bash
PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/mnt/ws/progs/bowtie2-2.1.0/:/mnt/ws/home/dsanders/Data/SRA_toolkit/bin:/mnt/ws/home/dsanders/Data/IGVTools:
export PATH
# remember you need a genome fasta file in the same directory with same index name such as TAIR10.Refseq.sm.fa

/mnt/ws/progs/bowtie2-2.1.0/bowtie2 -p 4 --qc-filter --local -D 15 -R 2 -N 1 -L 22 -i S,1,1.75 -X 1000 -x TAIR10_chrs -1 $1 -2 $2 -S ${1%.fq.gz}.sam  ;  # from paired end read

# notes: This particular alignment is for calling snps, thus we want to be loose with our alignment parameters
# I' using a custom mapping, i.e. --local to allow end trimming if incorrect match, -D 15 -R 2 -N 0 -L 20 -i S,1,1.75 -X 1000
# Presets:                 Same as:
#  For --end-to-end:
#   --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
#   --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
#   --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
#   --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

#  For --local:
#   --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
#   --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
#   --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
#   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

# Alignment:
#  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
#  -L <int>           length of seed substrings; must be >3, <32 (22)
#  -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)
#  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)

# sort and index bams 

samtools view -bS ${1%.fq.gz}.sam  > ${1%.fq.gz}_NOSORT.bam ;

samtools sort ${1%.fq.gz}_NOSORT.bam ${1%.fq.gz}.bam ;

samtools index ${1%.fq.gz}.bam ;

samtools flagstat ${1%.fq.gz}.bam > ${1%.fq.gz}.flagstat.txt ;

rm ${1%.fq.gz}_NOSORT.bam ${1%.fq.gz}.sam ;


#index fasta:
samtools faidx TAIR10_chrs.fa

#SNP calling steps:


#samtools mpileup -uD -f /mnt/ws/home/dsanders/Data/Sequencing_data/05_10_16_JJ_SNPs_Highalt_A_thaliana/processed_Neq1/TAIR10_chrs.fa $1 | bcftools view -bvcg - > ${1%.bam}.bcf
# -b binary calling format instead of vcf standard out from bcftools
# -v variant sites only
# -g call genotype
# -c (call snps and indels)


#vcfutils.pl qstats ${1%.bam}.bcf  # give you stats of all snps you called

#vcfutils.pl varFilter 1hasa2.vcf > 1hasa2_filter_stock.vcf &

# filters a variety of parameters to ensure they are high enough quality


# prep for annovar I also needed some sed magic to get my files to have 1 2 and 3 instead of chr1, chr2 ...

# perl ~/Data/annovar/convert2annovar.pl  -format vcf4 1hasa2_hom_snps_ens.vcf > 1hasa2.anvinput

# prepping background libraries, both files obtained from 
# wget ftp://ftp.ensemblgenomes.org/pub/release-27/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.27.gtf.gz

#~/Data/gtfToGenePred -genePredExt Arabidopsis_thaliana.TAIR10.27.gtf AT_refGene.txt

#perl ~/Data/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile Arabidopsis_thaliana.TAIR10.27.dna.genome.fa AT_refGene.txt --out AT_refGeneMrna.fa

# the database and files must be named exactly as in kai wangs tutorial... I've moved all to the annovar folder.

#perl ~/Data/annovar/annotate_variation.pl --geneanno --buildver AT ./processed_Neq0/1hasa2.anvinput ./atdb  



# using grep to find the total numbers of all mutation types in exonic.variants file

#echo $1 ;

#echo "total lines" 

# wc -l $1 ; 

# echo "frameshift insertion " ;

# awk ' {print $2"\t"$3}' $1 | grep "^frameshift" | grep "insertion" | wc -l ; 

# echo "frameshift deletion" ;

# awk ' {print $2"\t"$3}' $1 | grep "^frameshift" | grep "deletion" | wc -l ;

# echo "frameshift block substitution" ;

# grep "frameshift" $1 | grep "block" | wc -l ;

# echo "stopgain" ;

# grep "stopgain" $1 | wc -l ;

# echo "stoploss" ;

# grep "stoploss" $1 | wc -l ;

# echo "nonframeshift insertion" ;

# grep "nonframeshift" $1 | grep "insertion" | wc -l ;

# echo "nonframeshift deletion" ;

# grep "nonframeshift" $1 | grep "deletion" | wc -l ;

# echo "nonframeshift block substitution " ;

# grep "nonframeshift" $1 | grep "block" | wc -l ;

# echo "nonsynonymous SNV" ;

# grep "nonsynonymous" $1 | grep "SNV" | wc -l ;

# echo "synonymous SNV" ;

# awk ' {print $2"\t"$3}' $1 | grep "^synonymous" | grep "SNV" | wc -l  ;

# echo "unknown"

# grep "unknown" $1 | wc -l ;

