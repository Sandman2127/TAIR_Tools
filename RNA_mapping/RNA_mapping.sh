#!/bin/bash
PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/mnt/ws/progs/bowtie2-2.1.0/:/mnt/ws/progs/tophat-2.0.8b/:/mnt/ws/progs/cufflinks-2.1.1/:/mnt/ws/home/dsanders/Data/SRA_toolkit/bin
export PATH
# remember you need a genome fasta file in the same directory with same index name such as TAIR10.Refseq.sm.fa
# also your .gff/.gtf file must have same chromosome label: both chr1 or Chr1

# USE BELOW TO GET DATA FROM SRA DATABANK
# fastq-dump.2.4.5 --split-files -A SRR1013878

# FOR MAPPING ALL .FASTQS IN FOLDER DO: 

/mnt/ws/progs/tophat-2.0.8b/tophat2 --transcriptome-only -m 1 -p 8 -o ${1%.fastq} -G TAIR10_GFF3_genes_transposons_chr.gff TAIR10.Refseq.sm $1 ;

#--no-novel-juncs   # for RNA mapping just to transcriptome already known
#--no-novel-indels
#--no-gtf-juncs
#-T/--transcriptome-only                    (map only to the transcriptome)
#-m/--splice-mismatches         <0-2>  default 0

/mnt/ws/progs/cufflinks-2.1.1/cufflinks -p 8 -o ${1%.fastq}_cl ${1%.fastq}/accepted_hits.bam ;

# -G/--GTF  :  quantitate against reference transcript annotation 
# Above option will not assemble novel transcripts, and the program will ignore alignments not structurally compatible with any reference transcript.
 
#  -g/--GTF-guide : use reference transcript annotation to guide assembly
# above option will include all reference transcripts as well as any novel genes and isoforms that are assembled. 



# Next do cuffmerge, unfortunately cannot be run at the same time as the entire script

#cuffmerge -g TAIR10_GFF3_genes_transposons.gff -p 8 assemblies.txt   
 
#assemblies.txt looks like:
#		./XZ10_cl/transcripts.gtf
#		./XZ11_cl/transcripts.gtf
#		./XZ12_cl/transcripts.gtf
#		./XZ13_cl/transcripts.gtf
#		./XZ14_cl/transcripts.gtf

# automatically outputs in the current directory a directory and filed called merged_asm/merged.gtf which is your merged transcriptome between all samples use it in cuffdiff


# finally cuffdiff on your samples

#cuffdiff -o K36_M36_Het_diff_exp -b Tair10_whole_genome_lgChr.fas -p 8 -L K36_2_6,M36_2_6 -u merged_asm/merged.gtf ./XZ11_ATGTCA/accepted_hits.bam,./XZ12_CCGTCC/accepted_hits.bam ./XZ13_GTCCGC/accepted_hits.bam,./XZ14_GTGAAA/accepted_hits.bam   
#must have a genome fasta on hand and the merged.gtf from cuffmerge

