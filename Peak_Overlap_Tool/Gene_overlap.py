#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse 

# you must set the local directories of TEs "ALL_TEs.txt" and Genes "TAIR10_genes_only.bed" variables before using this program unless they are in the current directory

parser = argparse.ArgumentParser(prog='Gene Overlap Map Program')
parser.add_argument('-in', dest='file1',action="store", type=argparse.FileType('r'),help=" in tab delimited bed format i.e. chrom,start,end,type,val,strand, see script for details")
parser.add_argument('-TE_gb', dest='TE_gb',nargs='?',default=0,action="store",type=int,help=" 1 or 0")
parser.add_argument('-TE_po', dest='TE_po',nargs='?',default=0,action="store",type=int,help=" 1 or 0")
parser.add_argument('-Gene_gb', dest='Gene_gb',nargs='?',default=0,action="store",type=int,help=" 1 or 0")
parser.add_argument('-Gene_po', dest='Gene_po',nargs='?',default=0,action="store",type=int,help=" 1 or 0")
parser.add_argument('-prom_length', dest='prom_len',nargs='?',default=0,action="store",type=int,help="Length of promoter you'd like to look for")

args = parser.parse_args()

names1 = ["chrom","start","end","type","val","strand"]  # standard bed format, must only be 6 columns though

df=pd.read_table(args.file1,sep="\t",names=names1)

TEs=pd.read_table("ALL_TEs.txt",sep="\t")  

chromnames = {"AT1":"chr1","AT2":"chr2","AT3":"chr3","AT4":"chr4","AT5":"chr5","ATM":"chrM","ATC":"chrC"}

TEs['chroms']=TEs['Transposon_Name'].str.split("TE").apply(lambda x: x[0]).map(chromnames) 

gene_header = ["chrom","start","end","zero","strand","type","gene_name"]

Genes=pd.read_table("TAIR10_genes_only.bed",sep="\t",names=gene_header)

#TE file
#strands = {"False":"-","True":"+"}
#Transposon_Name	orientation_is_5prime	Transposon_min_Start	Transposon_max_End	Transposon_Family	Transposon_Super_Family
#0	AT1TE52125	False	15827287	15838845	ATHILA2	LTR/Gypsy


#Gene file
#"chrom","start","end","zero","strand","type","gene_name"
#chr1    3631    5899            0       +       gene    AT1G01010
#chr1    5928    8737            0       -       gene    AT1G01020


# -in file:
#chr1    257701  257800  hypo    0       +
#chr1    696501  696600  hypo    0       +

listA = []

# working with TEs section:

def TE_gb_overlap(chrom,start,end): 
    A = TEs[(TEs['chroms'] == str(chrom)) & (int(start) >= TEs['Transposon_min_Start']) & (int(end) <= TEs['Transposon_max_End'])]
    listA.append(A.loc[:,'Transposon_Name'].tolist())
    B = TEs[(TEs['chroms'] == str(chrom)) & (int(start) >= TEs['Transposon_min_Start']) & (int(start) <= TEs['Transposon_max_End']) & (int(end) >= TEs['Transposon_max_End'])]
    listA.append(B.loc[:,'Transposon_Name'].tolist())
    C = TEs[(TEs['chroms'] == str(chrom)) & (int(start) <= TEs['Transposon_min_Start']) & (int(end) >= TEs['Transposon_min_Start']) & (int(end) <= TEs['Transposon_max_End'])]
    listA.append(C.loc[:,'Transposon_Name'].tolist())

def TE_po_overlap(chrom,start,end,strand,length): 
    #"+ strand TRUE"
    A = TEs[(TEs['chroms'] == str(chrom)) & (TEs['orientation_is_5prime'] == bool('True')) & (int(start) >= (TEs['Transposon_min_Start']-length)) & (int(end) <= TEs['Transposon_min_Start'])]  # in the middle
    listA.append(A.loc[:,'Transposon_Name'].tolist())
    B = TEs[(TEs['chroms'] == str(chrom)) & (TEs['orientation_is_5prime'] == bool('True')) & (int(start) >= (TEs['Transposon_min_Start']-length)) & (int(start) <= TEs['Transposon_min_Start']) & (int(end) >= TEs['Transposon_min_Start'])] # overlaps prom into TE
    listA.append(B.loc[:,'Transposon_Name'].tolist())
    C = TEs[(TEs['chroms'] == str(chrom)) & (TEs['orientation_is_5prime'] == bool('True')) & (int(start) <= (TEs['Transposon_min_Start']-length)) & (int(end) >= (TEs['Transposon_min_Start']-length)) & (int(end) <= TEs['Transposon_max_End'])] # overlaps into the promoter
    listA.append(C.loc[:,'Transposon_Name'].tolist())

    #"- strand FALSE"
    D = TEs[(TEs['chroms'] == str(chrom)) & (TEs['orientation_is_5prime'] == bool('False')) & (int(start) >= TEs['Transposon_max_End']) & (int(start) <= (TEs['Transposon_max_End']+length))]  # in the middle of promoter
    listA.append(D.loc[:,'Transposon_Name'].tolist())
    E = TEs[(TEs['chroms'] == str(chrom)) & (TEs['orientation_is_5prime'] == bool('False')) & (int(start) >= TEs['Transposon_max_End']) & (int(start) <= (TEs['Transposon_max_End']+length)) & (int(end) >= (TEs['Transposon_max_End']+length))]  # into beginning of promoter - direction
    listA.append(E.loc[:,'Transposon_Name'].tolist())
    F = TEs[(TEs['chroms'] == str(chrom)) & (TEs['orientation_is_5prime'] == bool('False')) & (int(start) <= TEs['Transposon_max_End']) & (int(end) >= TEs['Transposon_max_End']) & (int(end) <= (TEs['Transposon_max_End']+length))]
    listA.append(F.loc[:,'Transposon_Name'].tolist())

# working with Genes section:

def Gene_gb_overlap(chrom,start,end): 
    A = Genes[(Genes['chrom'] == str(chrom)) & (int(start) >= Genes['start']) & (int(end) <= Genes['end'])]
    listA.append(A.loc[:,'gene_name'].tolist())
    B = Genes[(Genes['chrom'] == str(chrom)) & (int(start) >= Genes['start']) & (int(start) <= Genes['end']) & (int(end) >= Genes['end'])]
    listA.append(B.loc[:,'gene_name'].tolist())
    C = Genes[(Genes['chrom'] == str(chrom)) & (int(start) <= Genes['start']) & (int(end) >= Genes['start']) & (int(end) <= Genes['end'])]
    listA.append(C.loc[:,'gene_name'].tolist())

def Gene_po_overlap(chrom,start,end,strand,length): 
    #"+ strand TRUE"
    A = Genes[(Genes['chrom'] == str(chrom)) & (Genes['strand'] == str("plus")) & (int(start) >= (Genes['start']-length)) & (int(end) <= Genes['start'])]  # in the middle
    listA.append(A.loc[:,'gene_name'].tolist())
    B = Genes[(Genes['chrom'] == str(chrom)) & (Genes['strand'] == str("plus")) & (int(start) >= (Genes['start']-length)) & (int(start) <= Genes['start']) & (int(end) >= Genes['start'])] # overlaps prom into TE
    listA.append(B.loc[:,'gene_name'].tolist())
    C = Genes[(Genes['chrom'] == str(chrom)) & (Genes['strand'] == str("plus")) & (int(start) <= (Genes['start']-length)) & (int(end) >= (Genes['start']-length)) & (int(end) <= Genes['end'])] # overlaps into the promoter
    listA.append(C.loc[:,'gene_name'].tolist())

    #"- strand FALSE"
    D = Genes[(Genes['chrom'] == str(chrom)) & (Genes['strand'] == str("minus")) & (int(start) >= Genes['end']) & (int(start) <= (Genes['end']+length))]  # in the middle of promoter
    listA.append(D.loc[:,'gene_name'].tolist())
    E = Genes[(Genes['chrom'] == str(chrom)) & (Genes['strand'] == str("minus")) & (int(start) >= Genes['end']) & (int(start) <= (Genes['end']+length)) & (int(end) >= (Genes['end']+length))]  # into beginning of promoter - direction
    listA.append(E.loc[:,'gene_name'].tolist())
    F = Genes[(Genes['chrom'] == str(chrom)) & (Genes['strand'] == str("minus")) & (int(start) <= Genes['end']) & (int(end) >= Genes['end']) & (int(end) <= (Genes['end']+length))]
    listA.append(F.loc[:,'gene_name'].tolist())

if __name__ == '__main__':
    if args.TE_gb == 1:
        df.apply(lambda row: TE_gb_overlap(row['chrom'],row['start'],row['end']),axis=1)
    elif args.TE_po == 1:
        df.apply(lambda row: TE_po_overlap(row['chrom'],row['start'],row['end'],row['strand'],args.prom_len),axis=1)
    elif args.Gene_gb == 1:
        df.apply(lambda row: Gene_gb_overlap(row['chrom'],row['start'],row['end']),axis=1)
    elif args.Gene_po == 1:
        df.apply(lambda row: Gene_po_overlap(row['chrom'],row['start'],row['end'],row['strand'],args.prom_len),axis=1)
    else:
        print("program operations not present")  

    listB = []
    # compile all the data into a non-repeated list
    for i in listA:
        if i:
            listB.extend(i)
        else:
            pass
    for A in set(listB):        # set only prints unique values in a list 
        print(A)
