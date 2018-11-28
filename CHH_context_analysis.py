#!/usr/bin/env python
from __future__ import division
import pandas as pd
import argparse

#Example methratio .outy file
#chr	pos	strand	context	ratio	total_C	methy_C
#chr1	8	+	AACCC	0.000	1	0
#chr1	9	+	ACCCT	1.000	1	1
#chr1	10	+	CCCTA	1.000	2	2
#chr1	116	-	CCGGT	0.963	27	26
#chr1	117	-	CGGTT	0.185	27	5

#CHH_codons = ("CAA","CAT","CAC","CTA","CTT","CTC","CCA","CCT","CCC")

parser = argparse.ArgumentParser(prog='CHH methylation context analyser')
parser.add_argument('-in', dest='file1',action="store", type=argparse.FileType('r'),help="tab delimited .outy file BSmap methratio_2.5.py")

args = parser.parse_args()

def parsing_file():
    f = args.file1
    df = pd.read_table(f, sep='\t')  # read in with native column names 'chr','pos','strand','context','ratio','total_C','methy_C'
    return df

def linebyline(df):

# create empty dfs to dump data into
    CHH_rev_comp_codons = {"TTG":"CAA","ATG":"CAT","GTG":"CAC","TAG":"CTA","AAG":"CTT","GAG":"CTC","TGG":"CCA","AGG":"CCT","GGG":"CCC"}
    colnames=['chr','pos','strand','context','ratio','total_C','methy_C']

    CAA_df=pd.DataFrame(columns=colnames)
    CAT_df=pd.DataFrame(columns=colnames)
    CAC_df=pd.DataFrame(columns=colnames)

    CTA_df=pd.DataFrame(columns=colnames)
    CTT_df=pd.DataFrame(columns=colnames)
    CTC_df=pd.DataFrame(columns=colnames)

    CCA_df=pd.DataFrame(columns=colnames)
    CCT_df=pd.DataFrame(columns=colnames)
    CCC_df=pd.DataFrame(columns=colnames)

    for row in df.itertuples():
        chrom = row[1]
        pos = row[2]
        strand = row[3]
        context = row[4].strip()
        ratio = row[5]
        total_C = row[6]
        methyl_C = row[7]
        if strand == "+" :
            if context[2:5] == "CAA":
                CAA_df=CAA_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif context[2:5] == "CAT":
                CAT_df=CAT_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif context[2:5] == "CAC":
                CAC_df=CAC_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)

            elif context[2:5] == "CTA":
                CTA_df=CTA_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif context[2:5] == "CTT":
                CTT_df=CTT_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif context[2:5] == "CTC":
                CTC_df=CTC_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)

            elif context[2:5] == "CCA":
                CCA_df=CCA_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif context[2:5] == "CCT":
                CCT_df=CCT_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif context[2:5] == "CCC":
                CCC_df=CCC_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            else:
                pass

        elif strand == "-" :

            def Rev_comp(context):
                for key,value in CHH_rev_comp_codons.iteritems():
                    if context[0:3].strip() == key.strip():
                        rev_compl=value.strip()
                        return rev_compl

            if Rev_comp(context) == "CAA":
                CAA_df=CAA_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif Rev_comp(context) == "CAT":
                CAT_df=CAT_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif Rev_comp(context) == "CAC":
                CAC_df=CAC_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)

            elif Rev_comp(context) == "CTA":
                CTA_df=CTA_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif Rev_comp(context) == "CTT":
                CTT_df=CTT_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif Rev_comp(context) == "CTC":
                CTC_df=CTC_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)

            elif Rev_comp(context) == "CCA":
                CCA_df=CCA_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif Rev_comp(context) == "CCT":
                CCT_df=CCT_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            elif Rev_comp(context) == "CCC":
                CCC_df=CCC_df.append({'chr':chrom,'pos':pos,'strand':strand,'context':context,'ratio':ratio,'total_C':total_C,'methy_C':methyl_C},ignore_index=True)
            else:
                pass

    list_o_dfs = [CAA_df,CAT_df,CAC_df,CTA_df,CTT_df,CTC_df,CCA_df,CCT_df,CCC_df]
    return list_o_dfs

def dataframe_method(df):
    CAA_df = df[((df['context'].str.endswith('CAA')) & (df['strand'] == "+")) | ((df['context'].str.startswith('TTG')) & (df['strand'] == "-"))]
    CAT_df = df[((df['context'].str.endswith('CAT')) & (df['strand'] == "+")) | ((df['context'].str.startswith('ATG')) & (df['strand'] == "-"))]
    CAC_df = df[((df['context'].str.endswith('CAC')) & (df['strand'] == "+")) | ((df['context'].str.startswith('GTG')) & (df['strand'] == "-"))]

    CTA_df = df[((df['context'].str.endswith('CTA')) & (df['strand'] == "+")) | ((df['context'].str.startswith('TAG')) & (df['strand'] == "-"))]
    CTT_df = df[((df['context'].str.endswith('CTT')) & (df['strand'] == "+")) | ((df['context'].str.startswith('AAG')) & (df['strand'] == "-"))]
    CTC_df = df[((df['context'].str.endswith('CTC')) & (df['strand'] == "+")) | ((df['context'].str.startswith('GAG')) & (df['strand'] == "-"))]

    CCA_df = df[((df['context'].str.endswith('CCA')) & (df['strand'] == "+")) | ((df['context'].str.startswith('TGG')) & (df['strand'] == "-"))]
    CCT_df = df[((df['context'].str.endswith('CCT')) & (df['strand'] == "+")) | ((df['context'].str.startswith('AGG')) & (df['strand'] == "-"))]
    CCC_df = df[((df['context'].str.endswith('CCC')) & (df['strand'] == "+")) | ((df['context'].str.startswith('GGG')) & (df['strand'] == "-"))]

    list_o_dfs = [CAA_df,CAT_df,CAC_df,CTA_df,CTT_df,CTC_df,CCA_df,CCT_df,CCC_df]
    return list_o_dfs


def calc_perc(df):

# figure out which df this is:
    plus_df = df[df['strand'] == "+"][1:200]
    context_d = plus_df['context'].iloc[0][2:]

# calc percents
    indexlen=len(df.index)
    ratio = df['ratio'].sum()
    total_perc_rat = ratio/indexlen

    methy_C = df['methy_C'].sum()
    tot_C = df['total_C'].sum()
    total_perc = methy_C/tot_C
    
# print final output
    print context_d,"\t",ratio,"\t",indexlen,"\t",100*total_perc_rat,"\t",methy_C,"\t",tot_C,"\t",100*total_perc

if __name__ == "__main__":

    df = parsing_file()

    #data1 = linebyline(df)

    data = dataframe_method(df)

    #print "context","\t","ratio","\t","index_length","\t","total_perc_rat","\t","methylated_C","\t","total_C","\t","total_perc"

    #for i in data1:
    #    calc_perc(i)

    print "context","\t","ratio","\t","index_length","\t","total_perc_rat","\t","methylated_C","\t","total_C","\t","total_perc"

    for i in data:
        calc_perc(i)
