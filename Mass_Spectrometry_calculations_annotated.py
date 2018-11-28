#!/usr/bin/env python

import sys
from math import sqrt

# purpose: this algorithm calculates peak area for each peptide type i.e. all species of K9 in a sample then outputs individual replicates peak area. It then averages and provides standard deviation of all replicates in the same sample_group for groups of up to 3 technical/biological replicates 

#usage: python Automating_MSMS_calculations_annotated.py outputfile_from_skyline.tsv

"""
	
	****************************** Automating MS/MS relative abundance calculations ******************************
				     		    				by: Dean Sanders
								   				  July 18, 2017    
				   				University of Wisconsin - Genetics Department   

"""

def get_DA():
	dictA={}

	with open(sys.argv[1], 'rb') as f:
		next(f)
		rows = (line.split('\t') for line in f)
		val = "precursor"
		for row in rows:
		#print row[13].strip()
			if row[13].strip() == val:
				dictA[row[9]] = row[1:]
			else:
				pass #print val, row[13].strip()
	return dictA

### your samples in groups: this algorithm calculates peak area for each peptide type i.e. all species of K9 then outputs individual rep peak area as and averages provides stdev of all replicates in your groups up to 3 technical/biological reps ###

sample_groups = [["WT0-5-1.raw","WT0-5-2.raw"],["WT4-5-1.raw","WT4-5-2.raw","WT4-5-3.raw"],["hda2-0-1.raw","hda2-0-2.raw","hda2-0-3.raw"],["hda2-4-1.raw","hda2-4-2.raw","hda2-4-3.raw"],["hda2-12-1.raw","hda2-12-2.raw","hda2-12-3.raw"]] #Fill in your sample groups here

### setup parameters, list of your peptide notes (they must be spelled exactly the same as they are output)

H3K4 = ['H3K4un','H3K4me1','H3K4me2','H3K4me3']  # H3K4ac is another possibility
H3K9 =  ['H3K9un_K14un','H3K9ac_K14un','H3K9un_K14ac','H3K9ac_K14ac','H3K9me1_K14un','H3K9me2_K14un']  # ,'H3K9me3_K14un' ,'H3K9me1_K14ac',
H3K18 = ['H3K18un_K23un','H3K18ac_K23ac','H3K18un_K23me1'] # 'H3K18un_K23ac','H3K18ac_K23un', commented out for this analysis
H3_1K27K36 = ['H3.1K27un_K36un','H3.1K27un_K36me1','H3.1K27me1_K36un','H3.1K27me1_K36me1','H3.1K27me2_K36un','H3.1K27un_K36me2','H3.1K27me3_K36un', 'H3.1K27me2_K36me1']#, 'H3.1K27me1_K36me2']
H3_3K27K36 =  ['H3.3K27un_K36un' ,'H3.3K27un_K36me1','H3.3K27me1_K36un','H3.3K27me1_K36me1','H3.3K27me2_K36un','H3.3K27un_K36me2','H3.3K27me3_K36un','H3.3K27un_K36me3','H3.3K27me2_K36me1','H3.3K27me1_K36me2']																												
H3K79 = ['H3K79un', 'H3K79ac','H3K79me1','H3K79me2']
H4 = ['H4_un','H4_1ac','H4_2ac','H4_3ac','H4_4ac'] 

# use only precursors of +2 or +3 for each peptide type i.e. K9. If you use both you'll keep get IndexErrors with what appear to be properly formatted names 

Peptide_list = [H3K4,H3K9,H3K18,H3_1K27K36,H3_3K27K36,H4]  			# you can modify to include any peptide you desire

Precursor_charge = [2,3]

def match_sample(s,precurs_charge,dictA):					# find correct sample 
	d = {}
	val = "precursor"
	for key in dictA:							# if precursor charge 2 or 3, sample name and "precursor" are all correct append to d various data
		if s.strip() == dictA[key][7].strip() and int(dictA[key][6]) == int(precurs_charge) and dictA[key][12].strip() == val:  
						d[key] = s, dictA[key][7].strip(), dictA[key][6].strip(), int(precurs_charge), dictA[key][13].strip(),dictA[key][15].strip()
	return d
	
def calc_total_pep_area(Pep_list,match_sample_out):
	total_peptide_area = 0								# start total peptide area at 0 for each peptide
	for l in match_sample_out:								
		for p in Pep_list:  							# For this particular peptide i.e. H3K4, analyze subtypes 'H3K4un','H3K4me1','H3K4me2','H3K4me3' to calc total pep area
			if match_sample_out[l][5] == p:  				# match correct from sample data to the peptide list
				total_peptide_area += int(match_sample_out[l][4].strip())   # add the area to that samples total peptide area for H3K4
	return total_peptide_area   							# spits out a list in the exact order that your samples were fed in of total peptide areas for each sample, then peptide bc higher level iteration

def final_it(s,total_peptide_area,match_samp_out): 				# function designed to aggregate our data: sample, total peptide area, and each individual sample and gather them into one list
	nlist = []
	for m in match_samp_out:
		stuff = []
		calc = float(match_samp_out[m][4])/float(total_peptide_area)  # calculates each samples relative abundance i.e. total H3K4un/H3K4_total_area_for_all_peptides in one replicate
		stuff = s, total_peptide_area, match_samp_out[m][0],match_samp_out[m][1], match_samp_out[m][2], match_samp_out[m][3], match_samp_out[m][4],calc, match_samp_out[m][5]
		nlist.append(stuff)
	return nlist

def aggreg3(rep1,rep2,rep3,Pep_list,precurs):
	listC = rep1
	listC.extend(rep2)
	listC.extend(rep3)						# combine your three replicates into one list by extending
	sample_name = "a"						# gets overwritten, but I need to define the sample name here to call it in the print statement
	for p in Pep_list:
		counter = 0
		devlist = []
		n_count = 0
		for i in listC:
			if i[8] == p :
				devlist.append(i[7]) 		# calc std dev
				counter += i[7]				# calculate total of your replicates by adding their peak areas at a particular peptide, used below in print statement for average
				n_count += 1				# number of samples for debugging function
				sample_name = i[0]  		# continuously overwritten until final sample is completed
			else:
				pass
		differences = [x - (float(counter)/float(3)) for x in devlist]  # remember to change nums for different nums of reps Dean
		sq_diffs = [ d ** 2 for d in differences]
		ssd = sum(sq_diffs)
		variance = float(ssd)/float(2)
		sd = sqrt(float(variance))
		try:
			print sample_name,"\t",devlist[0],"\t",devlist[1],"\t",devlist[2],"\t","sample_count=","\t",n_count,"\t","total=","\t",counter,"\t","avg=","\t",float(counter)/float(3),"\t","stdev=","\t",sd,"\t",p,"\t","precurs=","\t",precurs
		except IndexError:  # try except loop added 7/10/17, catches exceptions when you don't have a devlist because there is no sample for that peptide (i.e. you entered a peptide w/mod in the peptidelist that is not present in the data)
			pass
def aggreg2(rep1,rep2,Pep_list,precurs):
	listC = rep1
	listC.extend(rep2)
	sample_name = "a"
	for p in Pep_list:
		counter = 0
		devlist = []
		n_count = 0
		for i in listC:
			if i[8] == p:
				devlist.append(i[7]) 
				counter += i[7]
				n_count += 1
				sample_name = i[0]
			else:
				pass
		differences = [x - (float(counter)/float(2)) for x in devlist]  # remember to change nums for different nums of reps Dean
		sq_diffs = [ d ** 2 for d in differences]
		ssd = sum(sq_diffs)
		variance = float(ssd)/float(1)
		sd = sqrt(float(variance))
		try:
			print sample_name,"\t",devlist[0],"\t",devlist[1],"\t","no_sample3","\t","sample_count=","\t",n_count,"\t","total=","\t",counter,"\t","avg=","\t",float(counter)/float(2),"\t","stdev=","\t",sd,"\t",p,"\t","precurs=","\t",precurs
		except IndexError:
			pass
def aggreg1(rep1,Pep_list,precurs):
	listC = rep1
	sample_name = "a"
	for p in Pep_list:
		counter = 0
		devlist = []
		n_count = 0
		for i in listC:
			if i[8] == p:
				devlist.append(i[7])   # bc only one sample is appended only one devlist entry exists
				counter += i[7]
				n_count += 1
				sample_name = i[0]
			else:
				pass
		try:
			print sample_name,"\t",devlist[0],"\t","no_sample2","\t","no_sample3","\t","sample_count=","\t",n_count,"\t","total=","\t",counter,"\t","avg=","\t",counter,"\t","stdev=","\t","null_one_sample","\t",p,"\t","precurs=","\t",precurs
		except IndexError:
			pass
for precurs in Precursor_charge:			# Iterate through both precursor charges (2,3)
	for peptide in Peptide_list:			# for each precursor charge iterate through all peptides types in the Peptide_list , defined in global variable above
		for sg in sample_groups:			# Do this for every sample in the sample group, defined in global variable
				a = len(sg)					# count the number of samples in that group and aggregate data based on this
				if a == 3:
					try:
						aggreg3(final_it(sg[0],calc_total_pep_area(peptide,match_sample(sg[0],precurs,get_DA())),match_sample(sg[0],precurs,get_DA())), final_it(sg[1],calc_total_pep_area(peptide,match_sample(sg[1],precurs,get_DA())),match_sample(sg[1],precurs,get_DA())), final_it(sg[2],calc_total_pep_area(peptide,match_sample(sg[2],precurs,get_DA())),match_sample(sg[2],precurs,get_DA())),peptide, precurs)
					except ZeroDivisionError:
						pass
				elif a == 2:
					try:
						aggreg2(final_it(sg[0],calc_total_pep_area(peptide,match_sample(sg[0],precurs,get_DA())),match_sample(sg[0],precurs,get_DA())), final_it(sg[1],calc_total_pep_area(peptide,match_sample(sg[1],precurs,get_DA())),match_sample(sg[1],precurs,get_DA())),peptide,precurs)
					except ZeroDivisionError:
						pass
				elif a == 1:
					try:
						aggreg1(final_it(sg[0],calc_total_pep_area(peptide,match_sample(sg[0],precurs,get_DA())),match_sample(sg[0],precurs,get_DA())),peptide,precurs)
					except ZeroDivisionError:
						pass
	
