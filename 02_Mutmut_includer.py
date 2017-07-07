#!/usr/bin/python
#-*- coding: UTF-8

#Script by Pierrick Wainschtein : w.pierrick@gmail.com
#This script takes as inputs: the output of CooVar "categorized-gvs.gvf" and the output of the first Python script 01_Mutmut_finder.py "double_mutations_input_second_script.txt" and compares the amino acid predicted by CooVar to the amino acid predicted (if predicted) by Refiner Genome. It outputs a warning if there's a difference between the amino acid predicted by refiner genome and the amino acid predicted by CooVar allowing for easy spotting of wrong amino acid prediction due to co-occuring variants. Please note that it works only for Non-synonymous SNPs as no prediction is being made by Refiner Genome for Synonymous SNPs.

#Example command line :  python 02_Mutmut_includer.py categorized-gvs.gvf double_mutations_input_second_script.txt

import sys
import re
import os
import csv

def real_codon_read (input_file):
#This function takes the output of Coovar and creates a dictionnary (variable : output) in the following format : {'chr5_1344087_SNP': ['A', 'A'], 'chr9_1750183_SNP': ['N', 'P'],...}. The first amino acid is the "original" false one, without taking co-occuring variants into account, and the second one is the truth one.
	SNP_falsecodon = []
	SNP_truecodon = []
	global output
	
	with open(input_file,'r') as data:
		a = 0
		SNP_list = []
		for i in csv.reader(data, delimiter='\t'):
			if a < 2:
				a = a+1
			else:
				SNP_list.append('chr'+str(i[0])+'_'+str(i[3])+'_SNP')
				a = a+1


	with open(input_file,'r') as data:
		a = 0
		b = 1
		for i in csv.reader(data, delimiter='>'):
			if a < 2:
				a = a+1
			else:
				SNP_falsecodon.append(i[1][-1])
				SNP_truecodon.append(i[2][:1])
				a = a+1

	SNP_dict_false = dict(zip(SNP_list, SNP_falsecodon))
	SNP_dict_true = dict(zip(SNP_list, SNP_truecodon))
	output = dict((k, [SNP_dict_false[k], SNP_dict_true.get(k)]) for k in SNP_dict_false)
	output.update((k, [None, SNP_dict_true[k]]) for k in SNP_dict_true if k not in SNP_dict_false)
	#print output

def aminoacid_addandcompare(mutants_list):
#This function compares the predictions from Coovar organised in a dictionnary in "output" variable with the predictions from Refiner Genome. For a non synonymous SNP, if the prediction is different, is prints the SNP position and predictions in the console. For 2 synonymous SNPs, it warns there's 2 synonymous nearby.
	#First comparing cooccuring with one non synonymous SNP (any case)
	print "\n List of SNP with different prediction, if a SNP is not paired, it is with a synonymous one:"
	with open (mutants_list) as mutmut_mutants:
		for j in csv.reader(mutmut_mutants, delimiter='\t'):
			i = j[:6-1]
			if i[2] == "NON_SYNONYMOUS":
				regexall = re.compile('\->([A-Z])')
				RGpredictedaa = (regexall.findall(i[3]))
				for key in output:
					if str(key) == i[0]:

						aacoovar = output.get(key)
						if not aacoovar[1] == str(RGpredictedaa[0]):
							print "Different prediction detected for: \t " + i[0] + "\t Prediction Refiner Genome:\t" + str(RGpredictedaa[0]) + "\t" + "Prediction CooVar:\t" + aacoovar[1]
	#Now dealing with 2 occuring synonymous SNPs
	print "\n List of double Synonymous SNPs:"
	mutmut_mutants = open(mutants_list, "r").readlines()
	# The number of lines
	n = len(mutmut_mutants)
	# Loop two by two (end at n-1 or n-2)
	for i in range(0, n-1, 2):
		first_line = mutmut_mutants[i]
		second_line = mutmut_mutants[i+1]

		syn1 = first_line.split("\t")[2]
		syn2 = second_line.split("\t")[2]
		if syn1 == syn2 and syn1 == "SYNONYMOUS":
			snp1 = first_line.split("\t")[0]
			snp2 = second_line.split("\t")[0]
			for key in output:
				if str(key) == snp1:
					aacoovar = output.get(key)
					print "WARNING, 2 consecutive SYNONYMOUS SNPs:\t" + snp1 + "\t" + snp2 + "\t Prediction Coovar :" + aacoovar[1]
			

def main(input_file):
	first_output = real_codon_read(input_file)
	second_output = aminoacid_addandcompare(mutants_list)

	
input_file = sys.argv[1]
mutants_list = sys.argv[2]

main(input_file)