#!/usr/bin/python
#-*- coding: UTF-8

#Script written by Pierrick Wainschtein : w.pierrick@gmail.com
#This script takes a raw file (type: 20160912_STAR_BOScalid_clipped_fichier_travail_RawData.txt) as input and performs 3 operations.
#First the script filters all the SNPs except the synonymous and non synonymous ones and writes them in another file.
#The second function look the pairs of SNPs located in close proximity (one or two base away) and store these pairs in another file.
#The last function look, for each pair of SNP, if the SNPs are present in the same mutant or not. If two adjacent SNPs dont have any common mutant, the pair will be discarded

#example with 5 mutants and a wild type: 20160912_STAR_BOScalid_clipped_fichier_travail_RawData.txt 6 

import sys
import re
import os
import csv

def syn_non_syn_filter (input_file):
#This function takes a raw input file a and select only the synonymous and non synonymous SNPS to output them in another file called a_SNSFiltered.txt
	with open(input_file,'r') as data:
		with open (input_file+'_SNSFiltered.txt', 'w+') as adding_data:
			print "Running script, filtering Synonymous and Non Synonymous..."
			lignes=data.readlines()
			for i in lignes:
			#regular search for the terms for each line and write them in another file
				if i.find("\tNON_SYNONYMOUS\t") > 0 or i.find("\tSYNONYMOUS\t") > 0 :
					adding_data.write(i)

				elif i.find("effet_set") > 0:
					adding_data.write(i)
				else:
					pass

	cove_SNS_file = input_file+'_SNSFiltered.txt'
	return cove_SNS_file
	print "Filtering done"
	data.close()

def SNP_sorting (cove_SNS_file):
#Use the SNP filtered list as input and returns a file containing SNPs Synonymous or Non Synonymous that are consecutive
	regex_all = '^((chr.*)\t(.*))$'
	regex = re.compile(regex_all)
	print "Running script, selecting nearby SNPs:"
	with open(cove_SNS_file,'r') as SNP_list:
		mutmut_list = open('SNS_double_SNP_no_same_mutant_filter.txt', 'w')
		lines=SNP_list.readlines()[1:]
		pos = 0
		pos_next = 3
		double_mut_count = 0
		pos_dual = []
		chrom = 1
		countline = 0
		for i in lines :
			countline += 1
			locallist = [int(countline)-2, int(countline)-1]
			pos_next  =  int(re.search('^chr[0-9]+_([0-9]+)_SNP\t(.*)$', i).group(1))
			chrom_next = int(re.search('^chr([0-9]+)_[0-9]+_SNP\t(.*)$', i).group(1))
			if chrom_next > chrom:
			#Check the chromosme so SNPs number are consecutive on the same chromosome
				pos = 0
				chrom = chrom_next
			else:
				if pos_next == pos + 1 or pos_next == pos + 2 :
					print "DOUBLE MUTATION DETECTEE\t" + str(pos) + "\t" + str(pos_next) + "\t" + str(chrom) + "\t" + str(chrom_next)
					pos = pos_next
					double_mut_count = double_mut_count + 1
					pos_dual.append(locallist)
				elif pos_next < pos:
					pos = 0
				else:
					pos = pos_next
		flattened = [val for sublist in pos_dual for val in sublist]
		#flattens the list (a list [[SNP, chr...], [SNP, chr,...]] ) within a list)
		for i in flattened:
			mutmut_list.write(lines[i])
#		return mutmut_list



def sample_count(sample_number):	
#This function takes the SNP list syn and non-syn, close together, and check if the pair of mutations is called for the same mutant. If not, it will discard the mutant
	with open ('double_mutations_input_second_script.txt', 'w') as mutmut_mutants, open ('SNS_double_SNP_no_same_mutant_filter.txt') as mut:
		genotypes = []
		print "Running script, selecting double mutations among same mutant..."
		genotypes_encode = genotypes
		dblemut_line_number = []
		indicelist = 0
		for i in csv.reader(mut, delimiter='\t'):
		#select only the columns with reference / alt and mutants alleles to compare them
			genotypes.append(i[4:int(sample_number)+6])
		print "Number of pairs detected :\t" + str(len(genotypes)/2)
		#print genotypes
		if len(genotypes)%2 !=0:
			print "WARNING, DOUBLE MUTATION LIST NOT EVEN NUMBER"

		for i in genotypes:
		#encode all the genotypes according to the reference allele in a list of [ 0 1 0 1 ] for each SNP and compares the list 2 per 2 (coocuring SNP)
			indicelist += 1
			for j in range(2,int(sample_number)+2):
				#print str(i)+ "\t" + str(i[0])+ str(i[j])
				a = str(i[0])
				if str(i[j]) == str(a)+"/"+str(a):
					genotypes_encode[indicelist-1][j] = 1 #1 non mute
				else:
					genotypes_encode[indicelist-1][j] = 0 #0 mute


		for i,(l1, l2) in enumerate(zip(genotypes[0::2], genotypes[1::2])):
			a = 0
			for j in range(2,int(sample_number)+2):
				similar_mut_count = 0
				#print l1[j], l2[j]
				if l1[j] == l2[j] == 0:
					a = 1
			if a == 1:
			#appends the line if the mutation is in the same mutant
				dblemut_line_number.append(i*2+1)
				dblemut_line_number.append(i*2+2)
			else:
				pass

		#print dblemut_line_number
		mut.close()
	with  open ('double_mutations_input_second_script.txt', 'w+') as mutmut_mutants, open('SNS_double_SNP_no_same_mutant_filter.txt', 'r') as mutants:
		lines = mutants.readlines()
		for i in dblemut_line_number:
			mutmut_mutants.write(lines[i-1])
		return mutmut_mutants

def coovar_input(sample_number):
#This function takes the list of double SNPs from previous function and outputs a file specifically formated for CooVar (CHRnumber Position Gene Ref Alt Mutants etc...) called "double_mutations_input_coovar.vcf"
	f = open('double_mutations_input_second_script.txt').readlines()

	f = [i.strip('\n').split() for i in f]

	new_data = []

	for i in f:
		data1 = i[0].split("_")
		new = data1[0][-1]+" "+data1[1]+" "

		new += i[1]+" "

		new += ' '.join(i[4:])
		new = new.replace(" ", "\t")

		new_data.append(new)

	with open('double_mutations_input_coovar.vcf', 'w') as fp:
		writer = csv.writer(fp, delimiter='\n')
		writer.writerow(new_data)

		
	
def main(input_file, mutant_number):
#Main, calling the functions
	first_output = syn_non_syn_filter(input_file)
	second_output = SNP_sorting(first_output)
	third_output = sample_count(mutant_number)
	fourth_output = coovar_input(mutant_number)
	
input_file = sys.argv[1]
mutant_number = sys.argv[2]
main(input_file, mutant_number)
print "Script over, all files outputed. Please use 'double_mutations_input_coovar.vcf' with Coovar and 'double_mutations_input_coovar.vcf' with second script."