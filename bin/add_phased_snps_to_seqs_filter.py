#!/usr/bin/env python

"""

Name: add_phased_snps_to_seqs.py

Author: Michael G. Harvey
Date: 1 June 2014


"""

import os
import sys
import re
import random
import argparse
import copy


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"seq_file",
			type=str,
			help="""The input file of fasta sequences"""
		)
	parser.add_argument(
			"phase_file",
			type=str,
			help="""The input file of phase data"""
		)
	parser.add_argument(
			"out_file",
			type=str,
			help="""The output file of phased sequences"""
		)
	parser.add_argument(
			"depth_filter",
			type=int,
			help="""The depth filter below which alleles will be removed"""
		)
	parser.add_argument(
			"--resolve",
			action="""store_true""",
			help="""Randomly resolve unphased SNPs?"""
		)
	return parser.parse_args()


def unphase(alleles):	
	a1 = alleles[0]
	a2 = alleles[1]
	if a1 == "A":
		if a2 == "C":
			code = "M"
		elif a2 == "G":
			code = "R"
		elif a2 == "T":
			code = "W"					
		elif a2 == "A":
			code = "A"					
	if a1 == "C":
		if a2 == "A":
			code = "M"
		elif a2 == "G":
			code = "S"
		elif a2 == "T":
			code = "Y"
		elif a2 == "C":
			code = "C"					
	if a1 == "G":
		if a2 == "A":
			code = "R"
		elif a2 == "C":
			code = "S"				
		elif a2 == "T":
			code = "K"
		elif a2 == "G":
			code = "G"					
	if a1 == "T":
		if a2 == "A":
			code = "W"
		elif a2 == "C":
			code = "Y"					
		elif a2 == "G":
			code = "K"	
		elif a2 == "T":
			code = "T"					
	if a1 == ".":
		code = "N"
	return str(code)


def main():
	args = get_args()
	seqfile = open("{0}".format(args.seq_file))
	phasefile = open("{0}".format(args.phase_file))
	sequences = seqfile.readlines()
	phasings = phasefile.readlines()
	outfile = open("{0}".format(args.out_file), 'wb')
	i = 0
	k = 0
	seq = None
	for sequence in sequences: # For each line in sequence file
		if sequence.startswith('>'): # If header			
			parts = re.split('\||\s', sequence)
			firstpart = parts[0].split('>')
			locus = str(firstpart[1]).rstrip() # Get locus name
			if i != 0:	# If not the first sequence
				seqAlist = copy.deepcopy(seq)
				seqBlist = copy.deepcopy(seq)
				prev = 'FALSE' 
				for phasing in phasings:
					parts = phasing.split()
					firstpart = parts[0].split('|')
					phaselocus = firstpart[0]
					if phaselocus == prev_locus: # If phase file locus same as sequence locus
						parts = phasing.split()
						pre_alleles = parts[3].split('/')
						depths = parts[6].split(',')
						alleles = list()
						for j, pre_allele in enumerate(pre_alleles):	
							if pre_allele != '.':
								if int(depths[j]) >= args.depth_filter:
									alleles.append(pre_allele)
									k += 1
							elif pre_allele == '.':
								alleles.append(pre_allele)	
						#print alleles
						hp = parts[5]
						if len(set(alleles)) > 1: # Is the site biallelic for this individual?
							if prev == 'TRUE': # If a biallelic site already exists at this locus
								if not hp == "NA": # If it was accurately phased
									print "Some phased SNPs at this locus"
									phasing = hp.split(',') 
									if phasing[0].endswith("-1"):
										seqAlist[int(parts[1])-1] = alleles[0] # Just insert alleles
										seqBlist[int(parts[1])-1] = alleles[1]
									elif phasing[0].endswith("-2"):
										seqAlist[int(parts[1])-1] = alleles[1] # Switch alleles
										seqBlist[int(parts[1])-1] = alleles[0]											
								elif hp == "NA": # If it was not phased
									print "Some unphased SNPs at this locus"
									if args.resolve: # Are we randomly resolving?
										rand = random.randint(0,1)
										if rand == 0:
											other = 1
										elif rand == 1:
											other = 0
										seqAlist[int(parts[1])-1] = alleles[rand]
										seqBlist[int(parts[1])-1] = alleles[other]	
									else: # If not randomly resolving, insert ambiguity code
										base = unphase(alleles)
										seqAlist[int(parts[1])-1] = base
										seqBlist[int(parts[1])-1] = base			
							elif prev == 'FALSE': # If no biallelic site already exists
								seqAlist[int(parts[1])-1] = alleles[0] # Just insert alleles
								seqBlist[int(parts[1])-1] = alleles[1]	
							prev = 'TRUE' # Register that a biallelic SNP has been found at this locus
						elif len(set(alleles)) == 1: # If the site is monomorphic for this individual				
							seqAlist[int(parts[1])-1] = alleles[0].replace('.','N') # Just insert alleles (replace errors with 'N')
							seqBlist[int(parts[1])-1] = alleles[0].replace('.','N')							
				outfile.write("{0}a".format(prev_seq))
				outfile.write("\n")
				outfile.write(''.join(seqAlist))
				outfile.write("\n")		
				outfile.write("{0}b".format(prev_seq))
				outfile.write("\n")
				outfile.write(''.join(seqBlist))
				outfile.write("\n")
				seq = None
			prev_locus = locus
			print prev_locus
			prev_seq = sequence.rstrip()
			#print locus
		else:
			if seq:
				seq = seq + list(sequence.lstrip().rstrip())
			else:
				seq = list(sequence.lstrip().rstrip())
		i += 1		
	seqAlist = copy.deepcopy(seq)
	seqBlist = copy.deepcopy(seq)
	prev = 'FALSE' 
	for phasing in phasings:
		parts = phasing.split()
		firstpart = parts[0].split('|')
		phaselocus = firstpart[0]
		if phaselocus == prev_locus: # If phase file locus same as sequence locus
			parts = phasing.split()
			pre_alleles = parts[3].split('/')
			depths = parts[6].split(',')
			alleles = list()
			for j, pre_allele in enumerate(pre_alleles):	
				if pre_allele != '.':
					if int(depths[j]) >= args.depth_filter:
						alleles.append(pre_allele)
						k += 1
				elif pre_allele == '.':
					alleles.append(pre_allele)	
			hp = parts[5]
			if len(set(alleles)) > 1: # Is the site biallelic for this individual?
				if prev == 'TRUE': # If a biallelic site already exists at this locus
					if not hp == "NA": # If it was accurately phased
						print "Some phased SNPs at this locus"
						phasing = hp.split(',') 
						if phasing[0].endswith("-1"):
							seqAlist[int(parts[1])-1] = alleles[0] # Just insert alleles
							seqBlist[int(parts[1])-1] = alleles[1]
						elif phasing[0].endswith("-2"):
							seqAlist[int(parts[1])-1] = alleles[1] # Switch alleles
							seqBlist[int(parts[1])-1] = alleles[0]											
					elif hp == "NA": # If it was not phased
						print "Some unphased SNPs at this locus"
						if args.resolve: # Are we randomly resolving?
							rand = random.randint(0,1)
							if rand == 0:
								other = 1
							elif rand == 1:
								other = 0
							seqAlist[int(parts[1])-1] = alleles[rand]
							seqBlist[int(parts[1])-1] = alleles[other]	
						else: # If not randomly resolving, insert ambiguity code
							base = unphase(alleles)
							seqAlist[int(parts[1])-1] = base
							seqBlist[int(parts[1])-1] = base			
				elif prev == 'FALSE': # If no biallelic site already exists
					seqAlist[int(parts[1])-1] = alleles[0] # Just insert alleles
					seqBlist[int(parts[1])-1] = alleles[1]	
				prev = 'TRUE' # Register that a biallelic SNP has been found at this locus
			elif len(set(alleles)) == 1: # If the site is monomorphic for this individual				
				seqAlist[int(parts[1])-1] = alleles[0].replace('.','N') # Just insert alleles (replace errors with 'N')
				seqBlist[int(parts[1])-1] = alleles[0].replace('.','N')							
	outfile.write("{0}a".format(prev_seq))
	outfile.write("\n")
	outfile.write(''.join(seqAlist))
	outfile.write("\n")		
	outfile.write("{0}b".format(prev_seq))
	outfile.write("\n")
	outfile.write(''.join(seqBlist))
	outfile.write("\n")
	print locus
	seqfile.close()
	phasefile.close()
	outfile.close()
	print k									
if __name__ == '__main__':
    main()