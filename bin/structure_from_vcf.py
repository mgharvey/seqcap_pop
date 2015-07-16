#!/usr/bin/env python

"""
Name: structure_from_vcf.py
Author: Michael G. Harvey
Date: 16 July 2015
Description: This script produces a STRUCTURE input file from a vcf file output using GATK
using the seqcap_pop pipeline (https://github.com/mgharvey/seqcap_pop). It includes linkage
information (the distance between SNPs at the same locus) and phase probabilities in cases
where multiple SNPs fall on the same locus. The phase line indicates the probability that 
the phase is correct relative to the previous site (use set MARKOVPHASE=1) for a given 
individual. Sometimes STRUCTURE doesn't recognize the newline characters from this script, 
but if you copy and paste the text from the file into an existing text file this often 
corrects the issue.

"""

import os
import sys
import argparse


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_file",
			type=str,
			help="""The input vcf file from freebayes"""
		)
	parser.add_argument(
			"out_file",
			type=str,
			help="""The output fasta file"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	infile = open("{0}".format(args.in_file), 'r')
	outfile = open("{0}".format(args.out_file), 'wb')
	lines = infile.readlines()
	dists = list()
	snps = 0
	prev_locus = None
	for line in lines:
		if not line.startswith("#"):
			locus = line.split('|')[0]
			pos = line.split()[1]
			outfile.write("{0}_{1}\t".format(locus, pos))
			if locus == prev_locus:
				dists.append("{0}".format(int(pos)-int(prev_pos)))
			else:
				dists.append("-1")
			prev_locus = locus
			prev_pos = pos
			parts = line.split()
			individuals = parts[9:len(parts)]
	outfile.write("\n")
	for dist in dists:
		outfile.write("{0}\t".format(dist))
	outfile.write("\n")
	prev_locus = None
	for i in range(len(individuals)):
		genotypes = list()
		phases = list()
		snps = 0
		for line in lines:
			if line.startswith("#CHROM"):
				parts = line.split()
				this_ind = parts[i+9]
			if not line.startswith("#"):
				parts = line.split()
				locus = line.split('|')[0]
				ind = parts[i+9]
				ind_parts = ind.split(':')
				calls = ind_parts[0].split('/')
				if '.' in calls:
					genotypes.append(["-9", "-9"])					
					phases.append("0.5")	
				else:
					if "HP" not in parts[8]:
						genotypes.append([calls[0], calls[1]])
						phases.append("0.5")	
					elif "HP" in parts[8]:
						phasing = ind_parts[4].split(',')
						if '.' in phasing:
							genotypes.append([calls[0], calls[1]])
							phases.append("0.5")	
						else:
							if phasing[0].endswith("-1"):
								genotypes.append([calls[0], calls[1]])
								if locus == prev_locus:
									phases.append("1.0")
								else:
									phases.append("0.5")
							elif phasing[0].endswith("-2"):
								genotypes.append([calls[1], calls[0]])
								if locus == prev_locus:
									phases.append("1.0")
								else:
									phases.append("0.5")
				prev_locus = locus
				snps += 1		
		outfile.write("{0}\t".format(this_ind))
		for genotype in genotypes:			
			outfile.write("{0}\t".format(genotype[0]))
		outfile.write("\n")		
		outfile.write("{0}\t".format(this_ind))
		for genotype in genotypes:			
			outfile.write("{0}\t".format(genotype[1]))
		outfile.write("\n")		
		for phase in phases:
			outfile.write("{0}\t".format(phase))
		outfile.write("\n")	
		print "{0}: {1} SNPs processed".format(this_ind, len(genotypes))
	infile.close()
	outfile.close()
									
if __name__ == '__main__':
    main()	
	
