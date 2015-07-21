#!/usr/bin/env python

"""
Name: faststructure_from_vcf.py
Author: Michael G. Harvey
Date: 21 July 2015

Description: This script produces a an input file in the STRUCTURE format accepted by 
fastSTRUCTURE from a vcf file output using GATK with the seqcap_pop pipeline (https://github.com/mgharvey/seqcap_pop). Unlike the structure_from_vcf script, this 
script DOES NOT include linkage information (the distance between SNPs at the same locus) 
and phase probabilities, because these are not allowed by fastSTRUCTURE. It may be 
preferable to output rather than all SNPs, only the first SNP or a random SNP from each
locus. These subsets can be output using the "--first" and "--random" arguments.

Usage: python structure_from_vcf.py in_file out_file --all --first --random

"""

import os
import sys
import argparse
import random


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_file",
			type=str,
			help="""The input vcf file from freebayes"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The directory where output files will be placed"""
		)
	parser.add_argument(
			"prefix",
			type=str,
			help="""The prefix for output files (e.g. "Genus_species")"""
		)
	parser.add_argument(
			"--all",
			action="""store_true""",
			help="""Output file with all SNPs"""
		)
	parser.add_argument(
			"--first",
			action="""store_true""",
			help="""Output file with first SNP per locus"""
		)
	parser.add_argument(
			"--random",
			action="""store_true""",
			help="""Output file with random SNP per locus"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	infile = open("{0}".format(args.in_file), 'r')
	lines = infile.readlines()
	all_genotype_matrix = list()
	first_genotype_matrix = list()
	random_genotype_matrix = list()
	sites = list() # For first locus only
	all_sites = list()
	first_sites = list()
	random_sites = list()
	site_index = 0
	prev_locus = None
	first_line = True
	for line in lines:
		if not line.startswith("#"):
			locus = line.split('|')[0]
			pos = line.split()[1]
			all_sites.append("{0}_{1}\t".format(locus, pos))
	for line in lines:
		if line.startswith("#CHROM"):
			individuals = line.split()[9:]			
		elif not line.startswith("#"):
			locus = line.split('|')[0]
			
			# For every locus
			if locus != prev_locus and first_line != True:				

				# If all
				if args.all == True:
					for site in sites: # select all sites
						genotypes = list()
						for ind in site:
							if '.' in ind:
								genotypes.append(["-9", "-9"])					
							else:
								genotypes.append([ind[0], ind[1]])
						all_genotype_matrix.append(genotypes)
				
				# If first
				if args.first == True:
					site = sites[0] # only select first site
					genotypes = list()
					for ind in site:
						if '.' in ind:
							genotypes.append(["-9", "-9"])					
						else:
							genotypes.append([ind[0], ind[1]])
					first_genotype_matrix.append(genotypes)
					first_sites.append(all_sites[(site_index-(len(sites)))])

				# If random
				if args.random == True:
					rand = random.randint(0, len(sites)-1) # randomly choose a site
					site = sites[rand]
					genotypes = list()
					for ind in site:
						if '.' in ind:
							genotypes.append(["-9", "-9"])					
						else:
							genotypes.append([ind[0], ind[1]])
					random_genotype_matrix.append(genotypes)
					random_sites.append(all_sites[(site_index-(len(sites))+rand)])

				sites = list()
			
			# For every row
			parts = line.split()
			ind_parts = parts[9:]
			inds = list()
			for ind_part in ind_parts:
				ind_parts = ind_part.split(':')
				inds.append(ind_parts[0].split('/'))
			sites.append(inds)	
			first_line = False
			prev_locus = locus				
			site_index += 1
					
	# Make output
	if args.all == True:
		allfile = open("{0}{1}_STRUCTURE_all.txt".format(args.out_dir, args.prefix), 'wb')
		#for all_site in all_sites:
		#	allfile.write("{0}\t".format(all_site))
		#allfile.write("\n")		
		for i, individual in enumerate(individuals):
			allfile.write("{0}\t.\t.\t.\t.\t.\t".format(individual))
			for genotyped_locus in all_genotype_matrix:	
				ind_genotype = genotyped_locus[i]		
				allfile.write("{0}\t".format(ind_genotype[0]))
			allfile.write("\n")		
			allfile.write("{0}\t.\t.\t.\t.\t.\t".format(individual))
			for genotyped_locus in all_genotype_matrix:	
				ind_genotype = genotyped_locus[i]		
				allfile.write("{0}\t".format(ind_genotype[1]))
			allfile.write("\n")		
		print "{0} SNPs in all-SNPs file".format(len(all_genotype_matrix))
			
	if args.first == True:
		firstfile = open("{0}{1}_STRUCTURE_first.txt".format(args.out_dir, args.prefix), 'wb')
		#for first_site in first_sites:
		#	firstfile.write("{0}\t".format(first_site))
		#firstfile.write("\n")		
		for i, individual in enumerate(individuals):
			firstfile.write("{0}\t.\t.\t.\t.\t.\t".format(individual))
			for genotyped_locus in first_genotype_matrix:	
				ind_genotype = genotyped_locus[i]		
				firstfile.write("{0}\t".format(ind_genotype[0]))
			firstfile.write("\n")		
			firstfile.write("{0}\t.\t.\t.\t.\t.\t".format(individual))
			for genotyped_locus in first_genotype_matrix:	
				ind_genotype = genotyped_locus[i]		
				firstfile.write("{0}\t".format(ind_genotype[1]))
			firstfile.write("\n")		
		print "{0} SNPs in first-SNPs file".format(len(first_genotype_matrix))

	if args.random == True:
		randomfile = open("{0}{1}_STRUCTURE_random.txt".format(args.out_dir, args.prefix), 'wb')
		#for random_site in random_sites:
		#	randomfile.write("{0}\t".format(random_site))
		#randomfile.write("\n")		
		for i, individual in enumerate(individuals):
			randomfile.write("{0}\t.\t.\t.\t.\t.\t".format(individual))
			for genotyped_locus in random_genotype_matrix:	
				ind_genotype = genotyped_locus[i]		
				randomfile.write("{0}\t".format(ind_genotype[0]))
			randomfile.write("\n")		
			randomfile.write("{0}\t.\t.\t.\t.\t.\t".format(individual))
			for genotyped_locus in random_genotype_matrix:	
				ind_genotype = genotyped_locus[i]		
				randomfile.write("{0}\t".format(ind_genotype[1]))
			randomfile.write("\n")		
		print "{0} SNPs in random-SNPs file".format(len(random_genotype_matrix))

	infile.close()
									
if __name__ == '__main__':
    main()	
	
