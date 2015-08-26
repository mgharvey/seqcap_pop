#!/usr/bin/env python

"""
Name: gphocs_from_phy.py 
Author: Michael G. Harvey
Date: 1 June 2015

Description: Convert a bunch of phylip alignments into the input file format for G-PhoCS.

Usage: 	python gphocs_from_phy.py in_dir out_file

"""

import os
import sys
import argparse
from Bio import AlignIO


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_dir",
			type=str,
			help="""The input directory of nexus files"""
		)
	parser.add_argument(
			"out_file",
			type=str,
			help="""The desired output file name"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	os.chdir(args.in_dir)
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles: # Remove hidden files
		if not prefile.startswith('.'):
			files.append(prefile)
	out = open("{0}".format(args.out_file), 'wb')
	out.write("{0}\n".format(len(files)))
	for i, file in enumerate(files):
		print "{0}{1}".format(args.in_dir, file)
		alignment = AlignIO.read("{0}/{1}".format(args.in_dir, file), "phylip-relaxed")
		out.write("locus{0} {1} {2}\n".format(i+1, len(alignment), alignment.get_alignment_length()))
		for record in alignment:
			out.write("{0} {1}\n".format(record.id, record.seq))	
		out.flush()
	out.close()
	
if __name__ == '__main__':
    main()