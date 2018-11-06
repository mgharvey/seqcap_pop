#!/usr/bin/env python

"""

Name: collate_sample_fastas.py

Author: Michael G. Harvey
Date: 29 July 2014

Description: 

"""

import os
import sys
import re
import argparse


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_dir",
			type=str,
			help="""The input directory"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	parser.add_argument(
			"in_end",
			type=str,
			help="""The ending of fasta files to search for"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	files = list()
	prefiles = list()
	dirfiles = os.listdir("{0}".format(args.in_dir))
	for dirfile in dirfiles:
		if dirfile.endswith(args.in_end):
			prefiles.append("{0}".format(dirfile))
	for prefile in prefiles:
		if not prefile.startswith('.'):
			files.append(prefile)
	uces = list()
	for file in files:
		file_uces = list()
		lines = open("{0}/{1}".format(args.in_dir, file), 'r')
		for line in lines:
			if line.lstrip().startswith(">"):
				header = re.split(">|\||\s", line)
				uce = header[1]
				if uce not in uces:
					uces.append(uce)
		lines.close()
	for uce in uces:
		print uce
		out = open("{0}/{1}.fa".format(args.out_dir, str(uce).rstrip()), 'wb')
		for file in files:
			prev = False
			lines = open("{0}/{1}".format(args.in_dir, file), 'r')
			for line in lines:
				if line.lstrip().startswith(">"):
					header = re.split(">|\||\s", line)
					this_uce = header[1]
					if this_uce == uce:
						uce_seq = next(lines)
						name = file.split('_')
						if prev == False:
							out.write(">{0}_{1}_{2}a\n{3}".format(name[0], name[1], name[2], uce_seq))
							prev = True
						else:
							out.write(">{0}_{1}_{2}b\n{3}".format(name[0], name[1], name[2], uce_seq))
		lines.close()
		out.close()
	
if __name__ == '__main__':
    main()