#!/usr/bin/env python

"""

Name: run_mafft.py

Author: Michael G. Harvey
Date: 29 July 2014

Description: 

"""

import os
import sys
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
	return parser.parse_args()


def main():
	args = get_args()
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles: # Remove hidden files
		if not prefile.startswith('.'):
			files.append(prefile)
	for file in files:
		os.system("mafft --auto --adjustdirection {0}/{1} > {2}/{3}".format(args.in_dir, file, args.out_dir, file))
		
if __name__ == '__main__':
    main()