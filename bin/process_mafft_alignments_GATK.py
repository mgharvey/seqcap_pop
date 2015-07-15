"""

Name: process_stacks_alignments.py

Author: Michael G. Harvey
Date: 30 July 2014

Description: 

Usage:

python process_mafft_alignments.py in_dir out_dir

"""

import os 
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import argparse


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_dir",
			type=str,
			help="""The input directory of alignments"""
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
	for prefile in prefiles:
		if not prefile.startswith('.'):
			files.append(prefile)			
	for file in files:
		if os.stat("{0}/{1}".format(args.in_dir, file))[6] == 0: # Remove files with no sequences	
			print "{0} is empty".format(file)
		else:
			alignment = AlignIO.read("{0}/{1}".format(args.in_dir, file), "fasta")
			missing = ['n','-','N','?']
			data_align = "FALSE"
			new_alignment=list()
			for record in alignment:
				data_samp = "FALSE"
				for i in str(record.seq):
					if i not in missing:
						data_align = "TRUE"
						data_samp = "TRUE"
				if data_samp == "TRUE": # If data for this individual
					new_alignment.append(record)
			if data_align == "TRUE": # If some data in this alignment
				out = open("{0}/{1}.phy".format(args.out_dir, file.replace(".fa","")), 'wb')
				out.write(" {0} {1}\n".format(len(new_alignment), MultipleSeqAlignment(new_alignment).get_alignment_length()))
				for record in new_alignment:
					out.write("{0}  {1}\n".format(record.id, str(record.seq.upper())))
			else:
				print "{0} contains only missing data".format(file)
							
if __name__ == '__main__':
    main()