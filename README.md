INTRODUCTION
-------

seqcap_pop is a series of commands and scripts for processing sequence capture data from population-level samples using tools within the [Phyluce](http://phyluce.readthedocs.org/en/latest/index.html#) package (Faircloth 2015). This pipeline was used to obtain datasets for a number of papers:

- Harvey et al. in prep.
- Oswald et al. in prep.
- Sanchez et al. in prep. 

LICENSE
-------

The code within this repository is available under a 3-clause BSD license. See the License.txt file 
for more information.

CITATION
--------

If you use scripts from this repository for your own research, please provide the link to this software repository in your manuscript:

    https://github.com/mgharvey/seqcap_pop

USAGE
--------

Install Phyluce and dependencies following the instructions [here](http://phyluce.readthedocs.org/en/latest/index.html).

### 1.	Clean Raw Reads

Follow the instructions in the [Illumiprocessor documentation](http://illumiprocessor.readthedocs.org/en/latest/index.html) to make a configuration (.conf) file. The command will look something like this:

```
illumiprocessor --input /path/to/folder/of/raw/reads --output /path/to/folder/of/clean/reads --config illumiprocessor.conf --cores 8
```

### 2.	Assemble reads into contigs

This can be done using any of various assemblers (see Phyluce documentation). I typically use VelvetOptimiser, which is not part of Phyluce, but can be obtained [here](http://bioinformatics.net.au/software.velvetoptimiser.shtml). An example VelvetOptimiser command is:

```
VelvetOptimiser.pl --s 65 --e 75 --optFuncKmer=n50 --c=tbp -t 2 -a -f '-fastq.gz -shortPaired /path/to/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ1.fastq.gz /path/to/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ2.fastq.gz /path/to/Genus_species_2/split-adapter-quality-trimmed/Genus_species_2-READ1.fastq.gz /path/to/Genus_species_2/split-adapter-quality-trimmed/Genus_species_2-READ2.fastq.gz -short /path/to/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ-singleton.fastq.gz /path/to/Genus_species_2/split-adapter-quality-trimmed/Genus_species_2-READ-singleton.fastq.gz'
```

Other assemblers, however, should also work fine.

### 2.	Map contigs to probes

You next need to see which of the contigs that were assembled correspond to the loci targeted in your sequence capture procedure. You do this by mapping the contigs to the sequences from the probes you used:

```
python extract_uce_bypass_MGH.py \
	/Volumes/G-DRIVE/Amazon/3-Velvet_out/Campephilus_melanoleucos/contigs.fa \
	/Volumes/G-DRIVE/Amazon/uce-and-exon-probes.fasta \
	/Volumes/G-DRIVE/Amazon/4-match_to_probes/Campephilus_melanoleucos.fasta \
	/Volumes/G-DRIVE/Amazon/4-match_to_probes/Campephilus_melanoleucos.lastz
```

