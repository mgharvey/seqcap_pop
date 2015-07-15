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

Install Phyluce and dependencies following the instructions [here](http://phyluce.readthedocs.org/en/latest/index.html). The custom Python scripts required for this pipeline are in the "bin" folder, so also download those and put them in your project folder. Make sure you have plenty of available hard drive space (depending on read counts, you may need upwards of 2GB for each sample). You will probably want to organize the output of each step below into a series of folders. Within my project folder, I typically set up a series of output folders labeled with consecutive numbers followed by a brief text descriptor:

- 1_raw-reads
- 2_clean-reads
- 3_velvet-output
- 4_match-contigs-to-probes
- 5_mapping
- 6_picard
- 7_merge-bams
- 8_GATK
- 9_SNP-tables
- 10_sequences
- 11_fasta-parts
- 12_raw-alignments
- 13_processed-phylip

You can make all of these folder now.

### 1.	Clean Raw Reads

Follow the instructions in the [Illumiprocessor documentation](http://illumiprocessor.readthedocs.org/en/latest/index.html) to make a configuration (.conf) file. The command will look something like this:

```
illumiprocessor --input /path/to/1_raw-reads --output /path/to/2_clean-reads \
	--config illumiprocessor.conf --cores 8
```

### 2.	Assemble reads into contigs

For population-level studies, I will typically use multiple individuals to make my assembly. This requires inputting reads from all of the individuals into the assembler simultaneously. You could also make your assembly from just a single individual (perhaps the one with the most reads, or the one which provides the best reference for some biological reason). 

Any of various assembly programs can be used for this step (see Phyluce documentation). I typically use VelvetOptimiser, which is not part of Phyluce, but can be obtained [here](http://bioinformatics.net.au/software.velvetoptimiser.shtml). An example VelvetOptimiser command is:

```
VelvetOptimiser.pl --s 65 --e 75 --optFuncKmer=n50 --c=tbp -t 2 -a -f '-fastq.gz -shortPaired /path/to/2_clean-reads/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ1.fastq.gz /path/to/2_clean-reads/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ2.fastq.gz /path/to/2_clean-reads/Genus_species_2/split-adapter-quality-trimmed/Genus_species_2-READ1.fastq.gz /path/to/2_clean-reads/Genus_species_2/split-adapter-quality-trimmed/Genus_species_2-READ2.fastq.gz -short /path/to/2_clean-reads/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ-singleton.fastq.gz /path/to/2_clean-reads/Genus_species_2/split-adapter-quality-trimmed/Genus_species_2-READ-singleton.fastq.gz'
```

Other assemblers, however, should also work fine. Once you have a final contigs fasta file (this will be labeled "contigs.fa" if you used Velvet or VelvetOptimiser), put that in the "3_velvet-output" folder.

### 3.	Map contigs to probes (Phyluce)

You next need to see which of the contigs that were assembled correspond to the loci targeted in your sequence capture procedure. You do this by mapping the contigs to the sequences from the probes you used. When your assemblies are from multiple individuals, particularly in species containing high levels of genetic structure, you are likely to obtain multiple contigs from the same locus (because some individuals are very divergent from others at that genomic position). Typically, Phyluce would dispose of contigs in situations where multiple contigs map to the same locus, so you would lose a lot of loci at this step. I have tweaked a script from Brant Faircloth and Graham Derryberry (extract_uce_bypass_MGH.py) so as to save one contig in situations where multiple contigs map to a single locus. If you have used multiple individuals to make your assembly, you may want to use this script. If you have used only one, you will probably be fine using the original script (extract_uce_bypass.py). The commands will look like this:

If multiple individuals in assembly:

```
python extract_uce_bypass_MGH.py \
	/path/to/3_velvet-output/Genus_species/contigs.fa \
	/path/to/probe_sequences.fasta \
	/path/to/4_match-contigs-to-probes/Genus_species.fasta \
	/path/to/4_match-contigs-to-probes/Genus_species.lastz
```

If one individual in assembly:

```
python extract_uce_bypass.py \
	/path/to/3_velvet-output/Genus_species/contigs.fa \
	/path/to/probe_sequences.fasta \
	/path/to/4_match-contigs-to-probes/Genus_species.fasta \
	/path/to/4_match-contigs-to-probes/Genus_species.lastz
```

### 4.	Map reads to contigs (BWA)

Next, we map the reads back to the contigs to obtain a pileup. For individual 1:

```
bwa index -a is /path/to/4_match-contigs-to-probes/Genus_species.fasta
```

```
bwa aln /path/to/4_match-contigs-to-probes/Genus_species.fasta \
	/path/to/2_clean-reads/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ1.fastq.gz \ 
	> /path/to/5-mapping/Genus_species_1_read1.sa.sai
```

```
bwa aln /path/to/4_match-contigs-to-probes/Genus_species.fasta \
	/path/to/2_clean-reads/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ2.fastq.gz  \
	> /path/to/5-mapping/Genus_species_1_read2.sa.sai
```

```
bwa sampe /path/to/4_match-contigs-to-probes/Genus_species.fasta \
	/path/to/5-mapping/Genus_species_1_read1.sa.sai /path/to/5-mapping/Genus_species_1_read2.sa.sai \
	/path/to/2_clean-reads/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ1.fastq.gz \
	/path/to/2_clean-reads/Genus_species_1/split-adapter-quality-trimmed/Genus_species_1-READ2.fastq.gz \
	> /path/to/5-mapping/Genus_species_1-aln.sam
```

Then do the same thing for individual 2 (and any additional individuals).

### 5.	Convert .sam file to .bam file (samtools)

samtools view -bS /path/to/5-mapping/Genus_species_1-aln.sam \
	> /path/to/5-mapping/Genus_species_1-aln.bam

### 6.	

