INTRODUCTION
-------

seqcap_pop is a series of commands and scripts we have been using for processing sequence capture data from population-level samples using tools within the [Phyluce](http://phyluce.readthedocs.org/en/latest/index.html#) package (Faircloth 2015). This pipeline includes some additional scripts not available from Phyluce and dependencies (in the "bin" folder).

LICENSE
-------

The code within this repository is available under a 3-clause BSD license. See the License.txt file 
for more information.

CITATION
--------

If you use this pipeline for your own research, please cite Phyluce following the instructions here: 

	http://phyluce.readthedocs.org/en/latest/citing.html

And provide a link to this repository for the additional scripts:

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

### 1.	Clean raw reads (Illumiprocessor)

Follow the instructions in the [Illumiprocessor documentation](http://illumiprocessor.readthedocs.org/en/latest/index.html) to make a configuration (.conf) file. Then execute Illumiprocessor. The command will look something like this:

```
illumiprocessor --input /path/to/1_raw-reads --output /path/to/2_clean-reads \
	--config illumiprocessor.conf --cores 8
```

### 2.	Assemble reads into contigs (e.g., VelvetOptimiser)

For population-level studies, I will typically use multiple individuals to make my assembly. This requires inputting reads from all of the individuals into the assembler simultaneously. You could also make your assembly from just a single individual (perhaps the one with the most reads, or the one which provides the best reference for some biological reason). More individuals require more memory (I typically have to run this step on high-memory cluster nodes). It may be worth trying different strategies for reference assembly and comparing your assemblies to determine which method produces the best results (most contigs mapping to probes, longest contigs, etc.).

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

Next, we map the reads back to the contigs to obtain a pileup. This is the most finicky step of the pipeline, as BWA often fails. I have 3 different versions of BWA installed, and if a sample repeatedly fails in one version, I switch to a different version (algorithms differ slightly across versions). I typically have success with BWA 0.7.4, 0.7.3, or 0.7.0. 

For individual 1:

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

```
samtools view -bS /path/to/5-mapping/Genus_species_1-aln.sam \
	> /path/to/5-mapping/Genus_species_1-aln.bam
```

Then do the same thing for individual 2 (and any additional individuals).

### 6.	Clean the .bam file (Picard)

This step soft-clips reads at the end of the reference contigs.

```
java -jar ~/anaconda/jar/CleanSam.jar \
	I=/path/to/5-mapping/Genus_species_1-aln.bam \
	O=/path/to/6_picard/Genus_species_1-aln_CL.bam \
	VALIDATION_STRINGENCY=SILENT
```

Then do the same thing for individual 2 (and any additional individuals).

### 7.	Add read groups (Picard)

```
java -Xmx2g -jar ~/anaconda/jar/AddOrReplaceReadGroups.jar \
    I=/path/to/6_picard/Genus_species_1-aln_CL.bam \
    O=/path/to/6_picard/Genus_species_1-aln_RG.bam \
    SORT_ORDER=coordinate \
    RGPL=illumina \
    RGPU=TestXX \
    RGLB=Lib1 \
    RGID=Genus_species_1 \
    RGSM=Genus_species_1 \
    VALIDATION_STRINGENCY=LENIENT
```

Then do the same thing for individual 2 (and any additional individuals).

### 8.	Mark PCR duplicate reads (Picard)

```
java -Xmx2g -jar ~/anaconda/jar/MarkDuplicates.jar \
    I=/path/to/6_picard/Genus_species_1-aln_RG.bam \
    O=/path/to/6_picard/Genus_species_1-aln_MD.bam \
    METRICS_FILE=/path/to/6_picard/Genus_species_1.metrics \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=false
```

Then do the same thing for individual 2 (and any additional individuals).

### 9.	Merge the BAM files across individuals within your species (Phyluce)
    
```
java -Xmx2g -jar ~/anaconda/jar/MergeSamFiles.jar \
    SO=coordinate \
    AS=true \
    I=/path/to/6_picard/Genus_species_1-aln_MD.bam \
    I=/path/to/6_picard/Genus_species_2-aln_MD.bam \
    O=/path/to/7_merge-bams/Genus_species.bam 
```

### 10.	Index the merged .bam file (samtools)

```
samtools index /Volumes/G-DRIVE/Amazon/7_merge-bams/Genus_species.bam 
```

### 11.	Create a dictionary from the reference contigs (Picard)

```
java -Xmx2g -jar ~/anaconda/pkgs/picard-1.106-0/jar/CreateSequenceDictionary.jar \
    R=/path/to/4_match-contigs-to-probes/Genus_species.fasta \
    O=/path/to/4_match-contigs-to-probes/Genus_species.dict

```

### 12.	Index the reference (samtools)

```
samtools faidx /path/to/4_match-contigs-to-probes/Genus_species.fasta
```

### 13.	Call indels (GATK)

```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -I /path/to/7_merge-bams/Genus_species.bam  \
    --minReadsAtLocus 7 \
    -o /path/to/8_GATK/Genus_species.intervals
```

### 14.	Realign indels (GATK)

```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -I /path/to/7-merge/Genus_species.bam  \
    -targetIntervals /path/to/8_GATK/Genus_species.intervals \
    -LOD 3.0 \
    -o /path/to/8_GATK/Genus_species_RI.bam
```
    
### 15.	Call SNPs (GATK)

```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -I /path/to/8_GATK/Genus_species_RI.bam \
    -gt_mode DISCOVERY \
    -o /path/to/8_GATK/Genus_species_raw_SNPs.vcf \
    -ploidy 2 \
    -rf BadCigar
```

### 16.	Annotate SNPs (GATK)
    
```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -I /path/to/8_GATK/Genus_species_RI.bam \
    -G StandardAnnotation \
    -V:variant,VCF /path/to/8_GATK/Genus_species_raw_SNPs.vcf \
    -XA SnpEff \
    -o /path/to/8_GATK/Genus_species_SNPs_annotated.vcf \
    -rf BadCigar      
```
  
### 17.	Annotate Indels (GATK)

```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -I /path/to/8_GATK/Genus_species_RI.bam \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -o /path/to/8_GATK/Genus_species_SNPs_indels.vcf \
    -rf BadCigar         
```
    
### 18.	Mask indels (GATK)

```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -V /path/to/8_GATK/Genus_species_raw_SNPs.vcf \
    --mask /path/to/8_GATK/Genus_species_SNPs_indels.vcf \
    --maskExtension 5 \
    --maskName InDel \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "BadValidation" \
    --filterExpression "QUAL < 30.0" \
    --filterName "LowQual" \
    --filterExpression "QD < 5.0" \
    --filterName "LowVQCBD" \
    -o /path/to/8_GATK/Genus_species_SNPs_no_indels.vcf  \
    -rf BadCigar
```
    
### 19.	Restrict to high-quality SNPs (bash)

```
cat /path/to/8_GATK/Genus_species_SNPs_no_indels.vcf | grep 'PASS\|^#' > /path/to/8_GATK/Genus_species_SNPs_pass-only.vcf 
```

### 20.	Read-backed phasing (GATK)

The output of this step is a .vcf file containing the final set of phased SNPs for all individuals. This can then be used to produce input files for SNP-based analyses. The subsequent two steps (21 and 22) also produce output containing these same SNPs, but separated into a file for each individual. 

```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -T ReadBackedPhasing \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -I /path/to/8_GATK/Genus_species_RI.bam \
    --variant /path/to/8_GATK/Genus_species_SNPs_pass-only.vcf \
    -L /path/to/8_GATK/Genus_species_SNPs_pass-only.vcf \
    -o /path/to/8_GATK/Genus_species_SNPs_phased.vcf \
    --phaseQualityThresh 20.0 \
    -rf BadCigar
```

### 21.	Make a vcf for each sample (GATK)

```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -T SelectVariants \
    --variant /path/to/8_GATK/Genus_species_SNPs_phased.vcf \
    -o /path/to/8_GATK/Genus_species_1_SNPs.vcf \
    -sn Genus_species_1 \
    -rf BadCigar
```

Then do the same thing for individual 2 (and any additional individuals).

### 22.	Make a table of phased SNPs for each sample (GATK)
   
```
java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R /path/to/4_match-contigs-to-probes/Genus_species.fasta \
    -V /path/to/8_GATK/Genus_species_1_SNPs.vcf \
    -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
    -o /path/to/9_SNP-tables/Genus_species_1_SNPs_phased-table.txt \
    -rf BadCigar
```

Then do the same thing for individual 2 (and any additional individuals).

### 23.	Add phased SNPs to reference and optionally filter (seqcap_pop script)

Because GATK does not output sequences or SNP calls for invariant sites, we need to add the SNPs back into the reference sequences if we want alignments. This is not optimal, as sequence data will be present at invariant sites even if no reads inform that site for a given individual. However, this is currently the only way I can find to take advantage of read-backed phasing and get the phased calls into alignment format. Freebayes allows output of invariant and variant sites, which can then easily be used to produce complete alignments. However, the phasing options in Freebayes are much more limited, and there is no obvious way to use Freebayes to output data once SNPs have been phased in GATK. If you want to go the Freebayes route, a script to obtain alignments from .vcf Freebayes output is [here](https://github.com/mgharvey/misc_python/blob/master/bin/freebayes_vcf2fa.py). This part of the pipeline could obviously use some improvement!

```
python add_phased_snps_to_seqs_filter.py \
	/path/to/4_match-contigs-to-probes/Genus_species.fasta \
	/path/to/9_SNP-tables/Genus_species_1_SNPs_phased-table.txt \
	/path/to/10_sequences/Genus_species \
	Genus_species_1_sequences.txt \
	1
```

Then do the same thing for individual 2 (and any additional individuals).

The final argument ("1") above filters out any alleles not supported by a particular number of reads. This can be increased in order to set a hard filter on the minimum number of reads for allele calls (generally not recommended).

### 24.	Collate sequences from all individuals into files by UCE (seqcap_pop script)

```
python collate_sample_fastas_GATK.py \
	/path/to/10_sequences/Genus_species/ \
	/path/to/11_fasta-parts/Genus_species/ \
	sequences.txt
```

### 25.	Align the sequences (MAFFT)

Although the sequences should line up, this serves to guarantee that there won't be any issues.

```
python run_mafft.py \
	/path/to/11_fasta-parts/Genus_species/ \
	/path/to/12_raw-alignments/Genus_species/
```

### 26.	Process the alignments (seqcap_pop script)

This makes phylip alignments from the raw fasta files from MAFFT.

```
python process_mafft_alignments_GATK.py \
	/path/to/12_raw-alignments/Genus_species \
	/path/to/13_processed-phylip/Genus_species
```

### File format conversion

For programs that don't accept .vcf files of SNPs or phylip sequence alignments, there are some additional scripts in the "bin" folder that permit data format conversion (I will be adding more here as I write them):

- **structure_from_vcf.py** should produce files of all SNPs containing linkage and phasing information for use in STRUCTURE.

```
python structure_from_vcf.py \
	/path/to/8_GATK/Genus_species_SNPs_phased.vcf \
	/path/to/desired/output/directory/Genus_species_STRUCTURE.txt
```

- **faststructure_from_vcf.py** should produce files of all SNPs, the first SNP per locus, or a random SNP from each locus for use in fastSTRUCTURE.

```
python faststructure_from_vcf.py \
	/path/to/8_GATK/Genus_species_SNPs_phased.vcf \
	/path/to/desired/output/directory/ \
	Genus_species \
	--all \
	--first \
	--random
```

- **gphocs_from_phy.py** reformats alignments for input into [G-PhoCS](http://compgen.bscb.cornell.edu/GPhoCS/).

```
python gphocs_from_phy.py \
	/path/to/13_processed-phylip \
	/path/to/desired/output/directory/Genus_species_GPHOCS.txt
```
