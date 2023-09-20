# Unio_Genome_Skimming
This is just a place where I've curated resources to help with assembling and annotating genome skim reads for phylogenetic analyses. 
# Purpose #

This workflow was created to reconstruct a phylogeny of Unionoida based on mitochondrial genomes, ribosomal repeats, and UCEs extracted from genome skim data. It does so as follows, but keep note that output/input from several of these steps can be somewhat interchangeably (i.e., annotating mtgenomes in MitoFinder from SPAdes assemblies from phyluce, etc.), therefor this is not necessarily a linear workflow:

1. Initial data QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [fastp](https://github.com/OpenGene/fastp) on paired raw Illumina reads. 
2. Trim adapters from raw paired raw Illumina reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [Illumiprocessor](https://illumiprocessor.readthedocs.io/en/latest/index.html).
3. Assemble and annotate mitogenomes from clean reads using [Mitofinder](https://github.com/RemiAllio/MitoFinder).
	3a. Assemble reads using megahit or metaspades implemented in MitoFinder. 
	3b. Extract gene orders from the assembled mitogenomes (specific to Unionoida mtgenomes).
	3c. Align individual genes using [MAFFT] (https://mafft.cbrc.jp/alignment/software/) or [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/) implemented in [AliView](https://ormbunkar.se/aliview/) or other sequence viewing software.  
4. Extract ribosomal repeat regions (ITS-1, 28S, etc.) from clean read data (can also do it from mitofinder output, I think) in [Geneious](https://www.geneious.com/). 
5. Concatenate alignments (if needed/desired) or export alignments in phylip or fasta format. 
6. Construct a phylogeny from the mitogenome and ribosomal gene alignments using [IQ-TREE](http://www.iqtree.org) and [RevBayes](https://revbayes.github.io/).
7. UCE Phylogenomics implemented in phyluce

These steps have been compiled from several resources listed below:

#### [Genome_Skimming_Lab_NMNH](https://github.com/trippster08/genome_skimming_LAB) ####
#### [Genome_skimming_Lab_repository](https://github.com/trippster08/genome_skimming_LAB/archive/refs/heads/main.zip) ####
#### [PhylOcto](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH#readme) ####
#### [MitoFinder](https://github.com/RemiAllio/MitoFinder) ####
#### [phyluce](https://phyluce.readthedocs.io/en/latest/index.html) ####

# STEPS #

## 1. Initial QC

### FastQC Raw Reads

A good first step is to run FastQC on all your raw read data to check the quality and help determine trimming parameters.

Run the FastQC script [fastqc.sh](https://github.com/trippster08/genome_skimming_LAB/blob/main/scripts/fastqc.sh), including the path to the directory containing your raw read files. If you are running on hydra, the command should look something like this: /scratch/genomics/<USERNAME>/<PROJECT>/data/raw.

```
sh fastqc.sh <path_to_raw_sequences>
```

There is the [fastqc.job](https://github.com/trippster08/genome_skimming_LAB/blob/main/jobs/fastqc.job) for running on hydra. 

If you do not enter the path to the raw sequences in the command, or enter a path to a directory that does not contain fastq.gz files, you will get the following error "Correct path to read files not entered (*.fastq.gz)". You may get additional errors, but they should stem from an incorrect or missing path, so adding that path should fix these errors.

Download the directory containing the FastQC results (it should be /data/raw/fastqc_analyses) to your computer. Open the html files using your browser to examine your read quality. Interpreting FastQC results can be tricky, and will not be discussed here. See LAB staff or others familiar with fastQC for help.

## 2. Trimming adapters and data filtering

### Illumiprocessor

This has basically been taken from the [QC step](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#clean-the-read-data) in the phyluce pipeline. 

If using Illumiprocessor (which outputs files that make it convenient to run with phyluce), do the following. 

You first need to set up you directory structure to look something like this:
<PROJECT_DIR>
	->illumiprocessor.conf
	->raw-fastq
    		->22113FL-01-01-173_L001_R1_001.fastq.gz
    		->22113FL-01-01-173_L001_R2_001.fastq.gz
    		->22113FL-01-01-174_L001_R1_001.fastq.gz
    		->22113FL-01-01-174_L001_R2_001.fastq.gz

Then your illumiprocessor.conf file needs to look like the following. You will likely need the Admire SampleKey document or the document that lists the IndexID/Barcode and IndexSequence for each sample. 

This is the section where you list the adapters you used.  the asterisk will be replaced with the appropriate index for the sample.
[adapters]
i7:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

This is the list of indexes you used
[tag sequences]
iTru7_101_01:GGTAACGT
iTru7_102_01:CGCCTTAT

This is how each index maps to each set of reads
[tag map]
22113FL-01-01-173:iTru7_101_01
22113FL-01-01-174:iTru7_102_01

You will want to rename your read files to something a bit more nice
This will rename 22113FL-01-01-173 and 22113FL-01-01-174 to whatever you want to use going forward. I recommend structuring like <Genus_species_accession>
[names]
22113FL-01-01-173:<new_sampleID>
22113FL-01-01-174:<new_sampleID>

Then run illumiprocessor

### fastp

If using fastp, do the following. 

Run the fastp script [fastp.sh](https://github.com/trippster08/genome_skimming_LAB/blob/main/scripts/fastp.sh), including the path to the directory containing your raw read files. For most, it should be something like: /scratch/genomics/<USERNAME>/<PROJECT>/data/raw.

```
sh fastp.sh <path_to_raw_sequences>
```

There is also the [fastp.job](https://github.com/trippster08/genome_skimming_LAB/blob/main/jobs/fastp.job) file that will run this on hydra. 

### Trimmomatic

Last, there is [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), but I recommend first trying the two methods before this. 

Use [adapter_trim.sh](https://github.com/alexfranzen/Unio_Genome_Skimming/blob/main/adapter_trim.sh) by supplying the path to the reads (a directory holding files named things like "[sampleID]**R1**[other stuff]" for forward reads and "[sampleID]**R2**[other stuff]" for reverse reads), and the path to a fasta file of your adapters (I think most of the genome skimming library preps so far use any of the TruSeq adapters).
Be sure to change the path (and version) of Trimmomatic if yours differs.

```
bash trim_adapters.sh dir-with-reads adapters-fasta
```

Sometimes the files will be named weirdly (ex. fastq.gz is in the sample name twice. Use command mv *[redundant_characters] [new_name] to batch rename files)

```
for f in *.fastq.gz_paired.fastq.gz; do mv -- "$f" "${f%.fastq.gz_paired.fastq.gz}_paired.fastq.gz"; done
for f in *.fastq.gz_unpaired.fastq.gz; do mv -- "$f" "${f%.fastq.gz_unpaired.fastq.gz}_unpaired.fastq.gz"; done
```

Any of these methods should sufficiently clean and trim all of you reads and make them read for the next steps. 

## 3. Assembling and annotating mitogenomes

You can do both assembly and annotation, or just assembly and vice versa, in [MitoFinder]((https://github.com/RemiAllio/MitoFinder)). It is convenient in some respects to do everything in MitoFinder, but I think you get somewhat better results if you do the assemblies outside of MitoFinder (This will come into play with phyluce). [SPAdes](https://github.com/ablab/spades) is a good assembler, but mitofinder doesn't run very well when you choose metaspades as the assembler. There is a way to take your SPAdes scaffolds and annotate them in MitoFinder (discussed later). Anyway, below are the commands and steps for both assembly and annotation. 

### mtgenome reference

A good first step is to get your mitochondrial genome reference sequences. Mitofinder requires the references to be in GenBank format (.gb), but this is really easy. 
First go to Genbank (https://www.ncbi.nlm.nih.gov/genbank/) and type in “Unionidae male complete mitochondrion” or “Unionidae female complete mitochondrion” in the search bar. 
Find and select the sequences you want to include in the reference file. I recommend keeping the male and female mtgenomes separate. 
Then click on “Send to” at the top right and choose “Complete Record”, then “File” under “Choose Destination”, and “GenBank (full)” under “Format” and hit the “Create File” button.
Then download it as a .gb file to your computer, and if working on hydra or HPC upload to project directory.

### 3a. assembly and annotation in Mitofinder

The basic command to run Mitofinder is: 
mitofinder -j [seqid] -1 [left_reads.fastq.gz] -2 [right_reads.fastq.gz] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]   

You can write individual job/scripts to run this command for each of your samples, or you can make a loop and run one job. Choose your own destiny. 

#### If you do individual commands, then you will do this ####

cd < PATH/TO/DIRECTORY/WITH/TRIMMED/READS > 
*Note files should contain '_paired.fastq.gz'

enter the following COMMANDS into the command line:

```
ls *_paired.fastq.gz|awk -F "_" '{print$1}'|sort|uniq|while read line  
do
bash mitofinder_PE.sh $line /Path/to/genbank/reference.gb
done
```

[mitofinder_PE.sh]

mitofinder -j ${1} -1 ${1}_R1_001_paired.fastq.gz -2 ${1}_R2_001_paired.fastq.gz -t trnascan -r $2 -o 5 -p 6 -m 10
*change the -p argument to however many cores/threads you want to use. It will be faster if you use more. 

[mitofinder_SE.sh]

mitofinder -j ${1} -s ${1}_R1_001_unpaired.fastq.gz -t trnascan -r $2 -o 5 -p 6 -m 10

#### To run mitofinder as a job on hydra ####

[mitofinder_annotate_spades.sh](https://github.com/trippster08/genome_skimming_LAB/blob/main/scripts/mitofinder_annotate_spades.sh)
[mitofinder_annotate_spades.job](https://github.com/trippster08/genome_skimming_LAB/blob/main/jobs/mitofinder_annotate_spades.job)

```
sh mitofinder_annotate_spades.sh <path_to_spades_contigs> <genetic_code>
```

## 3b. Extract Gene Orders from Unio Mitogenomes (Optional)
This is more or less taken step 3 of [PhylOcto](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH#3-extract-gene-orders-from-octocoral-mitogenomes)
Dependencies: [Biopython](https://biopython.org), and optionally the webserver at [http://trna.ucsc.edu/tRNAscan-SE/](http://trna.ucsc.edu/tRNAscan-SE/).

#### [gene_order_script.sh](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/gene_order_script.sh)
Recommended: If you'd like to make your life easier, first use [gene_order_script.sh](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/gene_order_script.sh) to aggregate all the .gb and .infos files you need from the MitoFinder results into a new directory.

Example usage:

```
bash gene_order_script.sh dir-with-mitofinder-results
```

#### [gene_order.py](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/gene_order.py)
Use [gene_order.py](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/gene_order.py) by supplying the path to a directory of .gb files you want to extract the gene orders from (also accepts the name of a single .gb file). You can also specify the name of your output file with --name, and whether or not to overwrite pre-existing files with --overwrite.

This will output a csv with a gene order match for every sample as well as information from their .infos file such as mitogenome length, GC content, and circularization. It grabs gene orders from [gene_order_params.py](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/gene_order_params.py) (more on that later).

Example usage:

```
python3 gene_order.py dir-with-gb-files --name output-name --overwrite
```

#### [gene_order_params.py](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/gene_order_params.py)
You can (and probably should) edit [gene_order_params.py](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/gene_order_params.py) to fit your preferences. 
For instance, here are a few valid gene orders.
```
class GeneOrder(Enum):
    A =     ('cox1', 'rrns', 'nad1', 'cob', 'nad6', 'nad3', 'nad4l', 'muts', 'rrnl', 'nad2', 'nad5', 'nad4', 'trnm', 'cox3', 'atp6', 'atp8', 'cox2')
    B =     ('cox1', 'rrns', 'nad1', 'cob', 'cox3', 'trnm', 'nad4', 'nad5', 'nad2', 'rrnl', 'muts', 'nad4l', 'nad3', 'nad6', 'atp6', 'atp8', 'cox2')
    C =     ('cox1', 'rrns', 'nad1', 'cob', 'cox2', 'atp8', 'atp6', 'cox3', 'trnm', 'nad4', 'nad5', 'nad2', 'rrnl', 'muts', 'nad4l', 'nad3', 'nad6')
```

You should also edit the expected_genes variable in this file, which deals with certain genes having multiple names. This variable is a tuple of tuples, where each 2nd order tuple contains every possible name for a gene, where the first name for a gene matches the one displayed in the gene orders. The expected genes should already include all names that mitofinder will use, but if you use a seperate annotating process, yours may differ.

For instance, here is what it should look like:
```
expected_genes = (('cox1',), ('rrns','rns'), ('nad1','nd1'), ('cob','cytb'), ('nad6','nd6'), ('nad3','nd3'), ('nad4l','nd4l'), ('muts',), ('rrnl','rnl'), ('nad2','nd2'), ('nad5','nd5'), ('nad4','nd4'), ('trnm','trna-met'), ('cox3',), ('atp6',), ('atp8',), ('cox2',))
```
So, the gene 'rns' will be forced to be named 'rrns' in order to match the gene order Enums, and 'nd6' becomes 'nad6'.

If you also wish, you can change the start gene in the parameters file. cox1 is widely used as the start for gene orders, but if you wish to have a different start be sure to change the gene order Enums to match.

## 4. Extract ribosomal repeat regions in [Geneious](https://www.geneious.com/)

There may be a better way to do this that I don't know about, but this is what I think most people in IZ are doing right now. 

Basically, you take your clean reads and then manually map those reads to a reference sequence [map_to_reference(s)]. You can download a reference from Genbank, or if you have Sanger Sequencing data for the region of interest that works well too. 

First, upload your PE reads and reference data to a folder in Geneious. 
I think you can only map one sequence at a time (i.e., cannot select multiple PE reads and map to the same reference) otherwise things get jumbled. 
Select you PE reads and the reference seq and click [map_to_reference(s)]. Then you get a list of options to select. I usually keep all the defaults except I change the 'fine tuning' parameter to 'Iterate up to 10 times'. 
Hit 'OK'.
This will take 2-5 minutes per sample, but it you run several at the same time it will make everything slower. 
Then you get an output formatted like {Reads Name} assembled to {Reference Name}.
Click on that file which will show your PE reads and the region you want [Contig View]. You need to manually select the region (bp) from the consensus sequence. 
After you have your region selected, click 'Extract' 
A window will appear for you to name the sequence with the default being 'Consensus extraction'. You will want to name it something intuitive. 
Finally you need to export the sequence(s) to fasta format. 

That is for one gene at a time, but I'm sure you could extract other regions or the flanking sequences if you want at the same time if you know the nucleotide position. 


## 5. Concatenate alignments (if needed/desired) or export alignments in phylip or fasta format

Usually when constructing a phylogeny you will align each gene individually and make individual gene trees, or concatenate alignments and make a gene tree from that. 

If you want to make your life easier, first use alignment_script.sh to aggregate all the .fasta files you need from MitoFinder output into a new directory.

[alignment_script.sh](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/alignment_script.sh)

```
bash alignment_script.sh dir-with-mitofinder-results
```

Next, use [alignment_from_mitofinder.py](https://github.com/Kenneth-Mitchell/PhylOcto_NMNH/blob/main/alignment_from_mitofinder.py) by supplying the path to a directory of .fasta files you want to create alignment files for. It will create a fasta file for each gene found (so, every gene in the mitochondrial genomes) ready for alignment. You can specify whether or not to overwrite pre-existing files with --overwrite.

```
python3 alignment_from_mitofinder.py dir-with-fasta-files --overwrite
```

Once you have your individual genes, you can align them using MAFFT in [AliView](https://ormbunkar.se/aliview/). I would recommend specifying or making sure the refinement method is L-INS-i and saving as a phylip file.


That will take care of the mitogenome genes, but if you want to include the ribosomal genes in the concatenated alignment, you will need to do it manually. There are more than likely better ways to do this, but this works I guess. 

First, you will need to name each of you taxa/samples that you have ribosomal genes for as the same naming scheme as you mitochondrial genes. That way when you concatenate (I do it in Geneious) you concatenate by name so that the indexing matches up. 

After concatenation, you export your alignment as a phylip file. 

## 6. Phylogenetic reconstruction

Once you have your alignment file(s) you can generate phylogenies. [IQ-TREE](http://www.iqtree.org/) is great for maximum-likelihood inference since it is fast, easy to use, and has great documentation. 

The commands are very simple

```
iqtree2 -s alignment.phy -m MFP -B 1000
```

This assumes that the executable should be added to your PATH enviroment variable so that IQ-TREE can be invoked by simply entering iqtree2 in the command-line. 

I won't get into Bayesian inference, but there are good tutorials for [RevBayes](https://revbayes.github.io/tutorials/) available. You can also try [BEAST](https://beast.community/first_tutorial). 


## 7. UCE Phylogenomics implemented in phyluce

I would more or less follow the phyluce [Tutorial I](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#). There are a couple of notes though that are helpful. 

When you get to the [Finding UCE loci](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#finding-uce-loci) step, you will need to modify the regex argument in your [phyluce_assembly_match_contigs_to_probes](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#finding-uce-loci) command.

If using the Unioverse probes, you will need the [Unioverse_probes_reform.fasta](https://github.com/alexfranzen/Unio_Genome_Skimming/blob/main/Unioverse_probes_reform.fasta) file and you will need to specify --regex 'RAPID_GENOMICS_8901_(\d+)'. 

Here's an example command

```
phyluce_assembly_match_contigs_to_probes --contigs /Users/Lampsilis/Documents/Alex_projects/unioverse_test/contigs --probes /Users/Lampsilis/Documents/Alex_projects/unioverse_test/Unioverse_probes_reform.fasta --output /Users/Lampsilis/Documents/Alex_projects/unioverse_test/uce-search-results --regex 'RAPID_GENOMICS_8901_(\d+)'
```

---

The project was funded under a NOAA OER grant to C. McFadden, A. Quattrini and S. Herrera.

Example files are provided for ease of use, these don't represent any actual findings from this project.

Feel free to use code for your own purposes, but please cite as according to the [citation file](CITATION.cff).

Please reach out to me [kmitchell@hmc.edu](mailto:kmitchell@hmc.edu) with anything you need help on.
