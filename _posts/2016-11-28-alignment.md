---
layout: default
title: "Setup, QC and Alignment"
author: "Radhika Khetani"
output: html_document
---

# Setting up Orchestra

To get on Orchestra you want to connect via ssh, and to turn X11 forwarding on. Turning X11 forwarding will let the Orchestra machines open windows on your local machine, which is useful for looking at the data.

> Note for Mac users: you can use the Terminal application/utility to connect to Orchestra as described below. In addition, you will also need to download and install [XQuartz](http://www.xquartz.org/) for X11 forwarding.
>
> &
> 
> Note for Windows users: you will need to download and install [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/) to connect to Orchestra, and [XMing](http://sourceforge.net/project/downloading.php?group_id=156984&filename=Xming-6-9-0-31-setup.exe) for X11 forwarding. [Follow [these instructions](http://www.geo.mtu.edu/geoschem/docs/putty_install.html).]  


Connect to Orchestra using X11 forwarding:

```bash
$ ssh -X your_user_name@orchestra.med.harvard.edu
```

Orchestra is set up with separate login nodes and compute nodes. You don't want to be doing any work on the login node, as that is set aside for doing non computationally intensive tasks and running code on there will make Orchestra slow for everybody. Here is how to connect to a compute node in interactive mode:

```bash
$ bsub -n 2 -Is -q training bash
```

*Note: if you are using your own account use the `interactive` queue for the above command*

Do that and you're ready to roll. You should see that you are now connected to a worker node named by an instrument like `clarinet` or `bassoon`.

Notice that we used the `-n 2` option to allow two cores to be used for the analysis (in general you can set this to larger numbers if required, but we'll leave it at 2 for today so as to avoid overloading the system). 

## Setting up to perform analysis 

The first thing we will do is to checkout a copy of the workshop material from Github into your home directories. 

```bash
$ git clone https://github.com/hms-dbmi/scw.git
$ cd scw/scw2016
```

This repository will remain accessible after the workshop, so you can download the code onto your own machines later. 

Next we will run a small setup script to create the correct environment variables for running the programs we will need on Orchestra 

```bash
$ source setup.sh
```
If you now list the directory contents with the `ls` command, you will see a directory for each of the four tutorials we will be covering this afternoon, containing `R` markdown code that you can use to re-run the analyses later. Today we will interactively go through all the steps outlined in these files.

Since the data files needed for the analyses are fairly large, these are not stored in the repository, and we must copy them over from another directory on Orchestra. 



## Introduction to the data

The raw data we will be using for this part of the workshop lives here `/groups/pklab/scw/scw2015/ES.MEF.data/subset`:

```bash
$ ls /groups/pklab/scw/scw2015/ES.MEF.data/subset
```

***
> [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files are the standard format for sequenced reads, and that is the format you will receive from the sequencing center after they sequence your cDNA libraries.

***

The 4 FASTQ files that we will be working with for this module contain 1000 reads each from ~100 samples. The samples are from a [single-cell RNA-seq experiment](http://genome.cshlp.org/content/21/7/1160.long) where researchers were looking at differences between expression in mouse embryonic fibroblasts (MEF) and embryonic stem (ES) cells from mice. 

***
***
**Note on mutiplexing and demultiplexing**

Libraries prepared from hundreds of single cells can be sequenced in the same lane (multiplexed). During the preparation of these libraries each cell is given a distinct "barcode", a short nucleotide sequence that will enable us to separate them out after sequencing. Illumina provides adaptors with a barcodes, and software that can demultiplex the data for you. So, if you use these Illumina barcodes, the sequencing center will return demultiplexed fastq files to you. 

However, Illumina offers only ~96 distinct barcode combinations (as of March 2015). For single cell work where we are interested in simultaneously sequencing more than 96 cells such as with many of the more recent droplet-based microfluidics approaches, we need additional barcodes. To this end, many groups design their own sets of barcodes; since Illumina's software is unable to use these to separate the samples, you will have to perform demultiplexing after receiving the data from the sequencing center. 

This is outside the scope of this workshop, but it is important to note that this will add an additional step prior to the three steps listed below.

***
***

We'll be taking this small subset of reads and performing the following steps:

1. looking at them to make sure they are of good quality 
* aligning them to the mouse genome 
* producing a table of number of reads aligning to each gene for each sample

The counts table generated from step 3 will be the starting point for the more interesting downstream functional analyses. We'll use these subsetted ES and MEF files to demonstrate the workflow; then, we'll look at pre-computed counts results on a full set of samples for the functional analyses.

## Copy across data into your own directories

We will now copy over the test data over into your alignment directory.

```bash
$ cd tutorials/alignment
$ cp -r /groups/pklab/scw/scw2015/ES.MEF.data/subset .
```

These commands mean:

* change directories to your home directory (`cd` without anything following it will always bring you to your home directory)
* change into the directory (`cd`) named 'alignment' 
* copy (`cp`) the folder `/groups/pklab/scw/scw2015/ES.MEF.data/subset` and everything underneath it (using the `-r`) to the current directory (denoted by a period `.`)


# Quality control

With the start of any data analysis it is important to poke around at your data to see where the warts are. We are expecting single-cell datasets to be extra messy; we should be expecting failures in preparing the libraries, failures in adequately capturing the transcriptome of the cell, libraries that sequence the same reads repeatedly and other issues. In particular we expect the libraries to not be very complex; these are just tags of the end of a transcript and so for each transcript there are a limited number of positions where a read can be placed. We'll see how this feature skews some of the more commonly used quality metrics now.

For RNA-Seq data many common issues can be detected right off the bat just by looking at some features of the raw reads. The most commonly used program to look at the raw reads is [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/). To run FastQC on the cluster we have to load the necessary module:

```bash
$ module load seq/fastqc/0.11.3
```
FastQC is pretty fast, especially on small files, so we can run FastQC on one of the full files (instead of just on the subset). First lets copy one of those files over:

```bash
$ mkdir ~/scw/scw2016/tutorials/alignment/fastq
$ cp /groups/pklab/scw/scw2015/ES.MEF.data/fastq/L139_ESC_1.fq \
	~/scw/scw2016/tutorials/alignment/fastq/
$ cd ~/scw/scw2016/tutorials/alignment/fastq
```

Now we can run FastQC on the file by typing:

```bash
$ fastqc L139_ESC_1.fq
```
And look at the nice HTML report it generates with Firefox:

```bash
$ firefox L139_ESC_1_fastqc.html
```    
Let's analyse some of the plots: 

* The **per base sequence quality** plot shows some major quality problems during sequencing; having degrading quality as you sequence further is normal, but this is a severe drop off. Severe quality drop offs like this are generally due to technical issues with the sequencer, it is possible it ran out or was running low on a reagent. The good news is it doesn't affect all of the reads, the median value still has a PHRED score > 20 (so 1 in 100 probability of an error), and most aligners can take into account the poor quality so this isn't as bad as it looks.
* More worrying is the non-uniform **per base sequence content** plot. It depends on the genome, but for the mouse if you are randomly sampling from the transcriptome then you should expect there to be a pretty even distribution of GC-AT in each sequence with a slight GC/AT bias. We can see that is not the case at all 
* In the next plot, the **per sequence GC content** plot, has a huge GC spike in the middle. Usually you see plots like these when the overall complexity of the reads that are sequenced is low; by that we mean you have tended to sequence the same sample repeatedly, and have very little diversity.
* That notion is reinforced looking at the **sequence duplication levels** plot. If we de-duplicate the reads, meaning remove reads where we have seen the same exact read twice, we'd throw out > 75% of the data. 
* It is also reinforced by the list of kmers that are more enriched than would be expected by chance; a kmer is just every possible k length mer that is seen in a sequence. We'd expect all of these features if we were sequencing the same sequences repeatedly. One thing we would not expect, however, is the big uptick at the end of the **kmer content** plot; the sequences at the end look like some kind of adapter contamination issue. We'd expect these reads to not align unless we trimmed the adapter sequence off.

What are those sequences? We can search for the reads that have one of those enriched sequences with `grep` (*g*lobally search a *r*egular *e*xpression and *p*rint) which print out every line in a file that matches a search string. grep the L139_ESC_1.fq file like this:

```bash
$ grep ACTTGAA L139_ESC_1.fq
```    
You should see a lot of sequences that look like this:
`TTGCTAGATATCAGTTGCCTTCTTTTGTTCAGCTAAGGAACTTGAA`

If we BLAST this sequence to the mouse genome, we come up empty, so it is some kind of contaminant sequence, it isn't clear where it comes from. The protocol listed here doesn't have too many clues either. If we could figure out what these sequences are, it would help troubleshoot the preparation protocol and we might be able to align more reads. As it is, these sequences are unlikely to align to the mouse genome, so they mostly represent wasted sequencing.

***
***
**FastQC reports and the implications for data quality (more information)**

* [This blog post](http://bioinfo-core.org/index.php/9th_Discussion-28_October_2010) has very good information on what bad plots look like and what they mean for your data.
* [This page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) has an example report from a good dataset.
* Please read [this note on evaluating fastqc results](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/2%20Basic%20Operations/2.2%20Evaluating%20Results.html), before being too alarmed by the red "X"s or the orange "!"s in the FastQC report.

***
***

# Alignment

For aligning RNA-seq reads it is necessary to use an aligner that is splice-aware; reads crossing splice junctions have gaps when aligned to the genome and the aligner has to be able to handle that possibility. There are a wide variety of aligners to choose from that handle spliced reads, e.g. STAR, HISAT2. STAR requires a lot of memory whereas HISAT2 is fast and does not require as much memory.

For this exercise we will use HISAT2, so let's load the 2 modules we will need. 
	
```bash
$ module load seq/hisat2/2.0.4
$ module load seq/samtools/1.3
```

To align reads with HISAT2 you need three things.

1. The genome sequence of the organism you are working with in FASTA format. You will use this to make an index.
1. The known splice junctions that can be obtained from the gene annotation or GTF file. You can use the GTF file to make a text file fo splice junctions.
1. The FASTQ file of reads you want to align.

First you must make an index of the genome sequence; this allows the HISAT2 algorithm to rapidly find regions of the genome where each read might align. In addition, if you are planning to use the known gene annotations as a guide for alignment, you will need to the splice junction text file. Both of these have already been created for use with this workshop, so **don't type in the following commands**:

```bash
## DO NOT RUN THIS
$ hisat2-build -p <num_cores> <genome.fa> <prefix_for_index>
	
$ hisat2_extract_splice_sites.py <genes.gtf> > <splicesites.txt>
```
We will be using precomputed indices for the mm10 genome today, and the path to those are listed below:

```bash
/n/scratch2/scw2016/reference_index/mm10_hisat2/
/n/scratch2/scw2016/reference_index/mm10_hisat2/splicesites.txt
```

Let's start by aligning the two ESC samples and two MEF samples we copied over earlier:

```bash
$ cd ~/scw/scw2016/tutorials/alignment/subset

$ bsub -J L139_ESC_1 -W 00:20 -n 3 -o %J.L139_ESC_1.out -e %J.L139_ESC_1.err -q training "hisat2 -p 2 \
-x /n/scratch2/scw2016/reference_index/mm10_hisat2/mm10_hisat -U L139_ESC_1.subset.fastq \
--known-splicesite-infile /n/scratch2/scw2016/reference_index/mm10_hisat2/splicesites.txt \
| samtools view -Sbo L139_ESC_1.subset.bam -"

$ bsub -J L139_ESC_2 -W 00:20 -n 3 -o %J.L139_ESC_2.out -e %J.L139_ESC_2.err -q training "hisat2 -p 2 \
-x /n/scratch2/scw2016/reference_index/mm10_hisat2/mm10_hisat -U L139_ESC_2.subset.fastq \
--known-splicesite-infile /n/scratch2/scw2016/reference_index/mm10_hisat2/splicesites.txt \
| samtools view -Sbo L139_ESC_2.subset.bam -"

$ bsub -J L139_MEF_50 -W 00:20 -n 3 -o %J.L139_MEF_50.out -e %J.L139_MEF_50.err -q training "hisat2 -p 2 \
-x /n/scratch2/scw2016/reference_index/mm10_hisat2/mm10_hisat -U L139_MEF_50.subset.fastq \
--known-splicesite-infile /n/scratch2/scw2016/reference_index/mm10_hisat2/splicesites.txt \
| samtools view -Sbo L139_MEF_50.subset.bam -"

$ bsub -J L139_MEF_49 -W 00:20 -n 3 -o %J.L139_MEF_49.out -e %J.L139_MEF_49.err -q training "hisat2 -p 2 \
-x /n/scratch2/scw2016/reference_index/mm10_hisat2/mm10_hisat -U L139_MEF_49.subset.fastq \
--known-splicesite-infile /n/scratch2/scw2016/reference_index/mm10_hisat2/splicesites.txt \
| samtools view -Sbo L139_MEF_49.subset.bam -"
```

Each of these should complete in about 2 minutes. Since we ran them all in parallel on the cluster, the whole set should take about 2 to 4 minutes instead of 10 - 15. Full samples would take longer.

***
`-J` names the job so you can see what it is when you run bjobs to check the status of the jobs. 

`-W 00:20` tells the scheduler the job should take about 20 minutes. 

`-q training` submits the job to the training queue. [Note: If you are using your own account, please use `-q short` or `-q priority` to submit the jobs to the short queue; the training queue can only be used by people using the training accounts.]

`-n 3` requests 3 cores.

`-o` & `-e` specify where to save the information that is output by the program, standard output and standard error respectively. 

The syntax of the HISAT2 command is:

```bash
## DO NOT RUN THIS

hisat2 -p <num_threads> -x <path_to_genome_index> \
-U <input_fastq_unpaired> --known-splicesite-infile <splicesites.txt> \
| samtools view -Sbo <name_of_ouput.bam> -
```

*Note: HISAT2 outputs alignment results in SAM format as standard output, so we pipe it to `samtools` and convert it to BAM in the same command.*

***
This method of submitting one job at a time is fine for a small number of samples, but if you wanted to run a full set of hundreds of cells, doing this manually for every sample is a waste of time and prone to errors. You can run all of these automatically by writing a `for` loop:

```bash
# don't type this in, it is here for future reference

module load seq/hisat2/2.0.4
module load seq/samtools/1.3

for file in *.fastq
do
    samplename=$(basename $file .fastq)
    
    bsub -J $samplename -W 00:20 -n 8 -o %J.$samplename.out \
    -e %J.$samplename.err -q short "hisat2 -p 7 -x genome_index \
    -U $samplename splicesites.txt | samtools view -Sbo $samplename.bam -"
    
done
```

This will loop over all of the files with a `.fastq` extension in the current directory and align them with HISAT2. We'll skip ahead now to doing some quality control of the alignments and finally counting the reads mapping to each feature.

## Quality checking the alignments

There are several tools to spot check the alignments, it is common to run [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc) on the alignment files to generate some alignment stats and to determine the overall coverage, how many genes were detected and so on. Another option for a suite of quality checking tools is RSeQC. For real experiments it is worth it to look at the output of these tools and see if anything seems amiss.

# Counting reads with featureCounts

The last step is to count the number of reads mapping to the features we are interested in. Quantification can be done at multiple levels; from the level of counting the number of reads supporting a specific splicing event, to the number of reads aligning to an isoform of a gene or the total reads mapping to a known gene. We'll be quantifying the latter, i.e. the total number of reads that can uniquely be assigned to a known gene; basically looking at the location of read alignment on the genome and putting it together with the location of the gene on the genome (this information is contained in the [GTF](http://mblab.wustl.edu/GTF2.html)/annotation file). There are several tools to do this, we will use [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) because it is very fast and accurate.

```bash
$ module load seq/subread/1.4.6-p3    #featureCounts is part of the subread package
    
$ featureCounts --primary -a /n/scratch2/scw2016/mm10_UCSC_genes.gtf \
	-o combined.featureCounts *.bam
```    
    
***
`--primary` tells featureCounts to only count the primary alignment for reads that map multiple times to the genome. This ensures we do not double count reads that map to multiple places in the genome.
***

We need to massage the format of this file so we can use it. We'll take the first column, which is the gene ids and every column after the 6th, which has the counts of each sample.

```bash
$ sed 1d combined.featureCounts | cut -f1,7- | sed s/Geneid/id/ > combined.counts
```

This command means *s*team *ed*it (`sed`) the file `combined.featureCounts` by 
***
`1d` deleting the first line

`-f1,7-` keeping the first column and everything from the 7th column on 

`s/Geneid/id/` changing the phrase "Geneid" to "id". 
***
This outputs a table with "I" rows of genes and the "J" columns of samples. Each entry in the table is the number of reads that can be uniquely assigned to the gene "i" in the sample "j". This file is of now ready and in the correct format for loading into R for differential expression analysis.

## Contact info:
email: [rkhetani@hsph.harvard.edu](mailto: rkhetani@hsph.harvard.edu)

webpage: [http://bioinformatics.sph.harvard.edu/](http://bioinformatics.sph.harvard.edu/)

twitter: [rs_khetani](https://twitter.com/rs_khetani)
