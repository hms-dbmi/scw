#!/bin/bash/

set -x 

#setup
	cd
	mkdir -p workshop
	cd workshop
#	cp -r /groups/pklab/scw/scw2015/ES.MEF.data/subset .

#alignment
	module load seq/tophat/2.1.0
	cd ~/workshop/subset
	    
	bsub -J L139_ESC_1 -W 00:20 -n 2 -q training "tophat -p 2 -o L139_ESC_1-tophat --no-coverage-search --transcriptome-index=/groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/tophat2_trans/genes /groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome L139_ESC_1.subset.fastq; mv L139_ESC_1-tophat/accepted_hits.bam L139_ESC_1-tophat/L139_ESC_1.bam"
	    
	bsub -J L139_ESC_2 -W 00:20 -n 2 -q training "tophat -p 2 -o L139_ESC_2-tophat --no-coverage-search --transcriptome-index=/groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/tophat2_trans/genes /groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome L139_ESC_2.subset.fastq; mv L139_ESC_2-tophat/accepted_hits.bam L139_ESC_2-tophat/L139_ESC_2.bam"
	    
	bsub -J L139_MEF_49 -W 00:20 -n 2 -q training "tophat -p 2 -o L139_MEF_49-tophat --no-coverage-search --transcriptome-index=/groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/tophat2_trans/genes /groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome L139_MEF_49.subset.fastq; mv L139_MEF_49-tophat/accepted_hits.bam L139_MEF_49-tophat/L139_MEF_49.bam"
	    
	bsub -J L139_MEF_50 -W 00:20 -n 2 -q training "tophat -p 2 -o L139_MEF_50-tophat --no-coverage-search --transcriptome-index=/groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/tophat2_trans/genes /groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome L139_MEF_50.subset.fastq; mv L139_MEF_50-tophat/accepted_hits.bam L139_MEF_50-tophat/L139_MEF_50.bam"
