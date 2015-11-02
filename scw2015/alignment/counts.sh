#!/bin/bash/

set -x 

#counts
	module load seq/subread/1.4.6-p3			#featureCounts is part of the subread package
	    
	featureCounts --primary -a /groups/shared_databases/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o combined.featureCounts L139_ESC_1-tophat/L139_ESC_1.bam L139_ESC_2-tophat/L139_ESC_2.bam L139_MEF_49-tophat/L139_MEF_49.bam L139_MEF_50-tophat/L139_MEF_50.bam
	sed 1d combined.featureCounts | cut -f1,7- | sed s/Geneid/id/ > combined.counts

