### This file contains command lines of mapping and genotyping for sample 1A, all other samples can be processed in the same way
 module add bwa/0.7.9a
 module add samtools/0.1.19
 module add gatk/2.5-2
 module add FastqMcf/r534
 ### step 1, filter raw sequencing reads
 /path/to/FastqMcf/fastq-mcf -S -x 5 -l 50  --qual-mean 25 --max-ns 5  -o 1A_1.fq.gz -o 1A_2.fq.gz  all_Illumina_adaptors.fasta  raw.1A_1.fq  raw.1A_2.fq

 ### step 2, bwa mapping
 bwa mem -R '@RG\tID:1A\tSM:bar' -t 8 dm3.fa 1A_1.fq.gz 1A_2.fq.gz | samtools view -bS - | samtools sort - 1A.sort

 ### step 3, genotyping by GATK
 samtools index  1A.sort.bam
 java -jar GenomeAnalysisTK.jar  -R dm3.fa -T UnifiedGenotyper  -I 1A.sort.bam  --output_mode EMIT_ALL_SITES --heterozygosity 0.01 -indelHeterozygosity 0.00125 --annotateNDA --out 1A.genotype.vcf.gz

