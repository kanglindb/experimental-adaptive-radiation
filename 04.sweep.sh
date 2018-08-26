### This file contains command lines of sweep detection for selection regime CS of generation 65, all other samples can be processed in the same way
  module add samtools/0.1.19
  module add popoolation/1.2.2
  module add epd/7.3-1
  module add poolhmm/1.4.4
	  
### step 1, get and filter pileup file  
  samtools mpileup -Q 20 CS.g65.sort.bam > CS.g65.pileup
  perl /path/to/popoolation/identify-genomic-indel-regions.pl --input CS.g65.pileup --output CS.g65.indel.gtf --min-count 2 
  perl /path/to/popoolation/current/basic-pipeline/filter-pileup-by-gtf.pl  --input CS.g65.pileup  --gtf CS.g65.indel.gtf  --output CS.g65.filter.pileup 
       
### step 2, split pileup file by chromosome and run poolhmm
  perl splitPileupByChr.pl -in CS.g65.pileup.gz -l chr.list -out CS.g65
  for chr in chr2L chr2R chr3L chr3R chr4 chrX
  do
    python /path/to/poolhmm/bin/pool-hmm.py  -n 88 -c 5 -C 800 -q 20 -e sanger -p -k 0.0000000001 --theta 0.00335584 -S -P 12 -f CS.g65.$chr
  done

