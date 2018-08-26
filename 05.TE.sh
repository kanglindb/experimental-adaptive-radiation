### This file contains command lines of (transposoble element) TE identification for selection regime CS of generation 65
  module add "bio/bwa/0.7.9a"
  module add "bio/samtools/0.1.19"
  module add "bio/popoolationte/1.02" 

### step 1, map pair-end reads separately to genome reference with known TE sequences
  bwa mem -R '@RG\tID:10A\tSM:bar' -t 8 dm3.te.combined.fa  CS.g65_1.fq.gz > CS.g65.1.sam
  bwa mem -R '@RG\tID:10A\tSM:bar' -t 8 dm3.te.combined.fa  CS.g65_2.fq.gz > CS.g65.2.sam

### step 2, run popoolation TE
  perl /path/to/popoolationte/samro.pl --sam1 CS.g65.1.sam --sam2 CS.g65.2.sam --fq1 CS.g65_1.fq.gz --fq2 CS.g65_2.fq.gz  --output CS.g65.pe.sam
  rm CS.g65.1.sam CS.g65.2.sam
  samtools view -Sb CS.g65.pe.sam | samtools sort - CS.g65.pe.sorted
  samtools view CS.g65.pe.sorted > CS.g65.pe.sorted.sam
  rm CS.g65.pe.sam
  perl /path/to/popoolationte/identify-te-insertsites.pl   --input CS.g65.pe.sorted.sam --te-hierarchy-file te.hierarchy.txt  --te-hierarchy-level family --narrow-range 99 --min-count 3 --min-map-qual 15 --output CS.g65.te-fwd-rev.txt 
  perl /path/to/popoolationte/genomic-N-2gtf.pl  --input dm3.te.combined.fa > CS.g65.poly_n.gtf  
  perl /path/to/popoolationte/crosslink-te-sites.pl  --directional-insertions CS.g65.te-fwd-rev.txt --min-dist 99 --max-dist 500 --output CS.g65.te-inserts.txt --single-site-shift 100 --poly-n CS.g65.poly_n.gtf --te-hierarchy te.hierarchy.txt  --te-hier-level order

