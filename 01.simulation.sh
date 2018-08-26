### step 1, format the input file for mimicree using DGRP haplotype which can be found at http://dgrp2.gnets.ncsu.edu/
  perl formatHaplotypeFromDGRP.pl -i dgrp2.genotypes -n 108 -m 1967223  -o dgrp2.hap4sim.2M

### step 2, do simulation using mimicree
  mkdir sim.2M
  mkdir cmh.2M
  ### chromosome can't contain chr...
  sed 's/chr//g' dgrp2.hap4sim.2M > dgrp2.hap4sim.2M.2
  for i in {1..10}
  do
    java -Xmx10g -jar  /path/to/MimicrEE/MimicrEESummary.jar --seed $i --haplotypes-g0 ./dgrp2.hap4sim.2M.2  --recombination-rate ./dmel.rr.txt --output-mode 41 --replicate-runs 5 --output-format sync --threads 1  --output-file ./sim.2M/s01-n$i.g41.sync
    java -Xmx10g -jar  /path/to/MimicrEE/MimicrEESummary.jar --seed $i --haplotypes-g0 ./dgrp2.hap4sim.2M.2  --recombination-rate ./dmel.rr.txt --output-mode 61 --replicate-runs 5 --output-format sync --threads 1  --output-file ./sim.2M/s01-n$i.g61.sync
    ### merge Sync files from generations 41 and 61
    perl mergeFile.pl -k 1,2,3 sim.2M/s01-n$i.g41.sync sim.2M/s01-n$i.g61.sync > sim.2M/s01-n$i.sync
    ### CMH test by using script from popoolaitons (https://sourceforge.net/p/popoolation2)
    perl /path/to/popoolation2/cmh-test.pl --min-count 1 --min-coverage 1 --max-coverage 100000  --population 1-2,3-4,5-6,7-8,9-10,11-12,13-14,15-16,17-18,19-20 --input sim.2M/s01-n$i.sync --output cmh.2M/s01-n$i.cmh
  done

### step 3, cat CMH results and determine the p-value for each of the five categories (<10%, 10-20%, 20-30%, 30-40%, and 40-50%)
  cat cmh.2M/s01-n*.cmh > All.cmh
  perl countPvalueFromSim.pl -i All.cmh -fdr 1 > pvalue.thresholds.fdr001.output

