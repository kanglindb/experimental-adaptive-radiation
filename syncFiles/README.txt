This folder contains information of PSIs and rig-PSIs in modified sync format with one additional 4th column (details about the sync 
file format please refer to https://code.google.com/archive/p/popoolation2/wikis/Tutorial.wiki) described as following:

Columns:
1        : chromosome
2        : position
3        : reference
4        : haplotype block ID; ID with "rHB" indicates haplotype block with 2 or more rig-PSIs; "-" indicates non-haplotype-block
5-(N+4)  : allele frequencies of samples 1 to N


Files were given by each selection rigime (CS, DS, HS, KS, and SS). Each file contains information of 15 samples for a given selection 
regime, the first 5 samples are the 5 replicates from the initial populaiton (UC4), the next 5 samples are the 5 replicates from 
generation 45 of that selection regime (e.g. SS45 in file "PSI.SS.sync"), and the last 5 samples are the 5 replicates from generation 65 
of that selection regime (e.g. SS65 in file "PSI.SS.sync").
