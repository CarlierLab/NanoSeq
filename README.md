# NanoSeq
NanoSeq is a pipeline used at the LIPME to assemble plasmid and PCR sequences with Nanopore. 
This pipeline is designed to enable routine sequencing of whole plasmids (<50bk) and PCR products (>1kb) on Nanopore Minion, Flongle or Promethion flow cells. It takes a set of basecalled and demultiplexed sequencing reads data in FASTQ or BAM format
The pipeline has been tested on data from R9.4 and R10.4 flow cells and performs best with >100x depth per sample, but lower depth is also possible (as low as 30x with the latest flow cells and chemistry). 
