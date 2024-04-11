# NanoSeq
NanoSeq is a pipeline used at the LIPME to assemble plasmid and PCR sequences with Nanopore. 
This pipeline is designed to enable routine sequencing of whole plasmids (<50bk) and PCR products (>1kb) on Nanopore Minion, Flongle or Promethion flow cells. It takes a set of basecalled and demultiplexed sequencing reads data in FASTQ or BAM format
The pipeline has been tested on data from R9.4 and R10.4 flow cells and performs best with >100x depth per sample, but lower depth is also possible (as low as 30x with the latest flow cells and chemistry). 

# Motivation
The analysis of nucleotide sequences from plasmids or PCR products is routinely done in many laboratories as part of workflows designed to validate synthetic constructs or analyze variants. This is traditionally done using Sanger sequencing and capillary electrophoresis, often outsourced to third-party providers to avoid large upfront instrument cost and maintenance. Oxford Nanopore Technology (ONT) sequencers are affordable, require low maintenance and infrastructure. Recent progress in ONT raw read accuracy contribute make ONT sequencing a viable alternative for in-lab routine DNA sequencing, but existing workflows lack flexibility, require DNA extracts of high quality and prior knowledge of sample properties (e.g. size). NanoSeq is designed to enable low- to medium scale sequencing of plasmid and PCR using ONT sample preparation kits and MinION or PromethION sequencers. The workflow generates consensus sequences from circular plasmids smaller than ~50kb and PCR products between >1 kb. NanoSeq is robust to partial degradation or contamination of sample nucleic acids, uneven coverage and does not require prior knowledge of size or repeats. 

# Requirements
- Linux PC with 8 cores and >16Gb RAM
- Anaconda or Miniconda

# Installation
## Software dependencies
- Python >3.10 with pandas, numpy, openpyxl, biopython >=1.83
- medaka https://github.com/nanoporetech/medaka
- Canu https://github.com/marbl/canu
- Flye https://github.com/fenderglass/Flye
- FastQC https://github.com/s-andrews/FastQC
- Seqtk https://github.com/lh3/seqtk
- Minimap2 https://github.com/lh3/minimap2
- Filtlong https://github.com/rrwick/Filtlong

The easiest way to install all software dependencies is to create a dedicated conda environment. 

```conda create -n NanoSeq python=3.10```

Activate the environment and install all dependencies

```conda activate NanoSeq```

If Flye is not available system-wide, install with

```conda install -c bioconda -c conda-forge -c defaults Flye```

If Canu is not available system-wide, install with

```conda install -c conda-forge -c bioconda -c defaults canu```

Install medaka with 

```pip install medaka```

Install the rest of the dependencies with

```conda install -c conda-forge -c defaults -c bioconda biopython=1.83 pandas openpyxl seqtk fastqc mash=2.2 filtlong```

Download the NanoSeq scripts

```git clone https://github.com/CarlierLab/NanoSeq.git```





