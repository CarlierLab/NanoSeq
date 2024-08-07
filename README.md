# NanoSeq
NanoSeq is a pipeline used at the LIPME to assemble plasmid and PCR sequences with Nanopore. 
This pipeline is designed to enable routine sequencing of whole plasmids (<50kb) and PCR products (>1kb) on Nanopore Minion, Flongle or Promethion flow cells. It takes a set of basecalled and demultiplexed sequencing reads data in FASTQ or BAM format
The pipeline has been tested on data from R9.4 and R10.4 flow cells and performs best with >100x depth per sample, but lower depth is also possible (as low as 30x with the latest flow cells and chemistry). 

# Motivation
The analysis of nucleotide sequences from plasmids or PCR products is routinely done in many laboratories as part of workflows designed to validate synthetic constructs or analyze variants. This is traditionally done using Sanger sequencing and capillary electrophoresis, often outsourced to third-party providers to avoid large upfront instrument cost and maintenance. Oxford Nanopore Technology (ONT) sequencers are affordable, require low maintenance and infrastructure. Recent progress in ONT raw read accuracy contribute make ONT sequencing a viable alternative for in-lab routine DNA sequencing, but existing workflows lack flexibility, require DNA extracts of high quality and prior knowledge of sample properties (e.g. size). NanoSeq is designed to enable low- to medium scale sequencing of plasmid and PCR using ONT sample preparation kits and MinION or PromethION sequencers. The workflow generates consensus sequences from circular plasmids smaller than ~50kb and PCR products >0.5 kb. NanoSeq is robust to partial degradation or contamination of sample nucleic acids, uneven coverage and does not require prior knowledge of size or repeats. 

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
- Miniasm and Minipolish https://github.com/rrwick/Minipolish for PCR assembly
- any2fasta to convert miniasm GFA output to fasta
- Mash: https://github.com/marbl/Mash

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

```conda install -c conda-forge -c defaults -c bioconda biopython=1.83 pandas openpyxl seqtk fastqc mash=2.2 filtlong any2fasta minipolish miniasm```

Finally, download the NanoSeq scripts

```git clone https://github.com/CarlierLab/NanoSeq.git```

# Usage
```
usage: NanoSeq [-h] [-i INPUT] [-x XLS] [-o OUTPUT] [-m MODEL]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input folder of subfolders containing demultiplexed ONT reads
  -x XLS, --xls XLS     Sample information as Excel spreadsheet
  -o OUTPUT, --output OUTPUT
                        Output directory
  -m MODEL, --model MODEL
                        medaka model. Default is r1041_e82_400bps_sup_v4.3.0
```

Required arguments:
- An Excel (.xlsx) spreadsheet containing sample information (see example). Mandatory fields include "Sample name", "DNA type", "Barcode". "Size (kb)" is a mandatory field but can be left empty.
- A directory containing sequencing reads in FASTQ format, structured according to barcode:
```
barcode_pass/
├── barcode01
│   └── AQP668_pass_barcode01_6bf36ca6_2d52d78f_0.fastq.gz
├── barcode02
│   ├── AQP668_pass_barcode02_6bf36ca6_2d52d78f_0.fastq.gz
│   └── AQP668_pass_barcode02_6bf36ca6_2d52d78f_1.fastq.gz
├── barcode03
│   ├── AQP668_pass_barcode03_6bf36ca6_2d52d78f_0.fastq.gz
│   └── AQP668_pass_barcode03_6bf36ca6_2d52d78f_1.fastq.gz
├── barcode04
│   └── AQP668_pass_barcode04_6bf36ca6_2d52d78f_0.fastq.gz
├── barcode61
│   └── AQP668_pass_barcode61_6bf36ca6_2d52d78f_0.fastq.gz
└── barcode72
    └── AQP668_pass_barcode72_6bf36ca6_2d52d78f_0.fastq.gz
```

# Example data
Sample data is provided as a directory containing .fastq.gz files and a metadata Excel spreadsheet. Launch the analysis on the example data with:

```python NanoSeq -i NanoSeq/example_data/data -x NanoSeq/example_data/sample_sheet_example.xlsx -o assemblies -m r1041_e82_400bps_sup_v4.3.0 ```
 
# Method

Assembly of ONT sequencing data with Flye or Canu (and many other assemblers) is relatively straightforward for large molecules or genomes. However, small, circular plasmids cause particular issues because because assemblers may artificially extend contig ends. Solutions exist to tackle this issue in the context of genome assembly (e.g. Trycycler https://github.com/rrwick/Trycycler), but they are cumbersome to use for a large number of samples. Other solutions (e.g. OnRAMP https://onramp.boylelab.org/ or the EPI2ME plasmid validation pipeline) require either reference FASTA files or to have prior knowledge of plasmid size. In a typical lab, this is not always known, or desirable if one wants to mix samples from plasmid or linear PCR products in the same library. 

Our solution to small plasmid assembly is to run the assembly using different subsets until circularization is detected. The size of the assembly is checked at each iteration against the calculated plasmid size. First, reads are assembles with Canu (with read correction). Only one contig is selected for analysis. If the sample is not heavily contaminated, this is often not a problem, otherwise the first contig in the output of the assembler is picked.

The *check_concatemer* function detects contigs that are more than 1.9x longer than the estimated molecule size (from read lengths distribution). Contigs that exceed this threshold are sliced in half and each half is compared to the other using Mash. If the 2 halves are >98% identical, the contig is sliced in half (leaving 20% overhangs to allow for circularization). 
The *circularize* function attempts to find overlaps of 30bp (allowing 1 mismatch) between the 5’end of the contig and an internal fragment. Upon detection of a match, the downstream sequence (until the contig end) is aligned to a sequence of identical length starting at the beginning of the contig. If both sequences are >95% identical, the contig is considered circular. The contig is considered non-circular if no overlap is detected, if multiple seed matches are detected, or if the full-length overlaps are <95% identical. Sequences are then rotated by 100 bp to improve correction of contig ends during polishing. Upon failure at the circularization step, assembly is repeated with Canu up to 4 times. Beyond that, the assembly step is repeated again with Flye until circularization is attained or the number of attempts exceed 9 times in total. If circularization still fails after 9 attempts, the contig is considered “non circular” and output as such.  
Final contigs are polished with Medaka. Consensus sequences are generated in FASTA and FASTQ formats, with positions with consensus qualities < Q10 masked. 

PCR (linear) products are identified as such in the metadata file and treated separately. Read data is subset to approximately 50x depth and assembled directly with Miniasm + Minipolish, and assemblies polished with Medaka to generate high quality consensus sequences. 

# Output
The output of the pipeline is written to the directory specified by the -o parameter. It is organized in 6 sub-directories as below:

```
output_dir
├── BAM_files
│   ├── barcodeXX_consensus.fai
│   ├── barcodeXX_consensus.fasta
│   ├── barcodeXX_consensus_mapping.sorted.bam
│   ├── barcodeXX_consensus_mapping.sorted.bam.bai
│   ├── barcodeXY_consensus.fai
│   ├── barcodeXY_consensus.fasta
│   ├── barcodeXY_consensus_mapping.sorted.bam
│   └── barcodeXY_consensus_mapping.sorted.bam.bai
├── PCR_assemblies
│   └── barcodeXX
│       ├── {sample_name}__consensus.fasta
│       └── {sample_name}_qualities.fastq
├── circular_assemblies
│   └── barcodeXY
│       ├── {sample_name}_consensus.fasta
│       └── {sample_name}_qualities.fastq
├── failed_to_assemble
├── failed_to_circularize
└── read_qualities
    ├── barcodeXX_reads_fastqc.html
    └── barcodeXY_reads_fastqc.html
```

- The output of FastQC in HTML format for each sample is stored in the `read_qualities` folder
- Plasmid consensus sequences that completed without failure are written to the `circular_assembly` folder.
- Plasmid consensus sequences that could not be circularized because unique sequence overlaps could not be detected at contig ends, or because assemblies contained more than one contig are written to `failed_to_circularize`
- PCR (linear) fragment assemblies are written to `PCR_assemblies`
- The `BAM_files` folder contains the alignments of all reads for each sample to the consensus sequence (.bam), as well as alignment indices (.bai), consensus sequence in fasta format and index in .fai format. Alignments can be viewed with e.g. the IGV browser.
- The `failed to assemble`folder is for troubleshooting only. It contains the sample data that failed to generate assemblies.

# Citation 
If you found this software useful, please cite:
Carlier, A., & Moreau, S. (2024). CarlierLab/NanoSeq: NanoSeq v0.1.3 (v0.1.3). Zenodo. https://doi.org/10.5281/zenodo.11366140
