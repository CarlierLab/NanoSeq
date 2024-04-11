"""
This module contains a class and methods for read alignments

Copyright 2024 Aurelien Carlier (aurelien.carlier@inrae.fr)
https://github.com/CarlierLab/NanoSeq

This file is part of NanoSeq. NanoSeq is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. NanoSeq is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with NanoSeq.
If not, see <http://www.gnu.org/licenses/>.

"""


import subprocess

class Mapping:
    def __init__(self, barcode, reference, read_location,output_folder):
        self.barcode = barcode
        self.read_location = read_location
        self.output_folder = output_folder
        self.reference = reference

    def minimap(self):
        print(f'mapping {self.barcode} with minimap')
        subprocess.call(f'minimap2 -ax map-ont {self.reference} {self.read_location} > \
                        {self.output_folder}/{self.barcode}_consensus_mapping.sam',shell=True)
        subprocess.call(f'samtools view -Sb {self.output_folder}/{self.barcode}_consensus_mapping.sam > \
                        {self.output_folder}/{self.barcode}_consensus_mapping.bam',shell=True)
        subprocess.call(f'samtools sort {self.output_folder}/{self.barcode}_consensus_mapping.bam -o \
                        {self.output_folder}/{self.barcode}_consensus_mapping.sorted.bam',shell=True)
        subprocess.call(f'samtools index {self.output_folder}/{self.barcode}_consensus_mapping.sorted.bam',shell=True)
        subprocess.call(f'samtools faidx {self.reference} -o {self.output_folder}/{self.barcode}_consensus.fai',shell=True)
        subprocess.call(f'cp {self.reference} {self.output_folder}/{self.barcode}_consensus.fasta',shell=True)
        print("Deleting temporary BAM and SAM files")
        subprocess.call(f'rm {self.output_folder}/{self.barcode}_consensus_mapping.sam  \
                        {self.output_folder}/{self.barcode}_consensus_mapping.bam',shell=True)
        BAM_path = f"{self.output_folder}/{self.barcode}_consensus_mapping.sorted.bam"
        fai_path = f'{self.output_folder}/{self.barcode}_consensus.fai'
        fasta_path = f'{self.output_folder}/{self.barcode}_consensus.fasta'
        BAI_path = f"{self.output_folder}/{self.barcode}_consensus_mapping.sorted.bam.bai"
        return (BAM_path,BAI_path,fai_path,fasta_path)