#!/usr/bin/env python3
"""
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
import argparse
import subprocess
import os
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from assembly import Assembly, PCR_Assembly
from mapping import Mapping

argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input", help="Input folder of subfolders containing demultiplexed ONT reads")
argParser.add_argument("-x", "--xls", help="Sample information as Excel spreadsheet")
argParser.add_argument("-o","--output",help="Output directory")
argParser.add_argument("-m","--model",help="medaka model. Default is r1041_e82_400bps_sup_v4.2.0")
args = argParser.parse_args()
print(args)
folder = os.listdir(args.input)
#Creating output folders

subprocess.call(f"mkdir {args.output}", shell=True)
subprocess.call(f'mkdir {args.output}/failed_to_assemble',shell=True)
subprocess.call(f'mkdir {args.output}/circular_assemblies',shell=True)
subprocess.call(f'mkdir {args.output}/failed_to_circularize',shell=True)
subprocess.call(f'mkdir {args.output}/PCR_assemblies',shell=True)


sample_df = pd.read_excel(args.xls,index_col=None, dtype={'Sample name':str, 'DNA type':str, 'Barcode':str})
#dropping unused sample slots
sample_df = sample_df.dropna(subset=['Barcode','Sample name'])
#converting barcodes and types to lowercase
sample_df["Barcode"] = sample_df.Barcode.apply(lambda x: x.lower())
barcodes = sample_df["Barcode"].values.tolist()

sample_df["DNA type"] = sample_df["DNA type"].apply(lambda x: x.lower())
bad_chars = [';', ':', '!', "*", " ","/",")","("]

# Removing subfolders for which there is no corresponding data in the sample sheet
folder = [i for i in folder if i in barcodes]
print(folder)

for subf in folder:
    #getting corresponding sample info from pandas df, assuming subf corresponds to the assigned barcode
    sub_df = sample_df[sample_df.Barcode == subf]
    sample_id = sub_df.iloc[0]['Sample name']
    sample_type = sub_df.iloc[0]['DNA type'].lower().replace(' ','')
    sample_length = sub_df.iloc[0]['Size (kb)']*1000

    #replacing space characters in sequence name
    sample_id = sample_id.replace(" ","_")
    #checking for illegal characters in sample_id
    for c in bad_chars:
        if c in sample_id:
            sample_id = subf
        else:
            pass

    path = f"{args.input}/{subf}"
    print(f'Assembling {sub_df} with ID {sample_id} and type = {sample_type}')
    try:
        if os.listdir(f"{args.input}/{subf}/final_assemblies") != []:
            print("Already assembled, nothing to do")
            continue
    except FileNotFoundError:
        pass
    if sample_type == "plasmid":
        assembly = Assembly(subf,sample_id,path,path)
        processed_reads = assembly.process_reads()
        QC = assembly.quality_control(processed_reads)
        try:
            length_raw,coverage_raw,total = assembly.getLength(processed_reads)
            print(f'Calculated coverage is {coverage_raw}')
        except IndexError:
            print('No reads. Moving on to next sample')
            continue
        
        if length_raw != 0:
            max_length = int(1.1*length_raw)
        else:
            max_length = 1000000
        #filtering reads with filtlong"  
        filtered_reads = assembly.filter_reads(processed_reads,max_length)
        #adjusting coverage to filtered reads
        length,coverage,total = assembly.getLength(filtered_reads)
        print(f'Calculated coverage after filtering is {coverage}')
        if coverage > 0:
            coverage50 = min(50/coverage,0.99)
            coverage150 = min(150/coverage,0.99)
        else:
            coverage50 = 0.99
            coverage150 = 0.99
            print(f'Unable to calculate coverage. Using all reads')
        print(f'calculated fractions coverage50 = {coverage50} and coverage150 = {coverage150}')

        
        i = 0
        success = False
        contigs = ""
        final_contigs = ""
        assembly_type = ""
        while (success == False) and (i<=10) :
            print(f"iteration number {i}")
            if i<5:
                i +=1
                reads_subset = assembly.subset_reads(filtered_reads,coverage50)
                contigs = assembly.assemble_w_Canu(reads_subset,length)
                if contigs != "":
                    contig_checked = assembly.check_concatemer(contigs,length)
                    circ_contig = assembly.circularize(contig_checked,30,1)
                    if circ_contig != "":
                        circ_seq = SeqRecord(circ_contig,id="%s"%sample_id)
                        subprocess.call(f"mkdir {path}/final_assemblies",shell=True)
                        SeqIO.write(circ_seq,f"{path}/final_assemblies/{sample_id}_circ.fasta","fasta")
                        final_contigs = f"{path}/final_assemblies/{sample_id}_circ.fasta"
                        assembly_type = "circular"
                        success = True
                    else:
                        pass                  

            elif 4<i<10:
                i +=1
                reads_subset = assembly.subset_reads(filtered_reads,coverage50)
                contigs = assembly.assemble_w_Flye(reads_subset)
                if contigs != "":
                    #Checking if contig is circular from assembly_info.txt
                    with open(f'{path}/assemblies/assembly_info.txt','r') as assembly_info:
                        data = assembly_info.read()
                        data = data.split('\n')[1]
                        is_circular = data.split('\t')[3]
                        contig_length = int(data.split('\t')[1])

                    if (is_circular == 'Y') and (contig_length < max_length):
                        print("contig is predicted circular, stitching and rotating")
                        contig_checked = assembly.check_concatemer(contigs,length)
                        circ_seq = SeqRecord(contig_checked,id="%s"%sample_id)
                        #rotating to the left by 100bp to allow polishing of contig 
                        circ_seq = circ_seq[100:] + circ_seq[0:100]
                        subprocess.call(f"mkdir {path}/final_assemblies",shell=True)
                        SeqIO.write(circ_seq,f"{path}/final_assemblies/{sample_id}_circ.fasta","fasta")
                        final_contigs = f"{path}/final_assemblies/{sample_id}_circ.fasta"
                        assembly_type = "circular"
                        success = True
                    else:
                        pass
            elif i == 10:
                i +=1
                print("unable to circularize. Running one last Flye assembly")
                reads_subset = assembly.subset_reads(filtered_reads,coverage50)
                contigs = assembly.assemble_w_Flye(reads_subset)
                if contigs != "":
                    subprocess.call(f"mkdir {path}/final_assemblies",shell=True)
                    subprocess.call(f'mv {contigs} {path}/final_assemblies/{sample_id}_non_circularized.fasta',shell=True)
                    final_contigs = f"{path}/final_assemblies/{sample_id}_non_circularized.fasta"
                    assembly_type = "linear"
                else:
                    print(f"unable to assemble sample {subf}. Moving on.")
                    final_contigs = ""
                    assembly_type = "failed"
            else:
                pass                            
      
    elif sample_type == "pcr":
        assembly = PCR_Assembly(subf,sample_id,path,path)
        processed_reads = assembly.process_reads()
        QC = assembly.quality_control(processed_reads)
        length,coverage,total = assembly.getLength(processed_reads)
        length = sample_length
        coverage = float(total/int(sample_length))
        print(f'Size of fragment from metadata is {length} and calculated coverage is {coverage}')

        if coverage > 0:
            coverage50 = min(50/coverage,0.99)
            coverage150 = min(150/coverage,0.99)
        else:
            coverage50 = 0.99
            coverage150 = 0.99
            print(f'Unable to estimate coverage. Using all reads')

        filtered_reads = assembly.filter_reads(processed_reads,20000)
        reads_subset = assembly.subset_reads(filtered_reads,coverage150)
        contigs = assembly.assemble(reads_subset)
        if contigs != "":
            subprocess.call(f"mkdir {path}/final_assemblies/",shell=True)
            subprocess.call(f'mv {contigs} {path}/final_assemblies/{sample_id}.fasta',shell=True)
            final_contigs = f"{path}/final_assemblies/{sample_id}.fasta"
            assembly_type = "PCR"
    else:
        print("Unknown DNA type: Not processing.")
        pass

    #Polishing final assemblies with Medaka
    if final_contigs != "":
        reads_subset150 = assembly.subset_reads(filtered_reads,coverage150)
        consensus, qualities = assembly.polish(reads_subset150,final_contigs,args.model)
        #Mapping reads back to consensus
        subprocess.call(f'mkdir {args.output}/BAM_files',shell=True)
        mappings = Mapping(subf,consensus,processed_reads, f'{args.output}/BAM_files')
        bam,bai,fai,fasta = mappings.minimap()

        if assembly_type == "circular":
            subprocess.call(f'mkdir {args.output}/circular_assemblies/{subf}',shell=True)
            subprocess.call(f'mv {consensus} {args.output}/circular_assemblies/{subf}',shell=True)
            subprocess.call(f'mv {qualities} {args.output}/circular_assemblies/{subf}',shell=True)
        elif assembly_type == "linear":
            subprocess.call(f'mkdir {args.output}/failed_to_circularize/{subf}',shell=True)
            subprocess.call(f'mv {consensus} {args.output}/failed_to_circularize/{subf}',shell=True)
            subprocess.call(f'mv {qualities} {args.output}/failed_to_circularize/{subf}',shell=True)
        elif assembly_type == "PCR":
            subprocess.call(f'mkdir {args.output}/PCR_assemblies/{subf}',shell=True)
            subprocess.call(f'mv {consensus} {args.output}/PCR_assemblies/{subf}',shell=True)
            subprocess.call(f'mv {qualities} {args.output}/PCR_assemblies/{subf}',shell=True)
    else:
        subprocess.call(f"mv {path} {args.output}/failed_to_assemble",shell=True)
print("All done!")

