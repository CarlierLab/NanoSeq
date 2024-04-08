# written by Aurélien Carlier 2023
import argparse
import subprocess
import random
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np

class Assembly:
    def __init__(self, barcode, name, read_location,output_folder):
        self.barcode = barcode
        self.name = name
        self.read_location = read_location
        self.output_folder = output_folder

    def process_reads(self):
        print("formatting read files")
        subprocess.call(f'zcat {self.read_location}/*.gz > {self.read_location}/{self.barcode}_reads.fastq', shell=True)
        path_to_reads = f'{self.read_location}/{self.barcode}_reads.fastq'
        return path_to_reads
    
    def quality_control(self,reads):
        """reads is the path to input read file"""
        print("performing QC with FastQC")
        self.reads = reads
        subprocess.call(f'mkdir {self.output_folder}/read_quality',shell=True)
        subprocess.call(f'fastqc {reads} -o {self.output_folder}/read_quality', shell=True)
        return True

    def getLength(self,reads_file):
        print("getting length of amplicon or plasmid")
        self.reads_file = reads_file
        reads = SeqIO.parse(reads_file,"fastq")
        out_size = []

        # getting length of each read
        read_len = np.array([len(read.seq) for read in reads])
        max_len = read_len.max()
        # creating the array with bin boundaries of width 100
        bin_array = np.arange(1,(max_len + 1),100)
        # returning the index in bin_array for each value in read_len
        bins = np.digitize(read_len,bin_array)
        # getting the number of occurences for each bin
        count_bins = [np.count_nonzero(bins == i) for i,j in enumerate(bin_array)]
        print(f'count bins {count_bins}')
        print(f'count_bins length = {len(count_bins)}')
        # some very long spurious reads may skew the distribution. Removing outliers
        # we do this by finding the index of the largest size bin containing more than 5 reads
        try:
            slice_index = [i for i,j in enumerate(count_bins) if j > 4][-1]
        except IndexError:
            slice_index = count_bins[-1]
        print(f'slice index {slice_index}')
        # then slicing the list
        count_bins = count_bins[:(slice_index + 1)]
        print(f'count_bins length after slicing = {len(count_bins)}')
        # considering reads in largest 2/3 of bin array to find index of local maximum
        offset = len(count_bins)//3
        sliced_count_bins = count_bins[len(count_bins)//3:]
        local_max = max(sliced_count_bins)
    
        print(f'local max {local_max}')
        if local_max > 4:
                # considering maximum valid if bin contains a minimum of 5 reads
                index_max = sliced_count_bins.index(local_max) + offset
                # getting the corresponding value in bin_array
                plsmd_size = int(bin_array[index_max])
        else:
            # assuming the last bin corresponds to the plasmid size
            index_max = len(count_bins)
            # getting the corresponding value in bin_array
            plsmd_size = int(bin_array[index_max])

        if plsmd_size > 2000:
            est_cov = np.sum(read_len)/plsmd_size
            out_size = [plsmd_size,est_cov,np.sum(read_len)]
        else:
            #approximating size and coverage
            plsmd_size = max_len
            est_cov = np.sum(read_len)/plsmd_size
            out_size = [plsmd_size,est_cov,np.sum(read_len)]
        print(f"Calculated size is {out_size[0]} and coverage is {int(out_size[1])}x")
        return out_size
    
    def filter_reads(self,processed_reads,max_length):
        min_length = 500
        path_to_filtered = f'{self.output_folder}/{self.barcode}_filt.fastq'
        print("filtering reads with filtong")
        print(f'minimum size is {min_length} and maximum allowed length is {max_length}')
        subprocess.call(f'filtlong --min_length {min_length} \
                        --max_length {max_length}\
                        --keep_percent 95 {processed_reads} > {path_to_filtered}',\
                        shell=True)
        print('Done')
        print(f"Filtered reads saved in {path_to_filtered}")
        return path_to_filtered

    def subset_reads(self,reads,coverage100):
        """coverage is calculated target coverage"""
        print(f'subsetting reads to achieve 50x with fraction = {coverage100}')

        n = random.randint(1,1000)

        path_to_subset = f'{self.output_folder}/{self.barcode}_subset.fastq'
                #cleaning up previous subset if it exists
        try:
            subprocess.call(f'rm {path_to_subset}')
        except FileNotFoundError:
            pass

        subprocess.call(f'seqtk sample -s{n} {reads} {coverage100} > \
                        {path_to_subset}', shell=True)
        
        return path_to_subset
    
    def assemble_w_Flye(self,reads):
        """where length is calculated size of plasmid"""
        print("Assembling reads with Flye")
        #cleaning up previous assemblies if they exist"""
        try:
            subprocess.call(f"rm -fR {self.output_folder}/assemblies",shell=True)
        except FileNotFoundError:
            pass

        subprocess.call(f'flye --nano-hq {reads} --threads 10 --out-dir {self.output_folder}/assemblies', shell=True,stderr=subprocess.DEVNULL)
        print("done")
        assembly_path = f'{self.output_folder}/assemblies/assembly.fasta'
        if not os.path.exists(assembly_path):
            assembly_path = ""

        return assembly_path
    
    def assemble_w_Canu(self,reads,length):
        """where length is calculated size of plasmid"""
        print("Assembling reads with canu")
        #cleaning up previous assemblies if they exist"""
        try:
            subprocess.call(f"rm -fR {self.output_folder}/assemblies",shell=True)
        except FileNotFoundError:
            pass

        subprocess.call(f'canu -p canu useGrid=False corOutCoverage=50 genomeSize={length}\
                         enableOEA=false -d {self.output_folder}/assemblies -nanopore {reads}', shell=True,stderr=subprocess.DEVNULL)
        print("done")
        assembly_path = f'{self.output_folder}/assemblies/canu.contigs.fasta'
        if not os.path.exists(assembly_path):
            assembly_path = ""

        return assembly_path

    def mash_dist(self,seq1,seq2):
    # calculates Mash distance between 2 sequences
        print("Running mash on 2 sequences...")
        subprocess.call("mkdir tmp",shell=True)
        rec1 = SeqRecord(seq1,id="contig_1")
        rec2 = SeqRecord(seq2,id="contig_2")
        SeqIO.write(rec1,"tmp/contig1","fasta")
        SeqIO.write(rec2,"tmp/contig2","fasta")
        mash_d = subprocess.check_output(["mash","dist","tmp/contig1","tmp/contig2"]).split(b'\t')[2]
        subprocess.call("rm -fR tmp",shell=True)
        return float(mash_d)

    def polish(self,reads,assembly,model):
        """polishing with medaka"""
        """model is basecalling model (string) to pass to medaka"""
        
        if model == None:
            model = "r1041_e82_400bps_sup_v4.2.0"
        else:
            model = model

        print(f'generating consensus with medaka using model {model}')
        #cleaning up data from previus runs
        try:
            subprocess.call(f'rm calls_to_draft.bam *.hdf qualities.fastq masked_qualities.fastq masked_consensus.fasta',shell=True)
            subprocess.call(f'mkdir {self.output_folder}/medaka',shell=True)
        except (FileNotFoundError, FileExistsError):
            pass

        subprocess.call(f'mini_align -r {assembly} -i {reads} -t 2 -m -p calls_to_draft'\
                        ,shell= True)
        subprocess.call(f'medaka consensus calls_to_draft.bam temp.hdf\
                        --model {model} --threads 8',shell=True)
        subprocess.call(f'medaka stitch --qualities *.hdf {assembly}\
                            qualities.fastq', shell=True)
        
        if os.path.isfile(f"qualities.fastq"):
            subprocess.call(f'fastq_masker -q 10 -i qualities.fastq \
                            -o masked_qualities.fastq',shell=True)
            subprocess.call(f'fastq_to_fasta -n -i masked_qualities.fastq \
                            -o masked_consensus.fasta',shell=True)
        else:
            print(f"medaka failed to generate a consensus for {self.barcode}")
        try:
            subprocess.call(f'mv masked_consensus.fasta {self.output_folder}/{self.name}_consensus.fasta',shell=True)
            subprocess.call(f'mv masked_qualities.fastq {self.output_folder}/{self.name}_qualities.fastq',shell=True)
        except FileNotFoundError:
            print('Output of Medaka not found')
            pass
        consensus_path = f'{self.output_folder}/{self.name}_consensus.fasta'
        qualities_path = f'{self.output_folder}/{self.name}_qualities.fastq'
        return (consensus_path,qualities_path)

    def check_concatemer(self,assembly_path,estimated_size):
        print("Checking if contig is a dimer")
        try:
            contig = list(SeqIO.parse(f"{assembly_path}","fasta"))[0] #extracting only 1st contig in case several were generated
            half = int(len(contig.seq)/2)
            # checking if contig is a duplicate. This often happens with long read assembly of circular DNA.
            if (1.9*estimated_size < len(contig.seq)):
                if self.mash_dist(contig.seq[:half],contig.seq[half:]) < 0.02:
                    print("contig is a duplicate - slicing in 2")
                    new_contig = contig.seq[:int(1.2*estimated_size-1)] #Leaving some overlapping ends to allow circularization
                    new_rec = SeqRecord(new_contig,id="contig_1")
                    subprocess.call("rm %s/assemblies/assembly_%d.fasta"%(path,i),shell = True)
                    SeqIO.write(new_rec,"%s/assemblies/assembly_%d"%(path,i),"fasta")
                else:
                    print("contig does not appear to be a dimer")
                    new_contig = contig.seq
            else:
                new_contig = contig.seq
        except (ValueError,IndexError):
            new_contig = ""
        return new_contig
    
    def circularize(self,contig_seq,len_match,max_mismatches):
        contig_seq = contig_seq
        len_match = int(len_match)
        max_mismatches = int(max_mismatches)
        print("trying to circularize...")
        sequence_string = contig_seq
        substring = sequence_string[0:len_match]
        string = sequence_string[len(substring):]
        indices = []
        for i in range(len(string) - len(substring) + 1):
            mismatches = sum(c1 != c2 for c1, c2 in zip(str(string[i:i+len(substring)]), str(substring)))
            if mismatches < max_mismatches:
                indices.append(i)
            else:
                pass
        if len(indices) == 1:
            index = indices[0] + len(substring)
            new_seq = string[0:index]
            overlap_seq = string[index:]
            if len(new_seq) or len(overlap_seq) < 100:
                print("overlaps at the contig ends < 100bp")
                sim_overlap = 0
            else:
                sim_overlap = self.mash_dist(new_seq, overlap_seq)

            print(f"sequence overlap mismatch is {sim_overlap*100/100}%")
            if sim_overlap < 0.05:
                print("contig ends overlap: Success!")
                # rotating sequence to allow for better correction of the joining site with medaka later on
                outseq = new_seq[100:] + new_seq[0:100]
            else:
                print("No overlap found: Cannot circularize")
                outseq = ""
        elif len(indices) < 1:
            print("No overlap found: Cannot circularize")
            outseq = ""
        elif len(indices) > 1:
            print("Multiple matches found: Cannot circularize")
            outseq = ""
        else:
            outseq = ""
            print("Unable to circularize")

        return Seq(outseq)
            

class PCR_Assembly(Assembly):
    def __init__(self, barcode, name, read_location,output_folder):
        Assembly.__init__(self, barcode, name, read_location, output_folder)
    
    def assemble(self,reads):
        print("Sample is a PCR product... Assembling directly with Canu")
        # reads is the path to  reads file
        """Assembling with Canu"""
        #subprocess.call(f'flye --nano-hq {reads} --threads 18 --out-dir {self.output_folder}/assemblies', shell=True)
        #assembly_path = f'{self.output_folder}/assemblies/assembly.fasta'
        subprocess.call(f'canu -p canu useGrid=False corOutCoverage=50 genomeSize=2000\
                         enableOEA=false -d {self.output_folder}/assemblies -nanopore {reads}', shell=True)
        print("done")
        assembly_path = f'{self.output_folder}/assemblies/canu.contigs.fasta'
        if not os.path.exists(assembly_path):
            assembly_path = ""

        return assembly_path
    
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


sample_df = pd.read_excel(args.xls,sheet_name='Feuil1',index_col=None, dtype={'Plasmid name':str, 'DNA type':str, 'Barcode':str})
#dropping unused sample slots
sample_df = sample_df.dropna(subset=['Barcode','Plasmid name'])
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
    sample_id = sub_df.iloc[0]['Plasmid name']
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

    print(subf,sample_id,sample_type)

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

        if coverage > 0:
            coverage50 = min(50/coverage,0.99)
            coverage150 = min(150/coverage,0.99)
        else:
            coverage50 = 0.99
            coverage150 = 0.99
            print(f'Unable to calculate coverage. Using all reads')

        filtered_reads = assembly.filter_reads(processed_reads,20000)
        reads_subset = assembly.subset_reads(filtered_reads,coverage50)
        contigs = assembly.assemble(reads_subset)
        if contigs != "":
            subprocess.call(f"mkdir {path}/final_assemblies/",shell=True)
            subprocess.call(f'mv {contigs} {path}/final_assemblies/{sample_id}.fasta',shell=True)
            final_contigs = f"{path}/final_assemblies/{sample_id}.fasta"
            assembly_type = "PCR"
    else:
        print("Unknown DNA type")
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







