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
import glob
import io
import logging
import os
import shutil
import subprocess
import sys
import traceback

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from assembly import Assembly, PCR_Assembly
from mapping import Mapping
from annotation import Annotation

# ──────────────────────────────────────────────────────────────────────────────
# In-memory log buffer — nothing is written to disk until the script finishes
# ──────────────────────────────────────────────────────────────────────────────

log_buffer = io.StringIO()

logging.basicConfig(
    stream=log_buffer,
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


class TeeStream:
    """Mirrors writes to both the original stream and the in-memory log buffer.

    Keeps stderr visible in the terminal during execution while simultaneously
    capturing it for the deferred log file.
    """

    def __init__(self, original: io.TextIOBase, buffer: io.StringIO) -> None:
        self.original = original
        self.buffer   = buffer

    def write(self, data: str) -> int:
        self.original.write(data)
        self.buffer.write(data)
        return len(data)

    def flush(self) -> None:
        self.original.flush()
        self.buffer.flush()

    def __getattr__(self, name):
        return getattr(self.original, name)


# Redirect stderr so subprocess output and Python tracebacks are captured too
_original_stderr = sys.stderr
sys.stderr = TeeStream(_original_stderr, log_buffer)

# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

argParser = argparse.ArgumentParser(description="NanoSeq — ONT plasmid assembly and annotation")

argParser.add_argument("-i", "--input",  required=True,
                       help="Input folder containing per-barcode subdirectories of ONT reads.")
argParser.add_argument("-x", "--xls",   required=True,
                       help="Sample information as Excel spreadsheet.")
argParser.add_argument("-o", "--output", required=True,
                       help="Output directory.")
argParser.add_argument("-m", "--model",  default="",
                       help="Medaka polishing model (e.g. r1041_e82_400bps_sup_v4.3.0). "
                            "Default: no polishing.")

ann_db = argParser.add_argument_group("Annotation databases")
ann_db.add_argument("-db", "--database", required=True,
                    help="Nucleotide feature database FASTA for blastn annotation.")
ann_db.add_argument("--protein-db", default=None, metavar="FASTA",
                    help="Protein database FASTA for blastp annotation of unannotated ORFs (optional).")

ann_tun = argParser.add_argument_group("Annotation parameters")
ann_tun.add_argument("--perc-identity",       type=float, default=80.0)
ann_tun.add_argument("--evalue",              type=float, default=1e-5)
ann_tun.add_argument("--word-size",           type=int,   default=11)
ann_tun.add_argument("--min-hit-coverage",    type=float, default=0.75)
ann_tun.add_argument("--min-orf-len",         type=int,   default=100)
ann_tun.add_argument("--blastp-evalue",       type=float, default=1e-5)
ann_tun.add_argument("--blastp-min-identity", type=float, default=30.0)
ann_tun.add_argument("--blastp-min-coverage", type=float, default=0.5)
ann_tun.add_argument("--threads",             type=int,   default=4)
ann_tun.add_argument("--figure-size",         type=int,   default=700)
ann_tun.add_argument("--restriction-enzymes", nargs="*",
                     default=["EcoRI", "BamHI", "HindIII"], metavar="ENZYME")

args = argParser.parse_args()

log_path = os.path.join(args.output, "nanoseq_run.log")


def _write_log(status: str) -> None:
    """Flush the in-memory buffer to disk exactly once."""
    sys.stderr = _original_stderr          # restore real stderr before writing
    os.makedirs(args.output, exist_ok=True)
    log_buffer.write(f"\n[STATUS] {status}\n")
    with open(log_path, "w") as fh:
        fh.write(log_buffer.getvalue())
    print(f"[NanoSeq] Log written to: {log_path}")


# ──────────────────────────────────────────────────────────────────────────────
# Main execution block
# ──────────────────────────────────────────────────────────────────────────────

try:
    logger.info("Starting NanoSeq with -i %s -x %s -o %s -m %s",
                args.input, args.xls, args.output, args.model)
    print(f"Starting NanoSeq with -i {args.input} -x {args.xls} "
          f"-o {args.output} -m {args.model}")

    # ── Output directory structure ────────────────────────────────────────────

    for subdir in ("failed_to_assemble", "circular_assemblies",
                   "failed_to_circularize", "PCR_assemblies", "read_qualities"):
        os.makedirs(f"{args.output}/{subdir}", exist_ok=True)

    # ── Sample sheet ──────────────────────────────────────────────────────────

    sample_df = pd.read_excel(
        args.xls, index_col=None,
        dtype={"Sample name": str, "DNA type": str, "Barcode": str},
    )
    sample_df = sample_df.dropna(subset=["Barcode", "Sample name"])
    sample_df["Barcode"]  = sample_df["Barcode"].str.lower()
    sample_df["DNA type"] = sample_df["DNA type"].str.lower()
    barcodes  = sample_df["Barcode"].tolist()
    bad_chars = [";", ":", "!", "*", " ", "/", ")", "(", "'"]

    folders = [f for f in os.listdir(args.input) if f in barcodes]
    logger.info("Folders to process: %s", folders)
    print(folders)

    # ── Annotation engine ─────────────────────────────────────────────────────

    annotator = Annotation(
        database=args.database,
        protein_db=args.protein_db,
        restriction_enzymes=args.restriction_enzymes,
        perc_identity=args.perc_identity,
        evalue=args.evalue,
        word_size=args.word_size,
        min_hit_coverage=args.min_hit_coverage,
        min_orf_len=args.min_orf_len,
        blastp_evalue=args.blastp_evalue,
        blastp_min_pident=args.blastp_min_identity,
        blastp_min_coverage=args.blastp_min_coverage,
        threads=args.threads,
        figure_size=args.figure_size,
    )

    # ── Per-sample assembly loop ──────────────────────────────────────────────

    for subf in folders:
        sub_df        = sample_df[sample_df["Barcode"] == subf]
        sample_id     = sub_df.iloc[0]["Sample name"]
        sample_type   = sub_df.iloc[0]["DNA type"].replace(" ", "")
        sample_length = sub_df.iloc[0]["Size (kb)"] * 1000

        sample_id = sample_id.replace(" ", "_")
        for c in bad_chars:
            if c in sample_id:
                sample_id = subf
                break

        path = f"{args.input}/{subf}"
        logger.info("Processing barcode=%s  sample_id=%s  type=%s", subf, sample_id, sample_type)
        print(f"Assembling {sub_df} with ID {sample_id} and type = {sample_type}")

        if glob.glob(f"{args.input}/{subf}/final_assemblies/*.fasta"):
            print("Already assembled, nothing to do")
            logger.info("Barcode %s already assembled — skipping.", subf)
            continue

        final_contigs  = ""
        assembly_type  = ""
        filtered_reads = ""
        coverage50     = 0.99

        match sample_type:
            case "plasmid":
                assembly        = Assembly(subf, sample_id, path, path)
                processed_reads = assembly.process_reads()
                QC              = assembly.quality_control(processed_reads)

                for report in QC:
                    try:
                        shutil.move(report, f"{args.output}/read_qualities")
                    except shutil.Error:
                        print("QC file already present")

                try:
                    length_raw, coverage_raw, total = assembly.get_length(processed_reads)
                    print(f"Calculated coverage is {coverage_raw}, predicted length is {length_raw}")
                    logger.info("Barcode %s — raw length=%s coverage=%s", subf, length_raw, coverage_raw)
                except IndexError:
                    logger.warning("Barcode %s — no reads found. Skipping.", subf)
                    print("No reads. Moving on to next sample")
                    continue

                max_length     = int(1.1 * length_raw) if length_raw != 0 else 1_000_000
                filtered_reads = assembly.filter_reads(processed_reads, max_length)
                length, coverage, total = assembly.get_length(filtered_reads)
                print(f"Calculated coverage after filtering is {coverage}")
                logger.info("Barcode %s — filtered coverage=%s", subf, coverage)

                if coverage > 0:
                    coverage50  = min(50  / coverage, 0.99)
                    coverage150 = min(150 / coverage, 0.99)
                else:
                    coverage50 = coverage150 = 0.99
                    logger.warning("Barcode %s — coverage is 0; using all reads.", subf)
                    print("Unable to calculate coverage. Using all reads")

                i, success = 0, False
                while not success and i <= 10:
                    i += 1
                    print(f"Iteration number {i}")
                    logger.debug("Barcode %s — assembly iteration %d", subf, i)

                    if i < 3:
                        reads_subset   = assembly.subset_reads(filtered_reads, coverage50)
                        contigs        = assembly.assemble_w_Canu(reads_subset, length)
                        if contigs:
                            contig_checked = assembly.check_concatemer(contigs, length)
                            circ_contig    = assembly.circularize(contig_checked, 30, 1)
                            if circ_contig:
                                circ_seq = SeqRecord(circ_contig, id=sample_id)
                                os.makedirs(f"{path}/final_assemblies", exist_ok=True)
                                SeqIO.write(circ_seq, f"{path}/final_assemblies/{sample_id}_circ.fasta", "fasta")
                                final_contigs = f"{path}/final_assemblies/{sample_id}_circ.fasta"
                                assembly_type = "circular"
                                success = True
                                logger.info("Barcode %s — circular assembly succeeded (Canu, iter %d).", subf, i)

                    elif i < 6:
                        reads_subset = assembly.subset_reads(filtered_reads, coverage50)
                        contigs      = assembly.assemble_w_Flye(reads_subset)
                        if contigs:
                            with open(f"{path}/assemblies/assembly_info.txt") as f:
                                data          = f.read().split("\n")[1].split("\t")
                                is_circular   = data[3]
                                contig_length = int(data[1])
                            if is_circular == "Y" and contig_length < max_length:
                                contig_checked = assembly.check_concatemer(contigs, length)
                                circ_seq       = SeqRecord(contig_checked, id=sample_id)
                                circ_seq       = circ_seq[100:] + circ_seq[0:100]
                                os.makedirs(f"{path}/final_assemblies", exist_ok=True)
                                SeqIO.write(circ_seq, f"{path}/final_assemblies/{sample_id}_circ.fasta", "fasta")
                                final_contigs = f"{path}/final_assemblies/{sample_id}_circ.fasta"
                                assembly_type = "circular"
                                success = True
                                logger.info("Barcode %s — circular assembly succeeded (Flye, iter %d).", subf, i)

                    elif i == 6:
                        print("Unable to circularize. Running one last Flye assembly")
                        logger.warning("Barcode %s — failed to circularize after 9 iterations; running final Flye.", subf)
                        reads_subset = assembly.subset_reads(filtered_reads, coverage50)
                        contigs      = assembly.assemble_w_Flye(reads_subset)
                        if contigs:
                            os.makedirs(f"{path}/final_assemblies", exist_ok=True)
                            final_contigs = f"{path}/final_assemblies/{sample_id}_non_circularized.fasta"
                            shutil.move(contigs, final_contigs)
                            assembly_type = "linear"
                            logger.info("Barcode %s — saved as linear (non-circularized).", subf)
                        else:
                            print(f"Unable to assemble sample {subf}. Moving on.")
                            logger.error("Barcode %s — assembly failed entirely.", subf)
                            assembly_type = "failed"

            case "large":
                assembly        = Assembly(subf, sample_id, path, path)
                processed_reads = assembly.process_reads()
                QC              = assembly.quality_control(processed_reads)
                for report in QC:
                    try:
                        shutil.move(report, f"{args.output}/read_qualities")
                    except shutil.Error:
                        print("QC file already present")

                filtered_reads = assembly.filter_reads(processed_reads, 1_000_000)
                reads_subset   = assembly.subset_reads(filtered_reads, 0.99)
                contigs        = assembly.assemble_w_Flye(reads_subset)

                if contigs:
                    with open(f"{path}/assemblies/assembly_info.txt") as f:
                        data          = f.read().split("\n")[1].split("\t")
                        is_circular   = data[3]
                        contig_length = int(data[1])

                    if is_circular == "Y" and contig_length < 1_000_000:
                        contig_checked = assembly.check_concatemer(contigs, contig_length)
                        circ_seq       = SeqRecord(contig_checked, id=sample_id)
                        circ_seq       = circ_seq[100:] + circ_seq[0:100]
                        os.makedirs(f"{path}/final_assemblies", exist_ok=True)
                        SeqIO.write(circ_seq, f"{path}/final_assemblies/{sample_id}_circ.fasta", "fasta")
                        final_contigs = f"{path}/final_assemblies/{sample_id}_circ.fasta"
                        assembly_type = "circular"
                        logger.info("Barcode %s (large) — circular assembly succeeded.", subf)
                    elif contig_length < 1_000_000:
                        print("Contig is predicted linear")
                        os.makedirs(f"{path}/final_assemblies", exist_ok=True)
                        final_contigs = f"{path}/final_assemblies/{sample_id}_non_circularized.fasta"
                        shutil.move(contigs, final_contigs)
                        assembly_type = "linear"
                        logger.info("Barcode %s (large) — linear assembly.", subf)
                    else:
                        print("Unable to assemble large plasmid. Moving on.")
                        logger.error("Barcode %s (large) — assembly failed.", subf)
                        assembly_type = "failed"

            case "pcr":
                assembly        = PCR_Assembly(subf, sample_id, path, path)
                processed_reads = assembly.process_reads()
                QC              = assembly.quality_control(processed_reads)
                for report in QC:
                    try:
                        shutil.move(report, f"{args.output}/read_qualities")
                    except shutil.Error:
                        print("QC file already present")

                length, coverage, total = assembly.get_length(processed_reads)
                if sample_length != "":
                    length   = sample_length
                    coverage = float(total / int(sample_length))
                else:
                    print("No length given. Guessing from distribution")
                    try:
                        coverage = float(total / int(length))
                    except ZeroDivisionError:
                        logger.warning("Barcode %s (PCR) — unable to guess length; defaulting to 1000 bp.", subf)
                        print("Unable to guess sample length. Setting to 1000 bp")
                        length   = 1000
                        coverage = float(total / int(length))

                print(f"Size of fragment from metadata is {length} and calculated coverage is {coverage}")
                logger.info("Barcode %s (PCR) — length=%s coverage=%s", subf, length, coverage)

                coverage     = float(total / int(sample_length))
                coverage150  = min(150 / coverage, 0.99) if coverage > 0 else 0.99
                filtered_reads = assembly.filter_reads(processed_reads, 20_000)

                k, contigs = 0, ""
                while k < 3 and not contigs:
                    k += 1
                    reads_subset = assembly.subset_reads(filtered_reads, coverage150)
                    contigs      = assembly.assemble(reads_subset)
                    if contigs:
                        os.makedirs(f"{path}/final_assemblies", exist_ok=True)
                        final_contigs = f"{path}/final_assemblies/{sample_id}.fasta"
                        shutil.move(contigs, final_contigs)
                        assembly_type = "PCR"
                        logger.info("Barcode %s (PCR) — assembly succeeded (attempt %d).", subf, k)
                    else:
                        print(f"Unable to assemble sample {subf}. Moving on.")
                        logger.error("Barcode %s (PCR) — assembly failed (attempt %d).", subf, k)
                        assembly_type = "failed"

            case _:
                print("Unknown DNA type: not processing.")
                logger.warning("Barcode %s — unknown DNA type '%s'. Skipping.", subf, sample_type)
                os.makedirs(f"{args.output}/failed_to_assemble/{subf}", exist_ok=True)
            subprocess.run(
                f"zcat {path}/*.fastq.gz > {args.output}/failed_to_assemble/{subf}/{subf}.fastq",
                shell=True,
            )
            logger.warning("Barcode %s — moved raw reads to failed_to_assemble/.", subf)

        # ── Post-assembly: polish → map → annotate ────────────────────────────

        if final_contigs and args.model:
            reads_subset50 = assembly.subset_reads(filtered_reads, coverage50)
            consensus, qualities = assembly.polish(reads_subset50, final_contigs, args.model)
            logger.info("Barcode %s — polishing complete.", subf)
            os.makedirs(f"{args.output}/BAM_files", exist_ok=True)
            mappings = Mapping(subf, consensus, processed_reads, f"{args.output}/BAM_files")
            mappings.minimap()
            logger.info("Barcode %s — mapping complete.", subf)

            match assembly_type:
                case "circular":
                    dest = f"{args.output}/circular_assemblies/{subf}"
                    os.makedirs(dest, exist_ok=True)
                    shutil.move(consensus, dest)
                    shutil.move(qualities, dest)
                case "linear":
                    dest = f"{args.output}/failed_to_circularize/{subf}"
                    os.makedirs(dest, exist_ok=True)
                    shutil.move(consensus, dest)
                    shutil.move(qualities, dest)
                case "PCR":
                    dest = f"{args.output}/PCR_assemblies/{subf}"
                    os.makedirs(dest, exist_ok=True)
                    shutil.move(consensus, dest)
                    shutil.move(qualities, dest)
                case _:
                    logger.error("Barcode %s — unknown assembly type '%s' after polishing.", subf, assembly_type)
                    print("Unknown assembly type. Check assembler output.")

        elif final_contigs and not args.model:
            os.makedirs(f"{args.output}/BAM_files", exist_ok=True)
            mappings = Mapping(subf, final_contigs, processed_reads, f"{args.output}/BAM_files")
            mappings.minimap()
            logger.info("Barcode %s — mapping complete (no polishing).", subf)

            match assembly_type:
                case "circular" | "linear":
                    dest_root = ("circular_assemblies" if assembly_type == "circular"
                                 else "failed_to_circularize")
                    dest = f"{args.output}/{dest_root}/{subf}"
                    os.makedirs(dest, exist_ok=True)
                    shutil.move(final_contigs, dest)

                    moved_fasta = f"{dest}/{os.path.basename(final_contigs)}"
                    ann_outdir  = f"{dest}/annotation"
                    print(f"Annotating {assembly_type} assembly: {moved_fasta}")
                    logger.info("Barcode %s — annotating %s assembly.", subf, assembly_type)
                    annotator.annotate(moved_fasta, ann_outdir)

                case "PCR":
                    dest = f"{args.output}/PCR_assemblies/{subf}"
                    os.makedirs(dest, exist_ok=True)
                    shutil.move(final_contigs, dest)

                case _:
                    logger.error("Barcode %s — unknown assembly type '%s' after mapping.", subf, assembly_type)
                    print("Unknown assembly type. Check assembler output.")

        else:
            os.makedirs(f"{args.output}/failed_to_assemble/{subf}", exist_ok=True)
            subprocess.run(
                f"zcat {path}/*.fastq.gz > {args.output}/failed_to_assemble/{subf}/{subf}.fastq",
                shell=True,
            )
            logger.warning("Barcode %s — moved raw reads to failed_to_assemble/.", subf)

    print("All done!")
    logger.info("NanoSeq finished successfully.")
    _write_log("COMPLETED SUCCESSFULLY")

except Exception:
    logger.critical("Unhandled exception — script aborted.\n%s", traceback.format_exc())
    _write_log("CRASHED — see traceback above")
    sys.exit(1)
