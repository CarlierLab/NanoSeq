"""
Plasmid annotation module — NanoSeq
====================================
Copyright 2026 Aurelien Carlier (aurelien.carlier@inrae.fr)
https://github.com/CarlierLab/NanoSeq

This file is part of NanoSeq. NanoSeq is free software: you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version. NanoSeq is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details. You should have received a copy of
the GNU General Public License along with NanoSeq.
If not, see <http://www.gnu.org/licenses/>.

Public API
----------
    from annotation import Annotation

    ann = Annotation(
        database="features.fasta",
        protein_db="proteins.fasta",          # optional
        restriction_enzymes=["EcoRI","BamHI"],# optional
    )
    annotated_records = ann.annotate("assembly.fasta", "results/sample1/")
"""

import csv
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from bokeh.io import output_file, save as bokeh_save
from bokeh.models import ColumnDataSource, HoverTool, Label
from bokeh.plotting import figure as bokeh_figure

log = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# Dependency checks
# ──────────────────────────────────────────────────────────────────────────────

REQUIRED_TOOLS       = ["blastn", "makeblastdb"]
REQUIRED_TOOLS_BLASTP = ["blastp"]


def check_dependencies(use_blastp: bool = False) -> None:
    tools   = REQUIRED_TOOLS + (REQUIRED_TOOLS_BLASTP if use_blastp else [])
    missing = [t for t in tools if not shutil.which(t)]
    if missing:
        sys.exit(
            f"Error: the following tools are not in PATH: {', '.join(missing)}\n"
            "Install BLAST+ from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        )


# ──────────────────────────────────────────────────────────────────────────────
# BLAST helpers
# ──────────────────────────────────────────────────────────────────────────────

BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length",
    "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "qlen", "slen",
]
BLAST_FMT = "6 " + " ".join(BLAST_FIELDS)


def make_blast_db(fasta_path: str, db_prefix: str) -> None:
    """Create (or reuse) a BLAST nucleotide database."""
    marker = db_prefix + ".nhr"
    if os.path.exists(marker):
        log.info("Nucleotide BLAST database already exists, skipping build (%s).", db_prefix)
        return
    cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl",
           "-out", db_prefix, "-parse_seqids"]
    log.info("Building nucleotide BLAST database …")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"makeblastdb failed:\n{result.stderr}")


def make_blast_protein_db(fasta_path: str, db_prefix: str) -> None:
    """Create (or reuse) a BLAST protein database."""
    marker = db_prefix + ".phr"
    if os.path.exists(marker):
        log.info("Protein BLAST database already exists, skipping build (%s).", db_prefix)
        return
    cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "prot",
           "-out", db_prefix, "-parse_seqids"]
    log.info("Building protein BLAST database …")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"makeblastdb (prot) failed:\n{result.stderr}")


def run_blastn(query_path: str, db_prefix: str, out_path: str,
               perc_identity: float, evalue: float,
               word_size: int, threads: int) -> None:
    """Run blastn with tabular output (format 6)."""
    cmd = [
        "blastn", "-query", query_path, "-db", db_prefix,
        "-out", out_path, "-outfmt", BLAST_FMT,
        "-perc_identity", str(perc_identity),
        "-evalue", str(evalue), "-word_size", str(word_size),
        "-num_threads", str(threads), "-dust", "no", "-soft_masking", "false",
    ]
    log.info("Running blastn …")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"blastn failed:\n{result.stderr}")


def run_blastp(query_path: str, db_prefix: str, out_path: str,
               evalue: float, threads: int) -> None:
    """Run blastp with tabular output (format 6)."""
    cmd = [
        "blastp", "-query", query_path, "-db", db_prefix,
        "-out", out_path, "-outfmt", BLAST_FMT,
        "-evalue", str(evalue), "-num_threads", str(threads), "-seg", "no",
    ]
    log.info("Running blastp …")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"blastp failed:\n{result.stderr}")


def parse_blast_table(path: str) -> list[dict]:
    """Parse BLAST tabular output into a list of hit dicts."""
    hits = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < len(BLAST_FIELDS):
                continue
            h = dict(zip(BLAST_FIELDS, parts))
            for int_field in ("length", "qstart", "qend", "sstart", "send", "qlen", "slen"):
                h[int_field] = int(h[int_field])
            for float_field in ("pident", "evalue", "bitscore"):
                h[float_field] = float(h[float_field])
            h["strand"] = 1 if h["sstart"] <= h["send"] else -1
            if h["qstart"] > h["qend"]:
                h["qstart"], h["qend"] = h["qend"], h["qstart"]
            hits.append(h)
    return hits


# ──────────────────────────────────────────────────────────────────────────────
# Circular-sequence helpers
# ──────────────────────────────────────────────────────────────────────────────

def doubled_id(seq_id: str) -> str:
    return seq_id + "__doubled__"


def write_doubled_fasta(records: list[SeqRecord], path: str) -> None:
    """Write a FASTA where every sequence is duplicated (for circularity)."""
    doubled = [
        SeqRecord(rec.seq + rec.seq, id=doubled_id(rec.id), description="")
        for rec in records
    ]
    with open(path, "w") as fh:
        SeqIO.write(doubled, fh, "fasta")


def remap_hit_to_circular(hit: dict, seq_len: int) -> dict | None:
    """Remap a hit from the doubled query back to circular coordinates."""
    qs = hit["qstart"] - 1
    qe = hit["qend"]
    if qs >= seq_len:
        return None
    h = hit.copy()
    h["_qs0"]   = qs % seq_len
    h["_qe0"]   = qe
    h["_wraps"] = qe > seq_len
    return h


def hits_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def select_non_overlapping(hits: list[dict], seq_len: int) -> list[dict]:
    """Greedy best-first selection of non-overlapping hits."""
    sorted_hits = sorted(hits, key=lambda h: h["bitscore"], reverse=True)
    selected:  list[dict]         = []
    occupied:  list[tuple[int,int]] = []

    def _overlaps_any(s: int, e: int) -> bool:
        for os_, oe_ in occupied:
            if hits_overlap(s, e, os_, oe_):
                return True
            if oe_ > seq_len:
                if hits_overlap(s, e, 0, oe_ - seq_len):
                    return True
        return False

    for h in sorted_hits:
        if _overlaps_any(h["_qs0"], h["_qe0"]):
            continue
        selected.append(h)
        occupied.append((h["_qs0"], h["_qe0"]))
    return selected


# ──────────────────────────────────────────────────────────────────────────────
# ORF finding (circular-aware)
# ──────────────────────────────────────────────────────────────────────────────

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_orfs(seq: str, min_orf_len: int = 100) -> list[dict]:
    """Find all ORFs on both strands of a circular sequence."""
    seq_len = len(seq)
    doubled = (seq + seq).upper()
    orfs: list[dict] = []

    def _scan_strand(nuc: str, strand: int) -> None:
        for frame in range(3):
            i = frame
            while i + 3 <= len(nuc):
                if nuc[i:i + 3] != START_CODON:
                    i += 3
                    continue
                j = i + 3
                while j + 3 <= len(nuc):
                    codon = nuc[j:j + 3]
                    if codon in STOP_CODONS:
                        if j + 3 - i >= min_orf_len:
                            if strand == 1:
                                raw_start, raw_end = i, j + 3
                            else:
                                raw_start = len(nuc) - (j + 3)
                                raw_end   = len(nuc) - i
                            if raw_start < seq_len:
                                orfs.append({
                                    "start":  raw_start % seq_len,
                                    "end":    raw_end,
                                    "strand": strand,
                                    "wraps":  raw_end > seq_len,
                                })
                        i = j + 3
                        break
                    j += 3
                else:
                    i += 3

    _scan_strand(doubled, 1)
    _scan_strand(str(Seq(doubled).reverse_complement()), -1)
    return orfs


def orf_contains_hit(orf: dict, hit_start: int, hit_end: int, hit_strand: int) -> bool:
    if orf["strand"] != hit_strand:
        return False
    return orf["start"] <= hit_start and orf["end"] >= hit_end


def find_containing_orf(orfs: list[dict], hit_start: int,
                        hit_end: int, hit_strand: int) -> dict | None:
    """Return the shortest ORF that fully contains the hit, or None."""
    candidates = [o for o in orfs if orf_contains_hit(o, hit_start, hit_end, hit_strand)]
    return min(candidates, key=lambda o: o["end"] - o["start"]) if candidates else None


# ──────────────────────────────────────────────────────────────────────────────
# Feature-type inference
# ──────────────────────────────────────────────────────────────────────────────

FEATURE_TYPE_KEYWORDS: dict[str, str] = {
    "promoter": "promoter", "prom": "promoter",
    "terminator": "terminator",
    "ori": "rep_origin", "origin": "rep_origin", "rep_origin": "rep_origin",
    "primer": "primer_bind",
    "ltr": "LTR",
    "rrna": "rRNA", "trna": "tRNA",
    "regulatory": "regulatory", "enhancer": "enhancer",
    "insulator": "insulator", "misc_rna": "misc_RNA",
}

CDS_HINT_KEYWORDS = {
    "gene", "cds", "protein", "resistance", "replication", "rep",
    "integrase", "recombinase", "transposase", "lacza", "gfp", "rfp",
    "fluorescent", "antibiotic", "ampicillin", "kanamycin", "chloramphenicol",
    "tetracycline", "neomycin", "hygromycin", "puromycin", "blasticidin",
    "zeocin", "spectinomycin", "gentamicin", "streptomycin", "bleomycin",
    "luciferase", "mcherry", "egfp", "cfp", "yfp", "his", "flag",
}


def infer_feature_type(description: str, in_orf: bool) -> str:
    name_lc   = description.lower()
    feat_type = "misc_feature"
    if in_orf:
        feat_type = "CDS"
    for kw, ftype in FEATURE_TYPE_KEYWORDS.items():
        if kw in name_lc:
            feat_type = ftype
    return feat_type


def label_from_id(seq_id: str) -> str:
    return seq_id.split("|")[-1].replace("_", " ").strip()


# ──────────────────────────────────────────────────────────────────────────────
# Build SeqFeature objects
# ──────────────────────────────────────────────────────────────────────────────

def _make_location(feat_start: int, feat_end: int,
                   strand: int, seq_len: int, wraps: bool) -> FeatureLocation | CompoundLocation:
    if wraps:
        p1e = min(feat_end, seq_len)
        p2e = feat_end - seq_len
        if strand == 1:
            return CompoundLocation([
                FeatureLocation(feat_start, p1e, strand=strand),
                FeatureLocation(0, p2e, strand=strand),
            ])
        else:
            return CompoundLocation([
                FeatureLocation(0, p2e, strand=strand),
                FeatureLocation(feat_start, p1e, strand=strand),
            ])
    return FeatureLocation(feat_start, min(feat_end, seq_len), strand=strand)


def build_feature(hit: dict, seq_len: int, orfs: list[dict],
                  min_hit_coverage: float,
                  feature_descriptions: dict[str, str]) -> SeqFeature | None:
    """Convert a BLAST hit to a SeqFeature, or None if coverage too low."""
    qs, qe, strand = hit["_qs0"], hit["_qe0"], hit["strand"]
    hit_coverage   = hit["length"] / hit["slen"] if hit["slen"] > 0 else 0
    if hit_coverage < min_hit_coverage:
        log.debug("Skipping %s (coverage %.1f%% < %.1f%%)",
                  hit["sseqid"], hit_coverage * 100, min_hit_coverage * 100)
        return None

    containing_orf = find_containing_orf(orfs, qs, qe, strand)
    in_orf         = containing_orf is not None
    label          = label_from_id(hit["sseqid"])
    description    = feature_descriptions.get(hit["sseqid"], "")
    feature_type   = infer_feature_type(description, in_orf)

    if feature_type == "CDS" and containing_orf is not None:
        fs, fe, wraps = containing_orf["start"], containing_orf["end"], containing_orf["wraps"]
    else:
        fs, fe, wraps = qs, qe, hit["_wraps"]

    loc = _make_location(fs, fe, strand, seq_len, wraps)

    qualifiers: dict[str, list] = {
        "label": [label],
        "note":  [f"hit={hit['sseqid']}; identity={hit['pident']:.1f}%; "
                  f"coverage={hit_coverage * 100:.1f}%; evalue={hit['evalue']:.2e}"],
    }
    if description:
        qualifiers["description"] = [description]
    if feature_type == "CDS":
        qualifiers["product"] = [label]
        if not in_orf:
            qualifiers["note"].append("WARNING: no enclosing ORF found")

    return SeqFeature(location=loc, type=feature_type, qualifiers=qualifiers)


# ──────────────────────────────────────────────────────────────────────────────
# blastp annotation of unannotated ORFs
# ──────────────────────────────────────────────────────────────────────────────

_ORF_ID_SEP = "__orf__"


def _orf_query_id(plasmid_id: str, orf: dict) -> str:
    return f"{plasmid_id}{_ORF_ID_SEP}{orf['start']}__{orf['end']}__{orf['strand']}"


def _parse_orf_query_id(qseqid: str) -> tuple[str, int, int, int] | None:
    if _ORF_ID_SEP not in qseqid:
        return None
    plasmid_id, coords = qseqid.split(_ORF_ID_SEP, 1)
    parts = coords.split("__")
    if len(parts) != 3:
        return None
    try:
        return plasmid_id, int(parts[0]), int(parts[1]), int(parts[2])
    except ValueError:
        return None


def get_unannotated_orfs(orfs: list[dict], features: list[SeqFeature],
                         seq_len: int) -> list[dict]:
    """Return ORFs that do not overlap any existing annotated feature."""
    annotated: list[tuple[int, int, int]] = []
    for feat in features:
        st = feat.location.strand if feat.location.strand is not None else 0
        for part in feat.location.parts:
            annotated.append((int(part.start), int(part.end), st))

    def _overlaps(orf: dict) -> bool:
        for fs, fe, fst in annotated:
            if fst != 0 and fst != orf["strand"]:
                continue
            if hits_overlap(orf["start"], min(orf["end"], seq_len), fs, fe):
                return True
            if orf["end"] > seq_len:
                if hits_overlap(0, orf["end"] - seq_len, fs, fe):
                    return True
        return False

    return [o for o in orfs if not _overlaps(o)]


def write_orf_protein_fasta(plasmid_id: str, orfs: list[dict],
                            seq: str, path: str) -> list[dict]:
    """Translate unannotated ORFs and write protein FASTA for blastp."""
    doubled = seq + seq
    written, records = [], []
    for orf in orfs:
        nuc = doubled[orf["start"]:orf["end"]]
        if orf["strand"] == -1:
            nuc = str(Seq(nuc).reverse_complement())
        try:
            prot = str(Seq(nuc).translate(to_stop=False))
        except Exception:
            continue
        if "*" in prot[:-1]:
            continue
        prot = prot.rstrip("*")
        if not prot:
            continue
        records.append(SeqRecord(Seq(prot), id=_orf_query_id(plasmid_id, orf), description=""))
        written.append(orf)
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "fasta")
    return written


def annotate_orfs_with_blastp(
    record: SeqRecord,
    orfs: list[dict],
    existing_features: list[SeqFeature],
    protein_db_prefix: str,
    workdir: str,
    evalue: float,
    min_pident: float,
    min_coverage: float,
    threads: int,
    protein_descriptions: dict[str, str] | None = None,
) -> list[SeqFeature]:
    """Search unannotated ORFs with blastp and return new CDS features."""
    seq_len  = len(record.seq)
    seq_str  = str(record.seq).upper()

    unannotated = get_unannotated_orfs(orfs, existing_features, seq_len)
    if not unannotated:
        log.info("  No unannotated ORFs to search with blastp.")
        return []
    log.info("  %d unannotated ORF(s) will be queried with blastp.", len(unannotated))

    prot_query  = os.path.join(workdir, f"{record.id}_orf_proteins.fasta")
    written_orfs = write_orf_protein_fasta(record.id, unannotated, seq_str, prot_query)
    if not written_orfs:
        log.info("  No valid ORF translations produced; skipping blastp.")
        return []

    blastp_out = os.path.join(workdir, f"{record.id}_blastp.tsv")
    run_blastp(prot_query, protein_db_prefix, blastp_out, evalue=evalue, threads=threads)

    best: dict[str, dict] = {}
    for hit in parse_blast_table(blastp_out):
        if hit["pident"] < min_pident:
            continue
        cov = hit["length"] / hit["slen"] if hit["slen"] > 0 else 0
        if cov < min_coverage:
            continue
        qid = hit["qseqid"]
        if qid not in best or hit["bitscore"] > best[qid]["bitscore"]:
            best[qid] = hit
            best[qid]["_coverage"] = cov

    new_features: list[SeqFeature] = []
    for hit in best.values():
        parsed = _parse_orf_query_id(hit["qseqid"])
        if parsed is None:
            continue
        _, orf_start, orf_end, orf_strand = parsed
        cov   = hit["_coverage"]
        loc   = _make_location(orf_start, orf_end, orf_strand, seq_len, orf_end > seq_len)
        label = label_from_id(hit["sseqid"])
        desc  = (protein_descriptions or {}).get(hit["sseqid"], "")
        qualifiers = {
            "label":   [label],
            "product": [label],
            "note":    [f"blastp hit={hit['sseqid']}; identity={hit['pident']:.1f}%; "
                        f"coverage={cov * 100:.1f}%; evalue={hit['evalue']:.2e}"],
        }
        if desc:
            qualifiers["description"] = [desc]
        new_features.append(SeqFeature(location=loc, type="CDS", qualifiers=qualifiers))

    log.info("  %d CDS feature(s) added from blastp.", len(new_features))
    return new_features


# ──────────────────────────────────────────────────────────────────────────────
# Per-record annotation pipeline
# ──────────────────────────────────────────────────────────────────────────────

def annotate_record(
    record: SeqRecord,
    raw_hits: list[dict],
    min_hit_coverage: float,
    min_orf_len: int,
    feature_descriptions: dict[str, str],
    protein_db_prefix: str | None = None,
    workdir: str | None = None,
    blastp_evalue: float = 1e-5,
    blastp_min_pident: float = 30.0,
    blastp_min_coverage: float = 0.5,
    threads: int = 4,
    protein_descriptions: dict[str, str] | None = None,
) -> SeqRecord:
    """Annotate a single plasmid SeqRecord (blastn + optional blastp pass)."""
    seq_len = len(record.seq)
    dbl_id  = doubled_id(record.id)

    record_hits = [h for h in raw_hits if h["qseqid"] == dbl_id]
    if not record_hits:
        log.info("  No BLAST hits for %s.", record.id)
        return _finalize_record(record)

    remapped = [r for h in record_hits
                for r in [remap_hit_to_circular(h, seq_len)] if r is not None]
    log.info("  %d hits after circular remapping.", len(remapped))

    selected = select_non_overlapping(remapped, seq_len)
    log.info("  %d non-overlapping hits selected.", len(selected))

    seq_str = str(record.seq).upper()
    orfs    = find_orfs(seq_str, min_orf_len=min_orf_len)
    log.info("  %d ORFs found (≥%d nt).", len(orfs), min_orf_len)

    features: list[SeqFeature] = []
    for hit in selected:
        feat = build_feature(hit, seq_len, orfs, min_hit_coverage, feature_descriptions)
        if feat is not None:
            features.append(feat)
    log.info("  %d features annotated (blastn).", len(features))

    if protein_db_prefix and workdir:
        features.extend(annotate_orfs_with_blastp(
            record=record, orfs=orfs, existing_features=features,
            protein_db_prefix=protein_db_prefix, workdir=workdir,
            evalue=blastp_evalue, min_pident=blastp_min_pident,
            min_coverage=blastp_min_coverage, threads=threads,
            protein_descriptions=protein_descriptions,
        ))

    return _finalize_record(record, features)


def _finalize_record(record: SeqRecord,
                     features: list[SeqFeature] | None = None) -> SeqRecord:
    return SeqRecord(
        seq=record.seq, id=record.id, name=record.name[:16],
        description=record.description, features=features or [],
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )


# ──────────────────────────────────────────────────────────────────────────────
# Restriction site detection
# ──────────────────────────────────────────────────────────────────────────────

def find_restriction_sites(record: SeqRecord,
                           enzyme_names: list[str]) -> list[SeqFeature]:
    """
    Detect restriction sites in a circular sequence using Bio.Restriction.

    Returns one misc_feature per cut site, tagged with ``bound_moiety``
    so the plotting code can distinguish them from annotation features.
    """
    import Bio.Restriction as Restriction

    enzymes = []
    for name in enzyme_names:
        try:
            enzymes.append(getattr(Restriction, name))
        except AttributeError:
            log.warning("Unknown restriction enzyme '%s' – skipped.", name)
    if not enzymes:
        return []

    rb      = Restriction.RestrictionBatch(enzymes)
    ana     = Restriction.Analysis(rb, record.seq, linear=False)
    results = ana.full()

    features: list[SeqFeature] = []
    for enzyme, positions in results.items():
        if not positions:
            continue
        site_len = len(enzyme.site)
        for pos in positions:
            pos0 = pos - 1
            end0 = min(pos0 + site_len, len(record.seq))
            features.append(SeqFeature(
                location=FeatureLocation(pos0, end0, strand=0),
                type="misc_feature",
                qualifiers={
                    "label":        [str(enzyme)],
                    "bound_moiety": [str(enzyme)],
                    "note":         [f"Restriction site; enzyme={enzyme}; "
                                     f"recognition_seq={enzyme.site}; position={pos}"],
                },
            ))
        log.info("  %s: %d site(s) → %s", enzyme, len(positions),
                 ", ".join(str(p) for p in positions))
    return features


# ──────────────────────────────────────────────────────────────────────────────
# Summary TSV
# ──────────────────────────────────────────────────────────────────────────────

def write_summary(records: list[SeqRecord], out_path: str) -> None:
    """Write a flat TSV summarising all annotated features."""
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["plasmid", "feature_type", "start_1based", "end_1based",
                         "strand", "label", "note", "description"])
        for rec in records:
            for feat in rec.features:
                strand = "+" if feat.location.strand == 1 else "-"
                writer.writerow([
                    rec.id, feat.type,
                    int(feat.location.start) + 1, int(feat.location.end),
                    strand,
                    feat.qualifiers.get("label",       [""])[0],
                    feat.qualifiers.get("note",        [""])[0],
                    feat.qualifiers.get("description", [""])[0],
                ])
    log.info("Summary written to %s", out_path)


# ──────────────────────────────────────────────────────────────────────────────
# Plasmid map plotting
# ──────────────────────────────────────────────────────────────────────────────

FEATURE_COLORS: dict[str, str] = {
    "CDS":          "#4A90D9",
    "promoter":     "#E67E22",
    "terminator":   "#E74C3C",
    "rep_origin":   "#2ECC71",
    "misc_feature": "#95A5A6",
    "primer_bind":  "#9B59B6",
    "LTR":          "#1ABC9C",
    "rRNA":         "#F39C12",
    "tRNA":         "#D35400",
    "regulatory":   "#C0392B",
    "enhancer":     "#8E44AD",
}
_DEFAULT_FEATURE_COLOR = "#7F8C8D"

_BACKBONE_INNER = 0.61
_BACKBONE_OUTER = 0.64
_FWD_INNER      = 0.64
_FWD_OUTER      = 0.84
_REV_INNER      = 0.41
_REV_OUTER      = 0.61
_RS_COLOR       = "#C0392B"
_RS_R_INNER     = _FWD_OUTER + 0.02
_RS_R_OUTER     = _FWD_OUTER + 0.14
_RS_R_LABEL     = _FWD_OUTER + 0.20


def _nt_to_angle(pos: int, seq_len: int) -> float:
    return math.pi / 2 - 2 * math.pi * pos / seq_len


def build_circular_bokeh_plot(record: SeqRecord, figure_size: int = 700):
    """Build a Bokeh circular plasmid map from an annotated SeqRecord."""
    seq_len = len(record.seq)

    p = bokeh_figure(
        width=figure_size, height=figure_size,
        x_range=(-1.5, 1.5), y_range=(-1.5, 1.5),
        title=f"{record.id}   ({seq_len:,} bp)",
        tools="pan,wheel_zoom,reset,save",
        toolbar_location="above",
    )
    p.axis.visible = False
    p.grid.visible = False
    p.outline_line_color = None

    # Backbone ring
    p.annular_wedge(x=0, y=0,
                    inner_radius=_BACKBONE_INNER, outer_radius=_BACKBONE_OUTER,
                    start_angle=0, end_angle=2 * math.pi,
                    color="#BBBBBB", line_color=None)

    # Position ticks at 0 %, 25 %, 50 %, 75 %
    for i in range(4):
        pos   = i * seq_len // 4
        angle = _nt_to_angle(pos, seq_len)
        ca, sa = math.cos(angle), math.sin(angle)
        p.line([_BACKBONE_OUTER * ca, (_BACKBONE_OUTER + 0.06) * ca],
               [_BACKBONE_OUTER * sa, (_BACKBONE_OUTER + 0.06) * sa],
               line_width=1.5, color="#666666")
        p.add_layout(Label(
            x=(_BACKBONE_OUTER + 0.13) * ca,
            y=(_BACKBONE_OUTER + 0.13) * sa,
            text=f"{pos:,} bp",
            text_align="center", text_baseline="middle",
            text_font_size="10px", text_color="#555555",
        ))

    # Feature wedges
    def _empty():
        return {k: [] for k in ("start_angle","end_angle","inner","outer",
                                "color","label","description","note","ftype")}

    fwd, rev = _empty(), _empty()
    for feat in record.features:
        if feat.qualifiers.get("bound_moiety"):
            continue                              # restriction sites drawn separately
        strand = feat.location.strand if feat.location.strand is not None else 1
        color  = FEATURE_COLORS.get(feat.type, _DEFAULT_FEATURE_COLOR)
        track  = fwd if strand >= 0 else rev
        inner  = _FWD_INNER if strand >= 0 else _REV_INNER
        outer  = _FWD_OUTER if strand >= 0 else _REV_OUTER
        for part in feat.location.parts:
            track["start_angle"].append(_nt_to_angle(int(part.start), seq_len))
            track["end_angle"].append(  _nt_to_angle(int(part.end),   seq_len))
            track["inner"].append(inner)
            track["outer"].append(outer)
            track["color"].append(color)
            track["label"].append(feat.qualifiers.get("label",       [""])[0])
            track["description"].append(feat.qualifiers.get("description", [""])[0])
            track["note"].append(  feat.qualifiers.get("note",        [""])[0])
            track["ftype"].append(feat.type)

    hover = HoverTool(tooltips="""
        <div style="max-width:360px;font-family:sans-serif;padding:4px;">
            <b>@label</b> <span style="color:#888;font-size:11px;">(@ftype)</span><br/>
            <span style="font-size:12px;">@description</span><br/>
            <i style="color:#666;font-size:11px;">@note</i>
        </div>
    """)
    p.add_tools(hover)

    for track in (fwd, rev):
        if not track["start_angle"]:
            continue
        p.annular_wedge(
            x=0, y=0, inner_radius="inner", outer_radius="outer",
            start_angle="start_angle", end_angle="end_angle", direction="clock",
            color="color", line_color="white", line_width=0.5,
            source=ColumnDataSource(track),
        )

    # Restriction site ticks
    rs_data = {k: [] for k in ("x0","y0","x1","y1","dot_x","dot_y","lx","ly","enzyme","position")}
    for feat in record.features:
        if not feat.qualifiers.get("bound_moiety"):
            continue
        name  = feat.qualifiers["bound_moiety"][0]
        pos   = int(feat.location.start)
        angle = _nt_to_angle(pos, seq_len)
        ca, sa = math.cos(angle), math.sin(angle)
        rs_data["x0"].append(_RS_R_INNER * ca);  rs_data["y0"].append(_RS_R_INNER * sa)
        rs_data["x1"].append(_RS_R_OUTER * ca);  rs_data["y1"].append(_RS_R_OUTER * sa)
        rs_data["dot_x"].append(_RS_R_OUTER * ca); rs_data["dot_y"].append(_RS_R_OUTER * sa)
        rs_data["lx"].append(_RS_R_LABEL * ca);  rs_data["ly"].append(_RS_R_LABEL * sa)
        rs_data["enzyme"].append(name)
        rs_data["position"].append(str(pos + 1))

    if rs_data["x0"]:
        rs_src  = ColumnDataSource(rs_data)
        p.segment(x0="x0", y0="y0", x1="x1", y1="y1",
                  line_color=_RS_COLOR, line_width=1.5, source=rs_src)
        rs_dots = p.circle(x="dot_x", y="dot_y", size=6,
                           color=_RS_COLOR, line_color="white", line_width=0.8,
                           source=rs_src)
        p.add_tools(HoverTool(renderers=[rs_dots],
                              tooltips=[("Enzyme","@enzyme"),("Position","@position bp")]))
        for i in range(len(rs_data["lx"])):
            p.add_layout(Label(
                x=rs_data["lx"][i], y=rs_data["ly"][i],
                text=rs_data["enzyme"][i],
                text_align="left" if rs_data["lx"][i] >= 0 else "right",
                text_baseline="middle", text_font_size="9px",
                text_color=_RS_COLOR, text_font_style="italic",
            ))
    return p


def plot_plasmid_map(record: SeqRecord, out_path: str, figure_size: int = 700) -> None:
    """Render a circular interactive Bokeh map and save as self-contained HTML."""
    plot = build_circular_bokeh_plot(record, figure_size=figure_size)
    output_file(out_path, title=f"Plasmid map – {record.id}")
    bokeh_save(plot)
    log.info("  Interactive map written → %s", out_path)


# ──────────────────────────────────────────────────────────────────────────────
# Annotation class — public API
# ──────────────────────────────────────────────────────────────────────────────

class Annotation:
    """
    Orchestrates plasmid sequence annotation.

    Databases are indexed once at construction time and reused across all
    ``annotate()`` calls, so it is efficient to create a single instance and
    call ``annotate()`` for every assembly in a run.

    Parameters
    ----------
    database : str
        Path to the nucleotide feature database FASTA.
    protein_db : str | None
        Path to the protein database FASTA (optional).
    restriction_enzymes : list[str]
        Enzyme names to map (default: EcoRI, BamHI, HindIII).
    perc_identity : float
        Minimum blastn percent identity (default: 80).
    evalue : float
        Maximum blastn E-value (default: 1e-5).
    word_size : int
        blastn word size (default: 11).
    min_hit_coverage : float
        Minimum fraction of database feature covered by a blastn hit (default: 0.75).
    min_orf_len : int
        Minimum ORF length in nt for CDS detection (default: 100).
    blastp_evalue : float
        Maximum blastp E-value (default: 1e-5).
    blastp_min_pident : float
        Minimum blastp percent identity (default: 30).
    blastp_min_coverage : float
        Minimum blastp subject coverage (default: 0.5).
    threads : int
        BLAST thread count (default: 4).
    figure_size : int
        Bokeh map width/height in pixels (default: 700).
    keep_tmp : bool
        If True, temporary BLAST working directories are not deleted.
    """

    def __init__(
        self,
        database: str,
        protein_db: str | None = None,
        restriction_enzymes: list[str] | None = None,
        perc_identity: float = 80.0,
        evalue: float = 1e-5,
        word_size: int = 11,
        min_hit_coverage: float = 0.75,
        min_orf_len: int = 100,
        blastp_evalue: float = 1e-5,
        blastp_min_pident: float = 30.0,
        blastp_min_coverage: float = 0.5,
        threads: int = 4,
        figure_size: int = 700,
        keep_tmp: bool = False,
    ) -> None:
        check_dependencies(use_blastp=bool(protein_db))

        self.database             = os.path.abspath(database)
        self.perc_identity        = perc_identity
        self.evalue               = evalue
        self.word_size            = word_size
        self.min_hit_coverage     = min_hit_coverage
        self.min_orf_len          = min_orf_len
        self.blastp_evalue        = blastp_evalue
        self.blastp_min_pident    = blastp_min_pident
        self.blastp_min_coverage  = blastp_min_coverage
        self.threads              = threads
        self.figure_size          = figure_size
        self.keep_tmp             = keep_tmp
        self.restriction_enzymes  = restriction_enzymes or ["EcoRI", "BamHI", "HindIII"]

        # Nucleotide database — persists beside the FASTA
        self._nucl_db_prefix = os.path.splitext(self.database)[0]
        make_blast_db(self.database, self._nucl_db_prefix)

        # Load nucleotide feature descriptions once
        self._feature_descriptions: dict[str, str] = {
            rec.id: rec.description
            for rec in SeqIO.parse(self.database, "fasta")
        }

        # Protein database — persists beside the FASTA
        self._protein_db_prefix: str | None = None
        self._protein_descriptions: dict[str, str] = {}
        if protein_db:
            protein_db = os.path.abspath(protein_db)
            self._protein_db_prefix = os.path.splitext(protein_db)[0]
            make_blast_protein_db(protein_db, self._protein_db_prefix)
            self._protein_descriptions = {
                rec.id: rec.description
                for rec in SeqIO.parse(protein_db, "fasta")
            }

    # ── public method ──────────────────────────────────────────────────────────

    def annotate(
        self,
        query_fasta: str,
        outdir: str,
    ) -> list[SeqRecord]:
        """
        Annotate all sequences in *query_fasta* and write results to *outdir*.

        For each sequence the following files are written:
          - ``<id>.gbk``       — annotated GenBank file
          - ``<id>_map.html``  — interactive circular Bokeh map

        A combined ``annotation_summary.tsv`` is written once at the end.

        Parameters
        ----------
        query_fasta : str
            Path to the plasmid sequences in FASTA format (assumed circular).
        outdir : str
            Output directory (created if absent).

        Returns
        -------
        list[SeqRecord]
            Annotated SeqRecord objects (one per input sequence).
        """
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        records = list(SeqIO.parse(query_fasta, "fasta"))
        if not records:
            log.warning("No sequences found in %s", query_fasta)
            return []
        log.info("Annotating %d sequence(s) from %s …", len(records), query_fasta)

        workdir = tempfile.mkdtemp(prefix="plasmid_annotator_")
        log.info("Temporary working directory: %s", workdir)

        try:
            # Write doubled query for circular BLAST
            doubled_fasta = os.path.join(workdir, "query_doubled.fasta")
            write_doubled_fasta(records, doubled_fasta)

            # blastn
            blast_out = os.path.join(workdir, "blast_hits.tsv")
            run_blastn(
                doubled_fasta, self._nucl_db_prefix, blast_out,
                perc_identity=self.perc_identity, evalue=self.evalue,
                word_size=self.word_size, threads=self.threads,
            )
            all_hits = parse_blast_table(blast_out)
            log.info("Parsed %d total blastn hits.", len(all_hits))

            annotated_records: list[SeqRecord] = []
            for rec in records:
                log.info("Annotating %s (len=%d bp) …", rec.id, len(rec.seq))

                ann_rec = annotate_record(
                    rec, all_hits,
                    min_hit_coverage=self.min_hit_coverage,
                    min_orf_len=self.min_orf_len,
                    feature_descriptions=self._feature_descriptions,
                    protein_db_prefix=self._protein_db_prefix,
                    workdir=workdir,
                    blastp_evalue=self.blastp_evalue,
                    blastp_min_pident=self.blastp_min_pident,
                    blastp_min_coverage=self.blastp_min_coverage,
                    threads=self.threads,
                    protein_descriptions=self._protein_descriptions,
                )

                # Restriction sites
                if self.restriction_enzymes:
                    rs_feats = find_restriction_sites(ann_rec, self.restriction_enzymes)
                    ann_rec.features.extend(rs_feats)
                    log.info("  %d restriction site(s) added.", len(rs_feats))

                # GenBank output
                gbk_path = outdir_path / f"{rec.id}.gbk"
                with open(gbk_path, "w") as fh:
                    SeqIO.write(ann_rec, fh, "genbank")
                log.info("  GenBank written → %s", gbk_path)

                # Interactive map
                map_path = outdir_path / f"{rec.id}_map.html"
                plot_plasmid_map(ann_rec, str(map_path), figure_size=self.figure_size)

                annotated_records.append(ann_rec)

            # Summary TSV
            write_summary(annotated_records, str(outdir_path / "annotation_summary.tsv"))

        finally:
            if self.keep_tmp:
                log.info("Keeping temporary directory: %s", workdir)
            else:
                shutil.rmtree(workdir, ignore_errors=True)

        log.info("Annotation complete. Results in %s", outdir)
        return annotated_records
