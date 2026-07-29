#!/usr/bin/env python3
"""
Soybean (or any plant) low-coverage bin-map pipeline → QTL IciMapping .bip

One-script rewrite of the classic rice binmap workflow (log.txt):
  parent gVCFs → joint genotype → pseudo genomes → dual-ref BWA →
  discriminative reads → per-sample bins/plots → consensus physical bins +
  cM map → .bip (with/without phenotypes)

Example:
  python3 binmap_pipeline.py \\
    --p1-gvcf parent1.g.vcf.gz --p2-gvcf parent2.g.vcf.gz \\
    --fq-list examples/fq.list --pheno-list examples/pheno.list \\
    --ref ZH13.v2.fasta --outdir binmap_out
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:  # pragma: no cover
    Image = ImageDraw = ImageFont = None  # type: ignore


LOG = logging.getLogger("binmap")

POP_TYPE_CODE = {
    "P1BC1F1": 1,
    "P2BC1F1": 2,
    "DH": 3,
    "F1DH": 3,
    "RIL": 4,
    "F1RIL": 4,
    "F2": 7,
    "F3": 8,
}

DEFAULT_BWA = "bwa"
DEFAULT_SEQTK = "seqtk"
DEFAULT_GATK = "gatk"
DEFAULT_REF = ""  # require --ref

CHR_RE = re.compile(r"(?:ZH13\.)?Chr0*(\d+)$", re.I)
MD_SNP_RE = re.compile(r"\^?([A-Za-z]+)|(\d+)")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def run(cmd: Sequence[str], cwd: Optional[Path] = None, log_path: Optional[Path] = None) -> None:
    LOG.info("RUN: %s", " ".join(map(str, cmd)))
    if log_path:
        with open(log_path, "ab") as lf:
            lf.write(("\n$ " + " ".join(map(str, cmd)) + "\n").encode())
            p = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=lf, stderr=subprocess.STDOUT)
    else:
        p = subprocess.run(cmd, cwd=str(cwd) if cwd else None)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(map(str, cmd))}")


def which_or(path: str, fallbacks: Sequence[str] = ()) -> str:
    if path and Path(path).exists():
        return path
    w = shutil.which(path) if path else None
    if w:
        return w
    for fb in fallbacks:
        if fb and Path(fb).exists():
            return fb
        w2 = shutil.which(fb) if fb else None
        if w2:
            return w2
    raise FileNotFoundError(f"Cannot find executable: {path}")


def normalize_chrom(name: str) -> Optional[str]:
    """ZH13.Chr01 / Chr01 / chromosome01 → Chr01; non-nuclear → None."""
    name = name.strip()
    m = CHR_RE.search(name.replace("chromosome", "Chr").replace("Chromosome", "Chr"))
    if not m:
        # try chromosome01 style
        m2 = re.search(r"(?:chr|chromosome)0*(\d+)$", name, re.I)
        if not m2:
            return None
        n = int(m2.group(1))
    else:
        n = int(m.group(1))
    if n < 1 or n > 20:
        return None
    return f"Chr{n:02d}"


def chrom_sort_key(c: str) -> int:
    m = re.search(r"(\d+)$", c)
    return int(m.group(1)) if m else 999


def read_fai(ref: Path) -> Dict[str, int]:
    fai = Path(str(ref) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(f"Missing fasta index: {fai}")
    lengths: Dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            chrom, length = line.split("\t")[:2]
            nc = normalize_chrom(chrom)
            if nc:
                lengths[nc] = int(length)
    if not lengths:
        raise RuntimeError("No nuclear Chr01–Chr20 found in reference .fai")
    return lengths


def ensure_ref_dict(ref: Path, gatk: str) -> None:
    """GATK accepts basename.dict for basename.fasta (e.g. ZH13.v2.dict)."""
    candidates = [
        ref.with_suffix(".dict"),  # ZH13.v2.fasta → ZH13.v2.dict
        Path(str(ref) + ".dict"),  # ZH13.v2.fasta.dict
        ref.parent / (ref.name.replace(".fasta", ".dict").replace(".fa", ".dict")),
    ]
    for c in candidates:
        if c.exists():
            LOG.info("Using sequence dictionary: %s", c)
            return
    out = ref.with_suffix(".dict")
    run([gatk, "CreateSequenceDictionary", "-R", str(ref), "-O", str(out)])


def parse_gt(gt_field: str) -> Optional[Tuple[int, int]]:
    gt = gt_field.split(":")[0]
    if gt in (".", "./.", ".|."):
        return None
    gt = gt.replace("|", "/")
    a, b = gt.split("/")
    if a == "." or b == ".":
        return None
    return int(a), int(b)


def md_snp_count(md: str) -> Optional[int]:
    """Count substitution mismatches in MD:Z tag; None if deletion present."""
    if md.startswith("MD:Z:"):
        md = md[5:]
    n = 0
    i = 0
    while i < len(md):
        if md[i].isdigit():
            while i < len(md) and md[i].isdigit():
                i += 1
        elif md[i] == "^":
            return None  # indel / deletion in MD
        elif md[i].isalpha():
            n += 1
            i += 1
        else:
            i += 1
    return n


def cigar_has_indel(cigar: str) -> bool:
    return bool(re.search(r"[ID]", cigar))


def fq_pair_from_r1(r1: Path) -> Path:
    name = r1.name
    for a, b in [
        ("_1.clean.fq.gz", "_2.clean.fq.gz"),
        ("_1.clean.fastq.gz", "_2.clean.fastq.gz"),
        ("_1.fq.gz", "_2.fq.gz"),
        ("_1.fastq.gz", "_2.fastq.gz"),
        ("_1.fq", "_2.fq"),
        ("_1.fastq", "_2.fastq"),
        ("_R1.fastq.gz", "_R2.fastq.gz"),
        ("_R1.fq.gz", "_R2.fq.gz"),
        ("_1.clean.fq", "_2.clean.fq"),
    ]:
        if name.endswith(a):
            return r1.with_name(name[: -len(a)] + b)
    # generic _1 → _2 before extension
    m = re.sub(r"(_1)(\.|$)", r"_2\2", name, count=1)
    if m != name:
        return r1.with_name(m)
    raise ValueError(f"Cannot derive R2 from R1 path: {r1}")


def sample_id_from_r1(r1: Path) -> str:
    name = r1.name
    name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", name, flags=re.I)
    name = re.sub(r"(_1|_R1)(\.clean)?$", "", name)
    name = re.sub(r"\.clean$", "", name)
    return name


def peek_read_name_style(fq: Path, n: int = 4) -> str:
    """Return 'slash' if /1 /2 present, else 'plain'."""
    with open_text(fq) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                if line.rstrip().endswith("/1") or line.rstrip().endswith("/2"):
                    return "slash"
                if " 1:" in line or " 2:" in line:
                    return "illumina"
            if i > n * 4:
                break
    return "plain"


def peek_read_len(fq: Path) -> int:
    with open_text(fq) as fh:
        for i, line in enumerate(fh):
            if i == 1:
                return len(line.rstrip("\n\r"))
    return 100


# ---------------------------------------------------------------------------
# Stage 1: parents → SNP table → pseudo genomes
# ---------------------------------------------------------------------------

def stage_genotype(
    p1_gvcf: Path,
    p2_gvcf: Path,
    ref: Path,
    work: Path,
    gatk: str,
    java_mem: str,
    bcftools: str = "bcftools",
    samtools: str = "samtools",
    parent_vcf: Optional[Path] = None,
    p1_name: Optional[str] = None,
    p2_name: Optional[str] = None,
) -> Path:
    work.mkdir(parents=True, exist_ok=True)
    ensure_ref_dict(ref, gatk)
    if not Path(str(ref) + ".fai").exists():
        run([samtools, "faidx", str(ref)])

    final_vcf = work / "parents.SNP.PASS.vcf.gz"
    logf = work / "gatk.log"

    if final_vcf.exists():
        LOG.info("Reuse existing %s", final_vcf)
        return final_vcf

    # Fast path: subset an existing joint/cohort VCF to the two parents
    if parent_vcf is not None:
        if not p1_name or not p2_name:
            raise ValueError("--parent-vcf requires --p1-name and --p2-name (VCF sample IDs)")
        LOG.info("Extracting parents %s,%s from %s (per-chromosome parallel)", p1_name, p2_name, parent_vcf)
        parts = work / "vcf_parts"
        parts.mkdir(exist_ok=True)
        procs = []
        part_files = []
        for i in range(1, 21):
            chrom = f"ZH13.Chr{i:02d}"
            part = parts / f"{chrom}.vcf.gz"
            part_files.append(part)
            cmd = [
                bcftools,
                "view",
                "-r",
                chrom,
                "-s",
                f"{p1_name},{p2_name}",
                "-v",
                "snps",
                "-f",
                "PASS",
                "-O",
                "z",
                "-o",
                str(part),
                str(parent_vcf),
            ]
            LOG.info("RUN: %s", " ".join(cmd))
            procs.append(subprocess.Popen(cmd, stdout=open(logf, "ab"), stderr=subprocess.STDOUT))
        rc = 0
        for p in procs:
            rc = max(rc, p.wait())
        if rc != 0:
            # fallback: contig names may be Chr01 without ZH13. prefix
            LOG.warning("ZH13.Chr* extract failed (rc=%s); retrying Chr*", rc)
            procs = []
            part_files = []
            for i in range(1, 21):
                chrom = f"Chr{i:02d}"
                part = parts / f"{chrom}.vcf.gz"
                part_files.append(part)
                cmd = [
                    bcftools, "view", "-r", chrom, "-s", f"{p1_name},{p2_name}",
                    "-v", "snps", "-f", "PASS", "-O", "z", "-o", str(part), str(parent_vcf),
                ]
                procs.append(subprocess.Popen(cmd, stdout=open(logf, "ab"), stderr=subprocess.STDOUT))
            rc = 0
            for p in procs:
                rc = max(rc, p.wait())
            if rc != 0:
                raise RuntimeError("Failed to extract parents from --parent-vcf")
        existing = [str(p) for p in part_files if p.exists() and p.stat().st_size > 0]
        if not existing:
            raise RuntimeError("No chromosome parts produced from --parent-vcf")
        run([bcftools, "concat", "-O", "z", "-o", str(final_vcf)] + existing, log_path=logf)
        run([bcftools, "index", "-f", str(final_vcf)], log_path=logf)
        return final_vcf

    tmp_g = work / "parents.combined.g.vcf.gz"
    total = work / "parents.total.vcf.gz"
    snp = work / "parents.SNP.vcf.gz"
    filt = work / "parents.SNP.filter.vcf.gz"
    gatk_cmd = [gatk, "--java-options", f"-Xmx{java_mem}"]

    run(
        gatk_cmd
        + [
            "CombineGVCFs",
            "-R",
            str(ref),
            "--variant",
            str(p1_gvcf),
            "--variant",
            str(p2_gvcf),
            "-O",
            str(tmp_g),
        ],
        log_path=logf,
    )
    run(
        gatk_cmd + ["GenotypeGVCFs", "-R", str(ref), "-V", str(tmp_g), "-O", str(total)],
        log_path=logf,
    )
    run(
        gatk_cmd + ["SelectVariants", "-V", str(total), "--select-type-to-include", "SNP", "-O", str(snp)],
        log_path=logf,
    )
    run(
        gatk_cmd
        + [
            "VariantFiltration",
            "-R",
            str(ref),
            "-V",
            str(snp),
            "--cluster-size",
            "3",
            "--cluster-window-size",
            "10",
            "--filter-expression",
            "QD < 10.0",
            "--filter-name",
            "lowQD",
            "--filter-expression",
            "FS > 15.0",
            "--filter-name",
            "highFS",
            "--genotype-filter-expression",
            "DP < 5 || DP > 200",
            "--genotype-filter-name",
            "InvalidDP",
            "-O",
            str(filt),
        ],
        log_path=logf,
    )
    run(
        [bcftools, "view", "-f", "PASS", "-O", "z", "-o", str(final_vcf), str(filt)],
        log_path=logf,
    )
    run([bcftools, "index", "-f", str(final_vcf)], log_path=logf)
    return final_vcf


def extract_parent_snps(vcf: Path, out_snp: Path) -> int:
    """
    Write TSV: chrom pos ref alt p1_base p2_base
    Only biallelic SNPs where both parents are homozygous and differ.
    """
    n = 0
    with open_text(vcf) as fh, open(out_snp, "w") as out:
        samples: List[str] = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                if len(samples) < 2:
                    raise RuntimeError("VCF needs two parent samples")
                continue
            cols = line.rstrip("\n").split("\t")
            chrom = normalize_chrom(cols[0])
            if not chrom:
                continue
            filt = cols[6]
            if filt not in ("PASS", "."):
                continue
            ref, alt = cols[3], cols[4]
            if len(ref) != 1 or len(alt) != 1 or alt == "*":
                continue
            fmt = cols[8].split(":")
            if "GT" not in fmt:
                continue

            def _gt(sample_field: str) -> Optional[Tuple[int, int]]:
                parts = sample_field.split(":")
                d = {k: parts[i] if i < len(parts) else "." for i, k in enumerate(fmt)}
                ft = d.get("FT", "PASS")
                if ft not in ("PASS", ".", ""):
                    return None
                return parse_gt(d.get("GT", "./."))

            g1 = _gt(cols[9])
            g2 = _gt(cols[10]) if len(cols) > 10 else None
            if g1 is None or g2 is None:
                continue
            a1, b1 = g1
            a2, b2 = g2
            if a1 != b1 or a2 != b2:
                continue  # skip het
            if a1 == a2:
                continue  # same genotype
            base = {0: ref, 1: alt}
            if a1 not in base or a2 not in base:
                continue
            p1b, p2b = base[a1], base[a2]
            out.write(f"{chrom}\t{cols[1]}\t{ref}\t{alt}\t{p1b}\t{p2b}\n")
            n += 1
    LOG.info("Differential homozygous SNPs: %d → %s", n, out_snp)
    return n


def build_pseudo_genome(ref: Path, snp_table: Path, which: str, out_fa: Path) -> None:
    """which: 'p1' or 'p2' → TSV columns p1_base / p2_base."""
    col = 4 if which == "p1" else 5
    subs: Dict[str, Dict[int, str]] = {}
    with open(snp_table) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            chrom, pos, base = p[0], int(p[1]), p[col]
            subs.setdefault(chrom, {})[pos] = base.upper()

    LOG.info("Building pseudo genome %s (%d chroms with SNPs) → %s", which, len(subs), out_fa)
    with open(ref) as fin, open(out_fa, "w") as fout:
        chrom: Optional[str] = None
        seq_chunks: List[str] = []

        def flush() -> None:
            nonlocal chrom, seq_chunks
            if chrom is None:
                return
            nc = normalize_chrom(chrom)
            raw = "".join(seq_chunks)
            if nc and nc in subs:
                arr = list(raw)
                for p, b in subs[nc].items():
                    if 1 <= p <= len(arr):
                        arr[p - 1] = b
                raw = "".join(arr)
            fout.write(f">{chrom}\n")
            for i in range(0, len(raw), 80):
                fout.write(raw[i : i + 80] + "\n")
            seq_chunks = []

        for line in fin:
            if line.startswith(">"):
                flush()
                chrom = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        flush()


def bwa_index(fa: Path, bwa: str) -> None:
    if Path(str(fa) + ".bwt").exists():
        LOG.info("bwa index exists: %s", fa)
        return
    run([bwa, "index", str(fa)])


# ---------------------------------------------------------------------------
# Stage 2: per-sample
# ---------------------------------------------------------------------------

@dataclass
class BinRecord:
    chrom: str
    start: int
    end: int
    origin: str  # parent1 / parent2 / hetero


@dataclass
class SampleResult:
    sample: str
    bins: List[BinRecord] = field(default_factory=list)
    n_reads: int = 0
    ok: bool = True
    error: str = ""


def subsample_and_trim(
    r1: Path,
    r2: Path,
    out_r1: Path,
    out_r2: Path,
    n_pairs: int,
    read_len: int,
    seqtk: str,
    seed: int = 42,
) -> None:
    """Subsample n_pairs then optionally hard-trim to read_len (0 = no trim)."""
    tmpdir = out_r1.parent
    s1 = tmpdir / f"{out_r1.stem}.sub.fq"
    s2 = tmpdir / f"{out_r2.stem}.sub.fq"

    with open(s1, "wb") as out:
        subprocess.run(
            [seqtk, "sample", "-s", str(seed), str(r1), str(n_pairs)],
            stdout=out,
            check=True,
        )
    with open(s2, "wb") as out:
        subprocess.run(
            [seqtk, "sample", "-s", str(seed), str(r2), str(n_pairs)],
            stdout=out,
            check=True,
        )

    if read_len and read_len > 0:
        _hard_trim_fq(s1, out_r1, read_len)
        _hard_trim_fq(s2, out_r2, read_len)
        s1.unlink(missing_ok=True)
        s2.unlink(missing_ok=True)
    else:
        # recompress to .fq.gz outputs
        _recompress_fq(s1, out_r1)
        _recompress_fq(s2, out_r2)
        s1.unlink(missing_ok=True)
        s2.unlink(missing_ok=True)


def _recompress_fq(inp: Path, outp: Path) -> None:
    with open_text(inp) as fin, gzip.open(outp, "wt") as fout:
        shutil.copyfileobj(fin, fout)


def _hard_trim_fq(inp: Path, outp: Path, read_len: int) -> None:
    with open_text(inp) as fin, gzip.open(outp, "wt") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            seq = fin.readline().rstrip("\n\r")
            plus = fin.readline()
            qual = fin.readline().rstrip("\n\r")
            fout.write(h if h.endswith("\n") else h + "\n")
            fout.write(seq[:read_len] + "\n")
            fout.write(plus if plus.endswith("\n") else plus + "\n")
            fout.write(qual[:read_len] + "\n")


def maybe_tag_headers(r1: Path, r2: Path, out_r1: Path, out_r2: Path, style: str) -> Tuple[Path, Path, bool]:
    """
    If read names cannot distinguish mates, append A/B comment for bwa -C.
    Returns (r1_path, r2_path, use_C).
    """
    if style in ("slash", "illumina"):
        return r1, r2, False
    for inp, outp, tag in ((r1, out_r1, "A"), (r2, out_r2, "B")):
        with open_text(inp) as fin, gzip.open(outp, "wt") as fout:
            i = 0
            for line in fin:
                if i % 4 == 0:
                    line = line.rstrip("\n\r")
                    if not line.endswith(f" {tag}"):
                        line = line + f" {tag}"
                    fout.write(line + "\n")
                else:
                    fout.write(line if line.endswith("\n") else line + "\n")
                i += 1
    return out_r1, out_r2, True


def align_to_pseudo(r1: Path, r2: Path, index: Path, out_sam: Path, bwa: str, threads: int, use_C: bool) -> None:
    cmd = [bwa, "mem", "-t", str(threads)]
    if use_C:
        cmd.append("-C")
    cmd += [str(index), str(r1), str(r2)]
    LOG.info("RUN: %s > %s", " ".join(cmd), out_sam)
    with open(out_sam, "w") as out:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise RuntimeError(f"bwa mem failed: {p.stderr.decode()[:500]}")


def parse_sam_md(sam: Path) -> Dict[str, Tuple[str, int, int, str, int]]:
    """
    key → (chrom_norm, start, end, md, n_snp)
    key = QNAME (unique per mate when /1 /2 present)
    """
    out: Dict[str, Tuple[str, int, int, str, int]] = {}
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            qname, flag, rname, pos, mapq, cigar = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
            if rname == "*" or mapq <= 0:
                continue
            if cigar_has_indel(cigar):
                continue
            chrom = normalize_chrom(rname)
            if not chrom:
                continue
            md = None
            for t in f[11:]:
                if t.startswith("MD:Z:"):
                    md = t
                    break
            if not md:
                continue
            n_snp = md_snp_count(md)
            if n_snp is None:
                continue
            # read length from SEQ or CIGAR M
            seq = f[9]
            rlen = len(seq) if seq != "*" else sum(int(x) for x in re.findall(r"(\d+)M", cigar))
            end = pos + rlen - 1
            out[qname] = (chrom, pos, end, md, n_snp)
    return out


def discriminative_reads(
    sam1: Path, sam2: Path, out_rlt: Path
) -> List[Tuple[str, str, int, int]]:
    """
    Classic rule (generalized to any L):
      P1 if (p1 perfect & p2 1snp) or (p1 1snp & p2 2snp)
      P2 if (p2 perfect & p1 1snp) or (p2 1snp & p1 2snp)
    Returns list of (origin P1/P2, chrom, start, end)
    """
    a = parse_sam_md(sam1)
    b = parse_sam_md(sam2)
    rows: List[Tuple[str, str, int, int]] = []
    with open(out_rlt, "w") as out:
        for qid, (c1, s1, e1, md1, n1) in a.items():
            if qid not in b:
                continue
            c2, s2, e2, md2, n2 = b[qid]
            if c1 != c2:
                continue
            # prefer coordinates from the better match later; use midpoint-ish start
            origin = None
            if (n1 == 0 and n2 == 1) or (n1 == 1 and n2 == 2):
                origin = "P1"
            elif (n2 == 0 and n1 == 1) or (n2 == 1 and n1 == 2):
                origin = "P2"
            if origin is None:
                continue
            start, end = s1, e1
            out.write(f"{origin}\t{qid}\t{c1}\t{start}\t{end}\t{md1}\t{md2}\n")
            rows.append((origin, c1, start, end))
    rows.sort(key=lambda x: (chrom_sort_key(x[1]), x[2]))
    return rows


def chrom_to_perl(chrom: str) -> str:
    """Chr01 → chromosome01 (Seq2Bin_F24.pl naming)."""
    n = chrom_sort_key(chrom)
    return f"chromosome{n:02d}"


def chrom_from_perl(name: str) -> str:
    m = re.search(r"(\d+)$", name)
    if not m:
        return name
    return f"Chr{int(m.group(1)):02d}"


def write_chrom_length_list(chrom_len: Dict[str, int], path: Path) -> None:
    with open(path, "w") as out:
        for chrom in sorted(chrom_len, key=chrom_sort_key):
            out.write(f"{chrom_to_perl(chrom)}\t{chrom_len[chrom]}\n")


def write_perl_rlt(reads: List[Tuple[str, str, int, int]], path: Path) -> None:
    """Write Seq2Bin-compatible filter.rlt (chromosome01 naming, tab-separated)."""
    with open(path, "w") as out:
        for origin, chrom, start, end in sorted(reads, key=lambda x: (chrom_sort_key(x[1]), x[2])):
            out.write(f"{origin}\tx\t{chrom_to_perl(chrom)}\t{start}\t{end}\n")


def seq2bin(
    reads: List[Tuple[str, str, int, int]],
    chrom_len: Dict[str, int],
    win_frac: float = 0.8,
    min_bin: int = 300_000,
    work_dir: Optional[Path] = None,
    p1_name: str = "P1",
    p2_name: str = "P2",
) -> List[BinRecord]:
    """
    Call original Seq2Bin_F24.pl for bin calling (edge / hetero logic identical).
    Falls back to simplified Python only if Perl path is missing.
    """
    script = Path(__file__).resolve().parent / "Seq2Bin_F24.pl"
    if not script.exists() or not reads:
        LOG.warning("Seq2Bin_F24.pl missing or no reads; using simplified seq2bin")
        return seq2bin_simple(reads, chrom_len, win_frac=win_frac, min_bin=min_bin)

    wd = Path(work_dir) if work_dir else Path(tempfile.mkdtemp(prefix="seq2bin_"))
    wd.mkdir(parents=True, exist_ok=True)
    chrom_list = wd / "chrom_length.list"
    rlt = wd / "sample.filter.rlt"
    write_chrom_length_list(chrom_len, chrom_list)
    write_perl_rlt(reads, rlt)

    # Seq2Bin must be run with basename args under cwd=wd; relative outdir paths break open().
    # Drop stale products so a failed run cannot reuse an old .bin.
    for pat in ("*.bin", "*.edge", "*win*.edge"):
        for pth in wd.glob(pat):
            pth.unlink(missing_ok=True)

    logf = wd / "seq2bin.log"
    if logf.exists():
        logf.unlink()
    cmd = [
        "perl",
        str(script),
        rlt.name,
        chrom_list.name,
        str(win_frac),
        p1_name,
        p2_name,
        str(int(min_bin)),
    ]
    LOG.info("RUN (cwd=%s): %s", wd, " ".join(cmd))
    with open(logf, "ab") as lf:
        p = subprocess.run(cmd, cwd=str(wd), stdout=lf, stderr=subprocess.STDOUT)

    # Seq2Bin writes sample.filter.bin then may die on GD/font at plot stage
    bin_candidates = [
        wd / "sample.filter.bin",
        wd / "sample.bin",
    ]
    bin_path = next((c for c in bin_candidates if c.exists() and c.stat().st_size > 0), None)
    if bin_path is None:
        found = list(wd.glob("*.bin"))
        bin_path = found[0] if found else None
    if bin_path is None:
        LOG.error("Seq2Bin failed (rc=%s); see %s; falling back to simplified seq2bin", p.returncode, logf)
        return seq2bin_simple(reads, chrom_len, win_frac=win_frac, min_bin=min_bin)

    bins: List[BinRecord] = []
    with open(bin_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.split()
            if len(cols) < 4:
                continue
            chrom = chrom_from_perl(cols[0])
            if chrom not in chrom_len:
                continue
            origin = cols[3]
            if origin not in ("parent1", "parent2", "hetero"):
                continue
            bins.append(BinRecord(chrom, int(cols[1]), int(cols[2]), origin))
    if not bins:
        LOG.error("Seq2Bin produced empty bins (rc=%s); falling back to simplified seq2bin", p.returncode)
        return seq2bin_simple(reads, chrom_len, win_frac=win_frac, min_bin=min_bin)
    LOG.info("Seq2Bin_F24.pl → %d bins from %s (perl_rc=%s)", len(bins), bin_path.name, p.returncode)
    return bins


def seq2bin_simple(
    reads: List[Tuple[str, str, int, int]],
    chrom_len: Dict[str, int],
    win_frac: float = 0.8,
    min_bin: int = 300_000,
) -> List[BinRecord]:
    """Simplified sliding-window bin caller (fallback only)."""
    by_chr: Dict[str, List[Tuple[str, int, int]]] = {}
    for origin, chrom, start, end in reads:
        by_chr.setdefault(chrom, []).append((origin, start, end))

    all_bins: List[BinRecord] = []
    for chrom in sorted(by_chr, key=chrom_sort_key):
        arr = sorted(by_chr[chrom], key=lambda x: x[1])
        n = len(arr)
        if n < 10:
            continue
        if n < 5000:
            win = 15
        elif n < 10000:
            win = 25
        elif n < 20000:
            win = 39
        elif n < 100000:
            win = 59
        else:
            win = 99
        half = win // 2
        dom = int(win * win_frac)
        calls: List[str] = ["hetero"] * n
        for i in range(n):
            left = max(0, i - half)
            right = min(n, i + half + 1)
            w = arr[left:right]
            p1 = sum(1 for o, _, _ in w if o == "P1")
            p2 = len(w) - p1
            if p1 > dom:
                calls[i] = "parent1"
            elif p2 > dom:
                calls[i] = "parent2"
            else:
                calls[i] = "hetero"

        bins: List[BinRecord] = []
        i = 0
        clen = chrom_len.get(chrom, arr[-1][2])
        while i < n:
            j = i
            while j + 1 < n and calls[j + 1] == calls[i]:
                j += 1
            start = 1 if i == 0 else arr[i][1]
            end = clen if j == n - 1 else arr[j][2]
            bins.append(BinRecord(chrom, start, end, calls[i]))
            i = j + 1

        for _ in range(100):
            if len(bins) <= 1:
                break
            changed = False
            new: List[BinRecord] = []
            k = 0
            while k < len(bins):
                if k < len(bins) - 1 and (bins[k].end - bins[k].start + 1) < min_bin:
                    if k > 0 and new and new[-1].origin == bins[k + 1].origin:
                        prev = new.pop()
                        nxt = bins[k + 1]
                        new.append(BinRecord(chrom, prev.start, nxt.end, prev.origin))
                        k += 2
                        changed = True
                        continue
                    if new:
                        prev = new.pop()
                        new.append(BinRecord(chrom, prev.start, bins[k].end, prev.origin))
                        k += 1
                        changed = True
                        continue
                new.append(bins[k])
                k += 1
            bins = new
            if not changed:
                break

        if bins:
            bins[0] = BinRecord(chrom, 1, bins[0].end, bins[0].origin)
            bins[-1] = BinRecord(chrom, bins[-1].start, clen, bins[-1].origin)
        all_bins.extend(bins)
    return all_bins


def find_arial_font(size: int = 16):
    """Prefer Liberation Sans (Arial-metric compatible), then DejaVu Sans."""
    if ImageFont is None:
        return None
    candidates = [
        "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
        "/usr/share/fonts/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            try:
                return ImageFont.truetype(path, size=size)
            except OSError:
                continue
    try:
        return ImageFont.load_default()
    except Exception:
        return None


def plot_bins(
    bins: List[BinRecord],
    chrom_len: Dict[str, int],
    reads: List[Tuple[str, str, int, int]],
    out_png: Path,
    p1_name: str,
    p2_name: str,
    sample: str,
) -> None:
    """
    Seq2Bin_F24-style figure.
    Right panel: sample / parents / proportions above color legend.
    Scale: 10 kb / pixel.
    """
    if Image is None:
        LOG.warning("PIL not available; skip plot %s", out_png)
        return

    bp_per_px = 10_000
    left = 280
    top = 160
    row_h = 120
    right_panel = 780
    chroms = sorted(chrom_len, key=chrom_sort_key)
    max_bp = max(chrom_len.values()) if chrom_len else 1
    track_w = int(max_bp / bp_per_px + 0.5)
    width = max(left + track_w + right_panel, 4200)
    height = top + len(chroms) * row_h + 50

    im = Image.new("RGB", (width, height), "white")
    dr = ImageDraw.Draw(im)
    font_chr = find_arial_font(52)
    font_mb = find_arial_font(52)
    font_leg = find_arial_font(48)

    red = (220, 20, 20)
    blue = (20, 60, 220)
    gold = (230, 180, 20)
    black = (0, 0, 0)
    gray = (200, 200, 200)

    def txt(xy, s, fill=black, fnt=None):
        dr.text(xy, s, fill=fill, font=fnt or font_leg)

    tot = sum(max(0, b.end - b.start + 1) for b in bins) or 1
    bp1 = sum(b.end - b.start + 1 for b in bins if b.origin == "parent1")
    bp2 = sum(b.end - b.start + 1 for b in bins if b.origin == "parent2")
    bph = sum(b.end - b.start + 1 for b in bins if b.origin == "hetero")

    # Right panel shifted left so long labels stay inside the image
    lx = left + track_w - 100
    info_lines = [
        (sample, black),
        (f"parent1={p1_name}", red),
        (f"parent2={p2_name}", blue),
        (
            f"P1={100 * bp1 / tot:.1f}%  P2={100 * bp2 / tot:.1f}%  hetero={100 * bph / tot:.1f}%",
            black,
        ),
    ]
    legend_items = [(p1_name, red), (p2_name, blue), ("hetero", gold)]
    leg_box = 52
    leg_gap = 84
    info_gap = 58
    block_h = len(info_lines) * info_gap + 24 + len(legend_items) * leg_gap
    y0_block = max(top, (height - block_h) // 2)
    # white backdrop so labels remain readable if they overlap track ends
    dr.rectangle(
        [lx - 16, y0_block - 12, width - 12, y0_block + block_h + 12],
        fill=(255, 255, 255),
        outline=None,
    )
    y = y0_block
    for lab, col in info_lines:
        txt((lx, y), lab, fill=col, fnt=font_leg)
        y += info_gap
    y += 16
    for lab, col in legend_items:
        dr.rectangle([lx, y, lx + leg_box, y + leg_box], fill=col, outline=black, width=2)
        txt((lx + leg_box + 16, y + 4), lab, fill=col, fnt=font_leg)
        y += leg_gap

    # Mb scale (font same size as Chr)
    scale_y0, scale_y1 = 120, 134
    dr.rectangle([left, scale_y0, left + track_w, scale_y1], fill=black)
    max_mb = int(max_bp / 1_000_000) + 1
    for mb in range(0, max_mb + 1):
        x = left + int(mb * 1_000_000 / bp_per_px)
        if x > left + track_w:
            break
        tick_h = 20 if mb % 5 == 0 else 10
        dr.line([x, scale_y0 - tick_h, x, scale_y1], fill=black, width=2)
        if mb % 5 == 0:
            txt((x - 36, 48), f"{mb} Mb", fnt=font_mb)

    def x_of(pos: int) -> int:
        return left + int(pos / bp_per_px)

    for chrom in chroms:
        n = chrom_sort_key(chrom)
        y0 = top + (n - 1) * row_h
        clen = chrom_len[chrom]
        x_end = x_of(clen)
        dr.rectangle([left, y0, x_end, y0 + 100], outline=gray, fill=gray)
        txt((6, y0 + 22), chrom, fnt=font_chr)

        for b in bins:
            if b.chrom != chrom:
                continue
            x0, x1 = x_of(b.start), x_of(b.end)
            if x1 <= x0:
                x1 = x0 + 1
            col = red if b.origin == "parent1" else blue if b.origin == "parent2" else gold
            dr.rectangle([x0, y0 + 21, x1, y0 + 40], fill=col, outline=col)

        for origin, c, start, end in reads:
            if c != chrom:
                continue
            x0, x1 = x_of(start), x_of(end)
            if x1 <= x0:
                x1 = x0 + 1
            if origin == "P1":
                dr.rectangle([x0, y0 + 51, x1, y0 + 70], fill=red, outline=red)
            elif origin == "P2":
                dr.rectangle([x0, y0 + 71, x1, y0 + 90], fill=blue, outline=blue)

    im.save(out_png)


def write_sample_bin(
    bins: List[BinRecord],
    path: Path,
    sample: str,
    p1_name: str,
    p2_name: str,
) -> None:
    with open(path, "w") as out:
        out.write(f"# parent1={p1_name}\tparent2={p2_name}\tsample={sample}\n")
        out.write("# chrom\tstart\tend\torigin\tsample\tlength\n")
        for b in bins:
            length = b.end - b.start + 1
            out.write(f"{b.chrom}\t{b.start}\t{b.end}\t{b.origin}\t{sample}\t{length}\n")


def load_bin_file(path: Path) -> List[BinRecord]:
    bins: List[BinRecord] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            bins.append(BinRecord(p[0], int(p[1]), int(p[2]), p[3]))
    return bins


def process_one_sample(args: Tuple) -> SampleResult:
    (
        sample,
        r1,
        r2,
        p1_fa,
        p2_fa,
        work_root,
        chrom_len,
        n_pairs,
        read_len,
        bwa,
        seqtk,
        threads,
        p1_name,
        p2_name,
        win_frac,
        min_bin,
        keep_align,
        force,
    ) = args
    sdir = Path(work_root) / "samples" / sample
    sdir.mkdir(parents=True, exist_ok=True)
    bin_path = Path(work_root) / "bins" / f"{sample}.bin"
    png = Path(work_root) / "plots" / f"{sample}.bin.png"
    try:
        # Resume: reuse finished sample unless --force
        if bin_path.exists() and not force:
            bins = load_bin_file(bin_path)
            LOG.info("Skip finished sample %s (%d bins)", sample, len(bins))
            return SampleResult(sample=sample, bins=bins, n_reads=0, ok=True)

        style = peek_read_name_style(Path(r1))
        fq1 = sdir / f"{sample}.r1.fq.gz"
        fq2 = sdir / f"{sample}.r2.fq.gz"
        tagged1 = sdir / f"{sample}.r1.use.fq.gz"
        tagged2 = sdir / f"{sample}.r2.use.fq.gz"

        subsample_and_trim(Path(r1), Path(r2), fq1, fq2, n_pairs, read_len, seqtk)
        use1, use2, use_C = maybe_tag_headers(fq1, fq2, tagged1, tagged2, style)

        sam1 = sdir / f"{sample}.p1.sam"
        sam2 = sdir / f"{sample}.p2.sam"
        align_to_pseudo(use1, use2, Path(p1_fa), sam1, bwa, threads, use_C)
        align_to_pseudo(use1, use2, Path(p2_fa), sam2, bwa, threads, use_C)

        rlt = sdir / f"{sample}.filter.rlt"
        reads = discriminative_reads(sam1, sam2, rlt)
        bins = seq2bin(
            reads,
            chrom_len,
            win_frac=win_frac,
            min_bin=min_bin,
            work_dir=sdir,
            p1_name=p1_name,
            p2_name=p2_name,
        )
        bin_path.parent.mkdir(parents=True, exist_ok=True)
        write_sample_bin(bins, bin_path, sample, p1_name, p2_name)
        png.parent.mkdir(parents=True, exist_ok=True)
        plot_bins(bins, chrom_len, reads, png, p1_name, p2_name, sample)

        for p in (fq1, fq2, tagged1, tagged2, rlt):
            p.unlink(missing_ok=True)
        if not keep_align:
            for p in (sam1, sam2):
                p.unlink(missing_ok=True)

        return SampleResult(sample=sample, bins=bins, n_reads=len(reads), ok=True)
    except Exception as e:
        LOG.exception("Sample %s failed", sample)
        return SampleResult(sample=sample, ok=False, error=str(e))


# ---------------------------------------------------------------------------
# Stage 3: consensus bins, cM, bip, physical table
# ---------------------------------------------------------------------------

def consensus_physical_bins(
    sample_bins: Dict[str, List[BinRecord]],
    chrom_len: Dict[str, int],
    round_bp: int = 100_000,
) -> List[Tuple[str, float, int, int]]:
    """
    Population-cleared bins.
    Returns list of (chrom, mid_100kb_unit, start_bp, end_bp) sorted.
    """
    edges: Dict[str, set] = {c: {0} for c in chrom_len}
    for bins in sample_bins.values():
        for b in bins:
            if b.chrom not in edges:
                continue
            edges[b.chrom].add(int(b.end / round_bp + 0.5))
    markers: List[Tuple[str, float, int, int]] = []
    for chrom in sorted(edges, key=chrom_sort_key):
        pts = sorted(x for x in edges[chrom] if x > 0)
        if not pts:
            continue
        # marker at mid between consecutive edge units (MCD logic)
        prev = 0
        units = [0] + pts
        # also ensure chrom end
        end_u = int(chrom_len[chrom] / round_bp + 0.5)
        if units[-1] != end_u:
            units.append(end_u)
        mids = []
        for i in range(1, len(units)):
            mid = (units[i] + units[i - 1]) / 2.0
            mids.append(mid)
        # convert consecutive mids to physical intervals
        for i, mid in enumerate(mids):
            start = int(units[i] * round_bp)
            end = int(units[i + 1] * round_bp) if i + 1 < len(units) else chrom_len[chrom]
            if start < 1:
                start = 1
            if end > chrom_len[chrom]:
                end = chrom_len[chrom]
            if end <= start:
                continue
            markers.append((chrom, mid, start, end))
    return markers


def genotype_at_markers(
    bins: List[BinRecord], markers: List[Tuple[str, float, int, int]]
) -> List[str]:
    """Return A/H/B/U for each marker (A=parent1, B=parent2, H=hetero, U=missing)."""
    by_chr: Dict[str, List[BinRecord]] = {}
    for b in bins:
        by_chr.setdefault(b.chrom, []).append(b)
    gts = []
    for chrom, mid, start, end in markers:
        probe = (start + end) // 2
        call = "U"
        for b in by_chr.get(chrom, []):
            if b.start <= probe <= b.end:
                if b.origin == "parent1":
                    call = "A"
                elif b.origin == "parent2":
                    call = "B"
                else:
                    call = "H"
                break
        gts.append(call)
    return gts


def genetic_positions(
    markers: List[Tuple[str, float, int, int]],
    geno_matrix: List[List[str]],
) -> List[float]:
    """Per-chromosome cumulative cM from recombination counts (MCD-style)."""
    n_mark = len(markers)
    n_sample = len(geno_matrix)
    if n_mark == 0 or n_sample == 0:
        return []
    cm = [0.0] * n_mark
    i = 0
    while i < n_mark:
        chrom = markers[i][0]
        j = i
        while j < n_mark and markers[j][0] == chrom:
            j += 1
        cm[i] = 0.0
        for k in range(i + 1, j):
            count = 0
            for s in range(n_sample):
                a, b = geno_matrix[s][k - 1], geno_matrix[s][k]
                if a in ("U",) or b in ("U",):
                    continue
                if a != b:
                    if a == "H" or b == "H":
                        count += 1
                    else:
                        count += 2
            cm[k] = cm[k - 1] + (100.0 * count) / (2.0 * n_sample)
        i = j
    return [round(x, 4) for x in cm]


def load_phenotypes(pheno_list: Optional[Path], sample_order: List[str]) -> Tuple[List[str], List[List[str]]]:
    """pheno.list: one path per trait file; each file: FID IID value."""
    if pheno_list is None or not pheno_list.exists():
        return [], []
    trait_names: List[str] = []
    trait_vals: List[List[str]] = []
    with open(pheno_list) as fh:
        paths = [Path(x.strip()) for x in fh if x.strip() and not x.startswith("#")]
    for p in paths:
        name = p.name
        for suf in (".pheno", ".txt", ".tsv"):
            if name.endswith(suf):
                name = name[: -len(suf)]
                break
        # shorten kedou35_F2_xxx → xxx
        name = re.sub(r"^.*_F2_", "", name)
        mp: Dict[str, str] = {}
        with open(p) as pf:
            for line in pf:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                sid, val = parts[0], parts[2]
                if val.upper() in ("NA", "NAN", ".", "-", "*"):
                    val = "NA"
                mp[sid] = val
        trait_names.append(name)
        trait_vals.append([mp.get(s, "NA") for s in sample_order])
    return trait_names, trait_vals


def write_physical_table(
    path: Path,
    markers: List[Tuple[str, float, int, int]],
    cm: List[float],
) -> None:
    with open(path, "w") as out:
        out.write("bin_id\tchrom\tstart\tend\tmid_bp\tposition_cM\n")
        for i, ((chrom, _mid, start, end), g) in enumerate(zip(markers, cm), 1):
            mid_bp = (start + end) // 2
            out.write(f"bin_{i:05d}\t{chrom}\t{start}\t{end}\t{mid_bp}\t{g}\n")


def write_bip(
    path: Path,
    markers: List[Tuple[str, float, int, int]],
    cm: List[float],
    sample_order: List[str],
    geno_matrix: List[List[str]],
    trait_names: List[str],
    trait_vals: List[List[str]],
    pop_type: str,
    map_function: int = 1,
) -> None:
    """QTL IciMapping .bip (format ref: WheatDH.bip). P1=2, H=1, P2=0, miss=-1."""
    pop_code = POP_TYPE_CODE.get(pop_type.upper(), POP_TYPE_CODE.get(pop_type, 7))
    chroms = sorted({m[0] for m in markers}, key=chrom_sort_key)
    chr_counts = {c: 0 for c in chroms}
    for m in markers:
        chr_counts[m[0]] += 1

    # bip coding
    code = {"A": "2", "H": "1", "B": "0", "U": "-1"}
    if pop_type.upper() in ("DH", "F1DH", "RIL", "F1RIL"):
        # collapse hetero to missing for inbred/DH
        code["H"] = "-1"

    n_chr = len(chroms)
    n_ind = len(sample_order)
    n_trait = len(trait_names)

    lines: List[str] = []
    lines.append("!**********************Note: lines staring with \"!\" are remarks and will be ignored in the program********************")
    lines.append("!****************************************** General Information ******************************************************")
    lines.append("!Generated by binmap_pipeline.py")
    lines.append(" 1            !Indicator: 1 for mapping; 2 for simulation")
    lines.append(f" {pop_code}       \t!Mapping Population Type")
    lines.append(f" {map_function}       \t!Mapping Function (1 for Kosambi Function; 2 for Haldane Function; 3 for Morgan mapping function)")
    lines.append(" 2       \t!Marker Space Type (1 for intervals; 2 for positions)")
    lines.append(" 1       \t!Marker Space Unit(1 for centiMorgans; 2 for Morgan)")
    lines.append(f" {n_chr}      \t!Number of Chromosomes (or Linkage Groups)")
    lines.append(f" {n_ind}     \t!Size of the mapping population")
    lines.append(f" {n_trait}           !Number of traits")
    lines.append("")
    lines.append("!***************************************** Chromosomes and Linkage Maps **********************************************")
    lines.append("!Chromosome;     NumMarkers")
    for i, c in enumerate(chroms, 1):
        lines.append(f" Chr{i}        {chr_counts[c]}")
    lines.append("")
    lines.append("!MarkerName; Chromosome; position or interval as indicated by General Information")
    chr_index = {c: i for i, c in enumerate(chroms, 1)}
    for i, (m, g) in enumerate(zip(markers, cm), 1):
        chrom, _mid, start, end = m
        lines.append(f"bin_{i:05d}          {chr_index[chrom]}       {g:.4f}")

    lines.append("")
    lines.append("!*********************************** Marker Types ********************************************************************")
    lines.append("!Co-dominant marker: When heterozygote F1 is present, P1 is coded as 2, F1 is coded as 1, and P2 is coded as 0; Otherwise, P1 is coded as 2, and P2 is coded as 0. ")
    lines.append("!Missing marker is coded as -1. ")

    # genotypes: one marker block — marker name then n_ind codes (may wrap)
    for i, _m in enumerate(markers):
        codes = [code.get(geno_matrix[s][i], "-1") for s in range(n_ind)]
        # first line starts with marker name
        row = codes
        chunk = 24
        first = True
        pos = 0
        while pos < len(row):
            part = row[pos : pos + chunk]
            if first:
                lines.append(f"bin_{i+1:05d}                           " + " ".join(f"{x:>2}" for x in part))
                first = False
            else:
                lines.append(" " + " ".join(f"{x:>2}" for x in part))
            pos += chunk

    lines.append("")
    lines.append("!****************************************** Phenotypic Data **********************************************************")
    lines.append("!\"*\", \"NA\", \"na\", \".\",or \"-\" are used for missing phenotype")
    if n_trait == 0:
        lines.append("!No phenotypes provided")
    else:
        for tname, vals in zip(trait_names, trait_vals):
            # trait name + values, wrap ~14 per line like example
            flat = [tname] + vals
            # write with tabs on first segment
            chunk = 14
            # first line: name + some values
            first_vals = vals[:13]
            lines.append(tname + "\t" + "\t".join(first_vals))
            pos = 13
            while pos < len(vals):
                lines.append("\t" + "\t".join(vals[pos : pos + 14]))
                pos += 14

    path.write_text("\n".join(lines) + "\n")
    LOG.info("Wrote bip: %s", path)


def cleanup_work(outdir: Path, keep_pseudo: bool) -> None:
    """Remove heavy intermediates; keep SNP table, bins, plots, bip, physical, logs, pseudo."""
    work = outdir / "work"
    for name in [
        "parents.combined.g.vcf.gz",
        "parents.combined.g.vcf.gz.tbi",
        "parents.total.vcf.gz",
        "parents.total.vcf.gz.tbi",
        "parents.SNP.vcf.gz",
        "parents.SNP.vcf.gz.tbi",
        "parents.SNP.filter.vcf.gz",
        "parents.SNP.filter.vcf.gz.tbi",
        "parents.SNP.PASS.vcf.gz",
        "parents.SNP.PASS.vcf.gz.tbi",
    ]:
        p = work / name
        p.unlink(missing_ok=True)
        Path(str(p) + ".idx").unlink(missing_ok=True)

    for sub in ("parts", "vcf_parts"):
        d = work / sub
        if d.exists():
            shutil.rmtree(d, ignore_errors=True)

    sroot = outdir / "samples"
    if sroot.exists():
        for sd in sroot.iterdir():
            if not sd.is_dir():
                continue
            for f in sd.iterdir():
                if (
                    f.suffix in {".sam", ".fq"}
                    or f.name.endswith(".fq.gz")
                    or f.name.endswith(".rlt")
                    or f.name.endswith(".Parent1")
                    or f.name.endswith(".Parent2")
                    or ".sub.fq" in f.name
                ):
                    f.unlink(missing_ok=True)

    if not keep_pseudo:
        shutil.rmtree(outdir / "pseudo", ignore_errors=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Low-coverage parental-bin map pipeline → IciMapping .bip",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--p1-gvcf", type=Path, default=None, help="Parent1 gVCF (.g.vcf or .gz)")
    p.add_argument("--p2-gvcf", type=Path, default=None, help="Parent2 gVCF")
    p.add_argument(
        "--parent-vcf",
        type=Path,
        default=None,
        help="Optional joint/cohort VCF with both parents (skips CombineGVCFs/GenotypeGVCFs)",
    )
    p.add_argument("--fq-list", required=True, type=Path, help="List of R1 fastq paths (one per line)")
    p.add_argument("--pheno-list", type=Path, default=None, help="Optional list of phenotype file paths")
    p.add_argument(
        "--ref",
        type=Path,
        required=True,
        help="Reference fasta (must match gVCF/VCF contigs; nuclear Chr01–Chr20 expected)",
    )
    p.add_argument("--outdir", type=Path, default=Path("binmap_out"))
    p.add_argument("--p1-name", default=None, help="Parent1 VCF sample ID / display name")
    p.add_argument("--p2-name", default=None, help="Parent2 VCF sample ID / display name")
    p.add_argument("--pop-type", default="F2", help="IciMapping population type: F2|DH|RIL|...")
    p.add_argument("--read-len", type=int, default=100, help="Hard-trim reads to this length; 0 = no trim")
    p.add_argument("--depth", type=float, default=0.1, help="Target coverage per sample for subsample")
    p.add_argument("--parallel", type=int, default=5, help="Max samples in parallel")
    p.add_argument("--bwa-threads", type=int, default=4, help="bwa threads per sample")
    p.add_argument("--win-frac", type=float, default=0.8, help="Window dominance fraction for bin call")
    p.add_argument("--min-bin", type=int, default=300000, help="Minimum bin size (bp)")
    p.add_argument("--java-mem", default="32g")
    p.add_argument("--bwa", default=DEFAULT_BWA)
    p.add_argument("--seqtk", default=DEFAULT_SEQTK)
    p.add_argument("--gatk", default=DEFAULT_GATK)
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--keep-align", action="store_true", help="Keep per-sample SAM files")
    p.add_argument("--keep-pseudo", action="store_true", default=True, help="Keep pseudo genomes (default on)")
    p.add_argument("--no-keep-pseudo", action="store_true", help="Delete pseudo genomes at end")
    p.add_argument(
        "--skip-genotype",
        action="store_true",
        help="Deprecated/no-op: SNP table is auto-reused when present",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Recompute samples even if bins/<sample>.bin already exists",
    )
    p.add_argument("--max-samples", type=int, default=0, help="Debug: only first N samples (0=all)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    work = outdir / "work"
    work.mkdir(exist_ok=True)

    bwa = which_or(args.bwa, [DEFAULT_BWA, "bwa"])
    seqtk = which_or(args.seqtk, [DEFAULT_SEQTK, "seqtk"])
    gatk = which_or(args.gatk, [DEFAULT_GATK, "gatk"])
    bcftools = which_or(args.bcftools, ["bcftools"])
    samtools = which_or(args.samtools, ["samtools"])

    if args.parent_vcf is None and (args.p1_gvcf is None or args.p2_gvcf is None):
        raise SystemExit("Provide either --parent-vcf, or both --p1-gvcf and --p2-gvcf")

    if args.p1_name:
        p1_name = args.p1_name
    elif args.p1_gvcf is not None:
        p1_name = Path(args.p1_gvcf).name.split(".g.vcf")[0]
    else:
        raise SystemExit("--p1-name is required with --parent-vcf")

    if args.p2_name:
        p2_name = args.p2_name
    elif args.p2_gvcf is not None:
        p2_name = Path(args.p2_gvcf).name.split(".g.vcf")[0]
    else:
        raise SystemExit("--p2-name is required with --parent-vcf")

    chrom_len = read_fai(args.ref)
    genome_bp = sum(chrom_len.values())
    LOG.info("Nuclear genome size: %d bp (%d chroms)", genome_bp, len(chrom_len))

    # --- parents ---
    snp_table = outdir / "parents.diff.snp.tsv"
    if snp_table.exists():
        LOG.info("Reuse SNP table %s (delete it or pass fresh outdir to rebuild)", snp_table)
    else:
        vcf = stage_genotype(
            args.p1_gvcf,
            args.p2_gvcf,
            args.ref,
            work,
            gatk,
            args.java_mem,
            bcftools=bcftools,
            samtools=samtools,
            parent_vcf=args.parent_vcf,
            p1_name=p1_name,
            p2_name=p2_name,
        )
        extract_parent_snps(vcf, snp_table)

    pseudo_dir = outdir / "pseudo"
    pseudo_dir.mkdir(exist_ok=True)
    p1_fa = pseudo_dir / f"{p1_name}.pseudo.fa"
    p2_fa = pseudo_dir / f"{p2_name}.pseudo.fa"
    if not p1_fa.exists() or not p2_fa.exists():
        build_pseudo_genome(args.ref, snp_table, "p1", p1_fa)
        build_pseudo_genome(args.ref, snp_table, "p2", p2_fa)
    bwa_index(p1_fa, bwa)
    bwa_index(p2_fa, bwa)

    # --- sample list ---
    r1_list: List[Path] = []
    with open(args.fq_list) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            r1_list.append(Path(line))
    if args.max_samples and args.max_samples > 0:
        r1_list = r1_list[: args.max_samples]
    LOG.info("Samples: %d", len(r1_list))

    # estimate pairs for depth
    # use first sample read length (after trim if set)
    rl = args.read_len if args.read_len > 0 else peek_read_len(r1_list[0])
    target_bases = args.depth * genome_bp
    n_pairs = max(1, int(target_bases / (2 * rl)))
    LOG.info("Subsample ~%d pairs/sample (read_len=%d, depth=%.2f)", n_pairs, rl, args.depth)

    jobs = []
    for r1 in r1_list:
        sid = sample_id_from_r1(r1)
        r2 = fq_pair_from_r1(r1)
        if not r1.exists() or not r2.exists():
            LOG.error("Missing fq for %s: %s / %s", sid, r1, r2)
            continue
        jobs.append(
            (
                sid,
                str(r1),
                str(r2),
                str(p1_fa),
                str(p2_fa),
                str(outdir),
                chrom_len,
                n_pairs,
                args.read_len,
                bwa,
                seqtk,
                args.bwa_threads,
                p1_name,
                p2_name,
                args.win_frac,
                args.min_bin,
                args.keep_align,
                args.force,
            )
        )

    results: Dict[str, SampleResult] = {}
    # parallel: each sample runs bwa twice — limit concurrency
    with ProcessPoolExecutor(max_workers=max(1, args.parallel)) as ex:
        futs = {ex.submit(process_one_sample, j): j[0] for j in jobs}
        for fut in as_completed(futs):
            sid = futs[fut]
            res = fut.result()
            results[sid] = res
            if res.ok:
                LOG.info("OK %s  informative_reads=%d  bins=%d", sid, res.n_reads, len(res.bins))
            else:
                LOG.error("FAIL %s: %s", sid, res.error)

    ok_samples = [s for s, r in results.items() if r.ok and r.bins]
    ok_samples.sort()
    if not ok_samples:
        LOG.error("No successful samples; abort bip")
        return 1

    sample_bins = {s: results[s].bins for s in ok_samples}
    markers = consensus_physical_bins(sample_bins, chrom_len)
    geno_matrix = [genotype_at_markers(sample_bins[s], markers) for s in ok_samples]
    cm = genetic_positions(markers, geno_matrix)

    phys_path = outdir / "bin_physical.tsv"
    write_physical_table(phys_path, markers, cm)
    LOG.info("Physical bin table: %s  (%d bins)", phys_path, len(markers))

    trait_names, trait_vals = load_phenotypes(args.pheno_list, ok_samples)
    bip_path = outdir / ("binmap.bip" if trait_names else "binmap.nopheno.bip")
    write_bip(
        bip_path,
        markers,
        cm,
        ok_samples,
        geno_matrix,
        trait_names,
        trait_vals,
        args.pop_type,
    )

    # also always write nopheno if phenotypes present
    if trait_names:
        write_bip(
            outdir / "binmap.nopheno.bip",
            markers,
            cm,
            ok_samples,
            geno_matrix,
            [],
            [],
            args.pop_type,
        )

    # summary
    with open(outdir / "samples.summary.tsv", "w") as out:
        out.write("sample\tok\tn_informative\tn_bins\terror\n")
        for s in sorted(results):
            r = results[s]
            out.write(f"{s}\t{int(r.ok)}\t{r.n_reads}\t{len(r.bins)}\t{r.error}\n")

    keep_pseudo = not args.no_keep_pseudo
    cleanup_work(outdir, keep_pseudo=keep_pseudo)

    LOG.info("Done. Outputs in %s", outdir)
    LOG.info("  bins/*.bin  plots/*.bin.png  %s  %s", phys_path.name, bip_path.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
