#!/usr/bin/env bash
# Example launcher — edit paths before running.
set -euo pipefail
cd "$(dirname "$0")"

python3 binmap_pipeline.py \
  --parent-vcf /path/to/parents.or.cohort.vcf.gz \
  --p1-name Parent1 \
  --p2-name Parent2 \
  --fq-list examples/fq.list \
  --pheno-list examples/pheno.list \
  --ref /path/to/reference.fasta \
  --outdir binmap_out \
  --depth 0.1 \
  --parallel 5 \
  --bwa-threads 4 \
  -v
