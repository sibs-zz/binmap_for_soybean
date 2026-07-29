# binmap_for_soybean

低覆盖重测序亲本 binmap 流程 → QTL IciMapping `.bip`  
Low-coverage parental bin-map pipeline → QTL IciMapping `.bip`

适用于大豆等作物双亲杂交群体（F2 / RIL / DH 等）。本流程改写自水稻全基因组重测序分型方法（Huang et al., *Genome Research*, 2009），用 Python 驱动 + `Seq2Bin_F24.pl` 完成 bin 划分，并输出可直接导入 [QTL IciMapping](https://www.isbreeding.net/) 的 `.bip`。

Designed for biparental mapping populations (F2 / RIL / DH, etc.) in soybean and other plants. This pipeline adapts the rice whole-genome resequencing genotyping approach (Huang et al., *Genome Research*, 2009): a Python driver plus `Seq2Bin_F24.pl` for bin calling, producing `.bip` files for [QTL IciMapping](https://www.isbreeding.net/).

仓库 / Repo: https://github.com/sibs-zz/binmap_for_soybean

---

## 目录 | Contents

- [中文说明](#中文说明)
- [English](#english)
- [Reference](#reference)

---

# 中文说明

## 1. Conda 环境安装（推荐）

需要已安装 [Miniconda](https://docs.conda.io/en/latest/miniconda.html) / [Mambaforge](https://github.com/conda-forge/miniforge)。

```bash
git clone https://github.com/sibs-zz/binmap_for_soybean.git
cd binmap_for_soybean

# 创建并激活环境（约含 Python、Pillow、bwa、seqtk、samtools、bcftools、GATK4、Perl、JDK）
conda env create -f environment.yml
conda activate binmap

# 检查
python -c "from PIL import Image; print('Pillow OK')"
bwa 2>&1 | head -1
seqtk 2>&1 | head -1
gatk --version | head -1
bcftools --version | head -1
samtools --version | head -1
perl -v | head -2
```

若 `gatk` 命令找不到，可试 `gatk4`，或在运行时指定：

```bash
python binmap_pipeline.py ... --gatk $(which gatk4 || which gatk)
```

### 仅用 pip（工具需自备）

若本机已有 `bwa` / `seqtk` / `gatk` / `bcftools` / `samtools` / `perl`：

```bash
pip install -r requirements.txt
```

### 仓库文件

```text
binmap_pipeline.py              # 主程序
Seq2Bin_F24.pl                  # bin calling（须与主程序同目录）
annotate_qtlout_physical.py     # 将 IciMapping qtlout 映射到物理坐标
environment.yml                 # conda 环境
requirements.txt                # Python 依赖（Pillow）
examples/                       # 输入格式示例、IciMapping bip 参考
example_run.sh                  # 示例启动脚本
```

## 2. 输入文件

### 亲本变异（二选一）

**A. 两个亲本 gVCF**

```bash
--p1-gvcf parent1.g.vcf.gz --p2-gvcf parent2.g.vcf.gz
```

**B. 已有联合 / 队列 VCF（含两个亲本）**

```bash
--parent-vcf cohort.vcf.gz --p1-name SAMPLE_P1 --p2-name SAMPLE_P2
```

`--ref` **必须**与 VCF contig 一致；核染色体期望为 `Chr01`…`Chr20`（或 `ZH13.Chr01` 等）。

### 子代 R1 列表 `--fq-list`

每行一个 R1；自动配对 R2（`_1`↔`_2`、`_R1`↔`_R2` 等）：

```text
/path/to/DY-100_1.clean.fq.gz
/path/to/DY-101_1.clean.fq.gz
```

### 表型列表 `--pheno-list`（可选）

每行一个表型文件。文件格式：`FID IID value`，ID 与样本名对应：

```text
DY-100 DY-100 12.3
DY-101 DY-101 NA
```

见 `examples/`。

## 3. 快速运行

```bash
conda activate binmap

python binmap_pipeline.py \
  --parent-vcf cohort.pass.vcf.gz \
  --p1-name Parent1 \
  --p2-name Parent2 \
  --fq-list examples/fq.list \
  --pheno-list examples/pheno.list \
  --ref /path/to/ref.fasta \
  --outdir binmap_out \
  --depth 0.1 \
  --parallel 5 \
  --bwa-threads 4
```

或从两个 gVCF：

```bash
python binmap_pipeline.py \
  --p1-gvcf P1.g.vcf.gz --p2-gvcf P2.g.vcf.gz \
  --fq-list examples/fq.list \
  --ref /path/to/ref.fasta \
  --outdir binmap_out
```

已有 `bins/<sample>.bin` 会跳过；`--force` 强制重算。

## 4. 主要参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--ref` | （必填） | 参考基因组 fasta |
| `--depth` | `0.1` | 每样本目标覆盖度 |
| `--read-len` | `100` | 硬截短；`0` = 不截短 |
| `--parallel` | `5` | 并行样本数 |
| `--bwa-threads` | `4` | 每样本 bwa 线程 |
| `--pop-type` | `F2` | IciMapping 群体类型 |
| `--win-frac` | `0.8` | 窗口主基因型比例 |
| `--min-bin` | `300000` | 最小 bin（bp） |
| `--force` | off | 重算已有 bin |
| `--max-samples` | `0` | 调试：前 N 个（0=全部） |
| `-v` | off | DEBUG 日志 |

可用 `--bwa` / `--seqtk` / `--gatk` / `--bcftools` / `--samtools` 指定路径。

## 5. 输出

```text
outdir/
  parents.diff.snp.tsv
  pseudo/
  bins/<sample>.bin
  plots/<sample>.bin.png
  bin_physical.tsv       # bin → 物理坐标 + cM
  binmap.bip             # IciMapping 输入（含表型）
  binmap.nopheno.bip
  samples.summary.tsv
```

### 对接 QTL IciMapping

1. 导入 `binmap.bip`（或 `binmap.nopheno.bip`）。
2. 连锁 / QTL 扫描 → 得到 `qtlout.txt`（含 `LeftMarker` / `RightMarker`，如 `bin_03082`）。
3. 用仓库脚本把 QTL 峰映射到参考基因组物理坐标：

```bash
python annotate_qtlout_physical.py \
  --qtlout qtlout.txt \
  --physical binmap_out/bin_physical.tsv \
  --out qtlout.with_physical.tsv
```

输出在原 QTL 表基础上增加：

| 列 | 含义 |
|----|------|
| PhysChrom | 物理染色体（如 Chr11） |
| PeakLeft_* / PeakRight_* | Left/RightMarker 对应的 start/end/mid |
| Peak_interval_start/end/Mb | 峰两侧 bin 合并区间 |
| CI_* / CI_Mb | 按 LeftCI–RightCI（cM）推算的物理置信区间 |

不写 `--out` 时，默认生成同目录下的 `qtlout.with_physical.tsv`。

bip 默认：F2、Kosambi、`Chr1`…`Chr20`；编码 P1=2 / H=1 / P2=0 / miss=-1。基因型与表型按同一样本 ID 顺序对齐。

---

# English

## 1. Conda install (recommended)

```bash
git clone https://github.com/sibs-zz/binmap_for_soybean.git
cd binmap_for_soybean
conda env create -f environment.yml
conda activate binmap
```

Includes Python, Pillow, Perl, bwa, seqtk, samtools, bcftools, GATK4, and JDK.  
Alternatively: `pip install -r requirements.txt` if bioinformatics tools are already on `PATH`.

Keep `Seq2Bin_F24.pl` next to `binmap_pipeline.py`. Helper script `annotate_qtlout_physical.py` maps IciMapping QTL peaks to physical coordinates.

## 2. Inputs

1. Parents: `--p1-gvcf` + `--p2-gvcf`, **or** `--parent-vcf` + `--p1-name` + `--p2-name`
2. `--fq-list` — R1 FASTQ paths (R2 inferred)
3. `--pheno-list` (optional) — trait files (`FID IID value`)
4. `--ref` (required) — FASTA matching VCF contigs (`Chr01`–`Chr20`)

## 3. Example

```bash
conda activate binmap
python binmap_pipeline.py \
  --parent-vcf cohort.pass.vcf.gz \
  --p1-name Parent1 --p2-name Parent2 \
  --fq-list examples/fq.list \
  --pheno-list examples/pheno.list \
  --ref /path/to/ref.fasta \
  --outdir binmap_out \
  --depth 0.1 --parallel 5
```

## 4. Outputs & IciMapping

- `bins/`, `plots/`, `bin_physical.tsv`
- `binmap.bip` → import into QTL IciMapping

After QTL scanning, annotate physical positions:

```bash
python annotate_qtlout_physical.py \
  --qtlout qtlout.txt \
  --physical binmap_out/bin_physical.tsv \
  --out qtlout.with_physical.tsv
```

Uses `LeftMarker` / `RightMarker` (`bin_xxxxx`) and CI (cM) against `bin_physical.tsv` to add chromosome bp intervals (`Peak_interval_*`, `CI_*`).

## License

Released for research use. Please cite the reference method below when applicable.

## Authors

Maintained under [sibs-zz](https://github.com/sibs-zz/).

---

# Reference

方法学参考 / Method reference（水稻全基因组重测序分型）：

> Huang X, Feng Q, Qian Q, Zhao Q, Wang L, Wang A, Guan J, Fan D, Weng Q, Huang T, Dong G, Sang T, Han B.  
> **High-throughput genotyping by whole-genome resequencing.**  
> *Genome Research*. 2009;19(6):1068–1076.  
> DOI: [10.1101/gr.089516.108](https://doi.org/10.1101/gr.089516.108) · PMID: [19420380](https://pubmed.ncbi.nlm.nih.gov/19420380/) · PMCID: [PMC2694477](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694477/)
