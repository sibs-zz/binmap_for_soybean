# binmap_for_soybean

低覆盖重测序亲本 binmap 流程 → QTL IciMapping `.bip`  
Low-coverage parental bin-map pipeline → QTL IciMapping `.bip`

适用于大豆等作物双亲杂交群体（F2 / RIL / DH 等）。核心脚本为单文件 Python 驱动，bin 划分调用同目录下的 `Seq2Bin_F24.pl`（与经典水稻流程逻辑一致）。

Designed for biparental mapping populations (F2 / RIL / DH, etc.) in soybean and other plants. One Python driver calls `Seq2Bin_F24.pl` for bin calling (same edge/hetero logic as the classic rice workflow).

仓库主页 / Homepage: https://github.com/sibs-zz/binmap_for_soybean

---

## 目录 | Contents

- [中文说明](#中文说明)
- [English](#english)

---

# 中文说明

## 依赖安装

### 软件

| 工具 | 用途 | 备注 |
|------|------|------|
| Python ≥ 3.8 | 主流程 | 标准库即可跑通；画图需 Pillow |
| Perl | 调用 `Seq2Bin_F24.pl` | 通常系统自带 |
| bwa | 比对到双亲伪基因组 | PATH 可找到即可 |
| seqtk | 抽样 / 截短 reads | |
| GATK | 亲本 gVCF → VCF（若不用 `--parent-vcf`） | |
| bcftools / samtools | VCF 处理、faidx | |
| Pillow (`pip install pillow`) | 绘制 `plots/*.bin.png` | 可选但推荐 |

```bash
pip install pillow
# 确保 bwa seqtk gatk bcftools samtools 在 PATH 中
```

本仓库文件：

```text
binmap_pipeline.py   # 主程序
Seq2Bin_F24.pl       # bin calling（须与主程序同目录）
examples/            # 输入格式示例、IciMapping bip 参考
```

## 输入文件

### 1) 亲本变异（二选一）

**A. 两个亲本 gVCF**

```bash
--p1-gvcf parent1.g.vcf.gz --p2-gvcf parent2.g.vcf.gz
```

**B. 已有联合 / 队列 VCF（含两个亲本样本）**

```bash
--parent-vcf cohort.vcf.gz --p1-name SAMPLE_P1 --p2-name SAMPLE_P2
```

参考基因组 `--ref` **必须**与 gVCF/VCF 的 contig 命名一致；脚本期望核染色体为 `Chr01`…`Chr20`（或 `ZH13.Chr01` 等形式）。

### 2) 子代 R1 列表 `--fq-list`

每行一个 **R1** fastq 路径；脚本自动配对 R2（`_1`↔`_2`、`_R1`↔`_R2` 等）：

```text
/path/to/DY-100_1.clean.fq.gz
/path/to/DY-101_1.clean.fq.gz
```

样本 ID 从文件名解析（如 `DY-100`）。

### 3) 表型列表 `--pheno-list`（可选）

每行一个表型文件路径。每个表型文件至少三列：`FID IID value`（空白分隔），用 **IID/FID 与样本 ID 匹配**：

```text
DY-100 DY-100 12.3
DY-101 DY-101 NA
```

见 `examples/pheno.list` 与 `examples/trait_example.pheno`。

## 快速运行

```bash
python3 binmap_pipeline.py \
  --parent-vcf cohort.pass.vcf.gz \
  --p1-name KD35-13 \
  --p2-name TZX-248-15 \
  --fq-list examples/fq.list \
  --pheno-list examples/pheno.list \
  --ref /path/to/ZH13.v2.fasta \
  --outdir binmap_out \
  --depth 0.1 \
  --parallel 5 \
  --bwa-threads 4
```

或从两个 gVCF 起步：

```bash
python3 binmap_pipeline.py \
  --p1-gvcf P1.g.vcf.gz --p2-gvcf P2.g.vcf.gz \
  --fq-list examples/fq.list \
  --ref /path/to/ref.fasta \
  --outdir binmap_out
```

断点续跑：已有 `bins/<sample>.bin` 会跳过该样本；加 `--force` 强制重算。

## 主要参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--ref` | （必填） | 参考基因组 fasta |
| `--depth` | `0.1` | 每样本目标覆盖度（抽样） |
| `--read-len` | `100` | 硬截短长度；`0` = 不截短 |
| `--parallel` | `5` | 并行样本数 |
| `--bwa-threads` | `4` | 每样本 bwa 线程 |
| `--pop-type` | `F2` | IciMapping 群体类型（写入 bip） |
| `--win-frac` | `0.8` | 窗口主基因型比例 |
| `--min-bin` | `300000` | 最小 bin 长度（bp），传给 Seq2Bin |
| `--force` | off | 重算已有 bin |
| `--max-samples` | `0` | 调试：只跑前 N 个（0=全部） |
| `-v` | off | DEBUG 日志 |

工具路径可用 `--bwa` / `--seqtk` / `--gatk` / `--bcftools` / `--samtools` 覆盖。

## 输出

```text
outdir/
  parents.diff.snp.tsv   # 双亲纯合差异 SNP
  pseudo/                # 双亲伪基因组 + bwa index
  bins/<sample>.bin      # 每样本 bin
  plots/<sample>.bin.png # bin 图
  bin_physical.tsv       # bin → 物理坐标 + cM
  binmap.bip             # 有表型时的 IciMapping 输入
  binmap.nopheno.bip     # 无表型 / 或额外一份无表型
  samples.summary.tsv
```

### 对接 QTL IciMapping

1. 打开 [QTL IciMapping](https://www.isbreeding.net/)（或同等 bip 兼容软件）。
2. 导入 `binmap.bip`（含表型）或 `binmap.nopheno.bip`。
3. 做连锁 / QTL 扫描，得到如 `qtlout.txt`。
4. 用 `LeftMarker` / `RightMarker`（`bin_xxxxx`）在 `bin_physical.tsv` 中查 `chrom / start / end`，对应到参考基因组物理位置。

bip 默认：群体类型 F2（码 7）、作图函数 **Kosambi**、染色体标签 `Chr1`…`Chr20`。基因型编码：P1=2，heterozygote=1，P2=0，缺失=-1。基因型和表型列按**同一套样本 ID 顺序**对齐。

---

# English

## Requirements

| Tool | Role |
|------|------|
| Python ≥ 3.8 | Driver |
| Perl | Runs `Seq2Bin_F24.pl` |
| bwa | Align to parental pseudo-genomes |
| seqtk | Subsample / trim reads |
| GATK | Parent gVCF → VCF (unless `--parent-vcf`) |
| bcftools / samtools | VCF / faidx |
| Pillow | Optional PNG plots (`pip install pillow`) |

Keep `Seq2Bin_F24.pl` next to `binmap_pipeline.py`.

## Inputs

1. **Parents** — either `--p1-gvcf` + `--p2-gvcf`, or `--parent-vcf` + `--p1-name` + `--p2-name`.
2. **`--fq-list`** — one R1 FASTQ path per line; R2 is inferred.
3. **`--pheno-list`** (optional) — paths to trait files (`FID IID value`).
4. **`--ref`** (required) — reference FASTA matching VCF contigs (`Chr01`–`Chr20`).

## Example

```bash
python3 binmap_pipeline.py \
  --parent-vcf cohort.pass.vcf.gz \
  --p1-name Parent1 --p2-name Parent2 \
  --fq-list examples/fq.list \
  --pheno-list examples/pheno.list \
  --ref /path/to/ref.fasta \
  --outdir binmap_out \
  --depth 0.1 --parallel 5
```

Resume skips finished `bins/<sample>.bin` unless `--force`.

## Key options

See the Chinese table above (`--depth`, `--read-len`, `--pop-type`, `--min-bin`, `--force`, …).

## Outputs

- `bins/`, `plots/` — per-sample bins and figures  
- `bin_physical.tsv` — physical ↔ cM map for markers  
- `binmap.bip` / `binmap.nopheno.bip` — **QTL IciMapping** input  

### Downstream IciMapping

Import `binmap.bip` into QTL IciMapping. Map QTL peaks to genome coordinates via `bin_physical.tsv` using `LeftMarker` / `RightMarker` from the QTL report.

## License

Code released for research use. Please cite your own study and acknowledge this pipeline when applicable.

## Authors

Maintained under [sibs-zz](https://github.com/sibs-zz/).
