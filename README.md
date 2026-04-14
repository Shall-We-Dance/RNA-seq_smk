# RNA-seq Snakemake Analysis Pipeline

This repository provides a reproducible Snakemake workflow for paired-end RNA-seq data. The pipeline follows a typical Smart-seq2 processing order: merge lane-level FASTQs to sample-level pairs, perform sample-level QC with fastp, align reads with STAR, quantify genes with featureCounts, and export both raw and normalized bigWig coverage tracks.

## Overview

### Inputs
- Paired-end raw FASTQ files (`*.fastq.gz`) for each sample.
- Multiple FASTQ pairs (lanes / technical replicates) may be assigned to the same sample.

### Outputs
- **QC**
  - `results/qc/fastp/<sample>/fastp.html/json` (per-sample)
  - `results/qc/multiqc/multiqc_report.html`
- **Alignment**
  - `results/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam` (optional, controlled by `output.keep_bam`)
  - `results/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam.bai` (optional, controlled by `output.keep_bam`)
- **Quantification**
  - `results/featurecount/totalRNA.counts.txt`
- **Signal tracks (bigWig)**
  - `results/bigwig/<sample>/<sample>.raw.bw`
  - `results/bigwig/<sample>/<sample>.normalized.bw`

## Pipeline steps

1. **Merge FASTQs (per-sample)**  
   Raw R1 files are concatenated into a single sample-level R1 file, and raw R2 files are concatenated into a single sample-level R2 file. (For gzip-compressed FASTQ files, stream concatenation with `cat` is valid.)

2. **fastp QC (per-sample)**  
   The merged FASTQs are processed by fastp once per sample. This produces sample-level HTML/JSON reports that are consistent with the reads used for alignment. The cleaned FASTQs from this step are temporary.

3. **STAR alignment**  
   Reads are aligned with STAR and a coordinate-sorted BAM is generated directly by STAR.

4. **featureCounts quantification**  
   Gene-level counts are generated from all sample BAM files.

5. **bigWig generation (raw + normalized)**  
   `bamCoverage` is run twice per sample:
   - raw signal (`--normalizeUsing None`)
   - normalized signal (`--normalizeUsing <bigwig.normalization>`, default `RPKM`)  
   Optional filtering removes chrM/scaffold/random/alt contigs from tracks.

## Requirements

- **Snakemake** (recommended: Snakemake >= 7)
- Conda / Mamba (recommended for environment management)
- STAR
- samtools
- fastp
- MultiQC
- deepTools (`bamCoverage`)
- subread (`featureCounts`)

## Installation

Recommended: create environments on-the-fly via Snakemake.

Example:
```bash
snakemake -s workflow/Snakefile --use-conda --cores 16
```

To speed up conda solves, consider using mamba:

```bash
snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --cores 16
```

## Configuration

Edit `config.yaml`.

Key fields:

* `reference.star_index`: prebuilt STAR genome index directory (must already contain index files such as `SA`)
* `reference.fasta`: reference FASTA
* `reference.gtf`: reference annotation GTF
* `reference.chrom_sizes`: chromosome sizes file (`.fai`-format 2-column file). If omitted, `<reference.fasta>.fai` is used.
* `samples`: mapping of sample name to lists of FASTQs for R1 and R2
* `star.quant_mode_gene_counts`: add STAR option `--quantMode GeneCounts` (default `true`)
* `output.keep_bam`: keep STAR BAM files (`false` by default to save disk)
* `bigwig.bin_size`: bigWig bin size
* `bigwig.normalization`: normalized bigWig method (`RPKM`, `CPM`, `BPM`, etc.)
* `bigwig.remove_chrM_and_scaffolds`: whether to mask chrM/scaffold-like contigs in bigWig

Example:

```yaml
reference:
  star_index: "/path/to/STAR/index"
  fasta: "/path/to/genome.fa"
  gtf: "/path/to/genes.gtf"
  chrom_sizes: "/path/to/genome.fa.fai"

samples:
  sampleA:
    R1:
      - "raw/sampleA_L001_R1.fastq.gz"
      - "raw/sampleA_L002_R1.fastq.gz"
    R2:
      - "raw/sampleA_L001_R2.fastq.gz"
      - "raw/sampleA_L002_R2.fastq.gz"

star:
  quant_mode_gene_counts: false

output:
  keep_bam: false
```

Notes:

* The workflow assumes paired-end reads and requires both R1 and R2 lists to be the same length per sample.
* STAR genome index construction is **not** performed in this workflow; build the STAR index in advance and point `reference.star_index` to that directory.
* For bigWig generation, ensure chromosome sizes are available (`reference.chrom_sizes` or `<reference.fasta>.fai`).

## Running the workflow

```bash
snakemake -s workflow/Snakefile --use-conda --cores 16
```

Dry run:

```bash
snakemake -s workflow/Snakefile -n
```

## Contact

For questions, please open an issue or contact the pipeline maintainer.
