# Histovar

Histovar is a somatic variant caller for cancer hotspots (e.g., K27M, G34R/V) in duplicated histone H3 genes. It uses statistical modeling and deep learning to call variants across paralogous regions.

## Key Algorithms

- **BAM/CRAM Scanning**: Uses indexed `pysam.AlignmentFile` queries for read-count checks. This skips empty coordinates to speed up execution.
- **Codon-Aware Discovery**: Resolves amino acid nomenclatures to 3-base codons. It scans all bases of the codon to detect mutations.
- **Expectation-Maximization (EM)**: Assigns multi-mapped reads to a locus using an EM algorithm.
  - Parses BWA `XA` tags for alternate mapping.
  - Uses ±500bp padding around transcript boundaries to link reads to the CDS.
  - Distributes reads with `MAPQ < 5` without `XA` tags across identical paralog clusters.
- **Deep Learning Filtering**: Differentiates variants from mapping artifacts using PyTorch (`DeepHistovarNet`) or Random Forest models based on MAPQ, strand bias, allele frequency, and EM posterior.
- **Histone Filtering**: Filters for genes following histone nomenclature (`H1-H4`, `HIST...`) during GTF parsing to limit the search space.

---

## Installation

### Docker (Recommended)
```bash
docker pull byoo/histovar:latest
```

### Apptainer / Singularity
Apptainer can run the Docker image directly.

```bash
# Pull the image
apptainer pull histovar.sif docker://byoo/histovar:latest

# Run the caller
apptainer run histovar.sif somatic \
    --tumor tumor.bam \
    --gtf MANE.GRCh38.v1.5.refseq_genomic.gtf.gz \
    --out results.tsv \
    --tsv \
    --genes all
```

### Local Environment
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install torch
pip install -e .
```

---

## Usage

### Basic Example (Docker)
```bash
docker run --rm -v $(pwd):/data byoo/histovar:latest somatic \
    --tumor /data/tumor.bam \
    --gtf /data/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz \
    --out /data/results.tsv \
    --tsv \
    --genes ""
```

### CLI Options
- `--tumor`: Path to tumor BAM/CRAM.
- `--normal`: (Optional) Path to matched normal BAM/CRAM for somatic subtraction.
- `--gtf`: Path to reference GTF.
- `--out`: Path to output VCF or TSV.
- `--tsv`: Output as TSV instead of VCF.
- `--genes`: Comma-separated list of genes. Default: `"H3C14,H3C15"`. Use `"all"` for all histones.
- `--fasta`: Path to reference FASTA. Required for CRAM files.
- `--min-reads`: Minimum alternate reads (Default: 2).
- `--min-vaf`: Minimum VAF (Default: 0.01).
- `--indels`: Enable indel discovery (experimental).
- `--model`: ML engine (`dl` or `rf`).
- `--build`: Genomic build (`hg38`, `hg19`).

---

## Output Columns (TSV)

* **CHROM**: Chromosome.
* **POS**: Genomic coordinate (1-based).
* **REF**: Reference base.
* **ALT**: Alternate base.
* **FILTER**: `PASS` if ML > 0.5, else `LowQual`.
* **AA**: Amino acid substitution (e.g., `K28M`).
* **GENES**: Gene name.
* **POST**: EM posterior confidence gap (0.0 to 1.0).
* **ML**: ML probability (0.0 to 1.0).
* **TUMOR_AF**: VAF in tumor.
* **NORMAL_AF**: VAF in normal.
* **SOMATIC_SCORE**: Calculated as `(TUMOR_AF - NORMAL_AF) * TUMOR_DEPTH`.
* **EVENT / MATE_POS**: Links ambiguous mutations across paralogs.

---

## Training Custom ML Models

1. **Simulate Data**: Generate synthetic BAMs with true variants and artifacts.
2. **Train Models**: Run training script with `MLScorer.train(X, y)`.
3. **Deploy**: Save `.pt` or `.pkl` files and reference them in the pipeline.

Benchmark performance:
```bash
histovar benchmark --truth data/simulated_truth.tsv --vcf calls.vcf
```
