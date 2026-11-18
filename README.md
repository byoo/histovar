# Histovar

Histovar is a somatic variant caller for cancer hotspots (e.g., K27M, G34R/V) in duplicated histone H3 genes. It uses Expectation-Maximization (EM) paralog resolution and statistical somatic modeling to call variants across highly homologous regions.

## How It Works

- **BAM/CRAM Scanning**: Uses indexed `pysam.AlignmentFile` queries for read-count checks. This skips empty coordinates to speed up execution.
- **Codon-Aware Discovery**: Resolves amino acid nomenclatures to 3-base codons. It scans all bases of the codon to detect mutations.
- **Expectation-Maximization (EM)**: Assigns multi-mapped reads (which are common in highly homologous histone gene families) to their most likely locus of origin.
  - Iteratively updates the prior probability (abundance) of each locus and the fractional assignment of reads based on mapping quality and match characteristics.
  - Parses BWA `XA` tags to identify alternate mapping candidates.
  - Uses configurable padding (default 500bp) around transcript boundaries to correctly map flanking reads to the coding sequences.
  - Distributes low-mapping-quality reads (`MAPQ < 5`) without `XA` tags across identical paralog clusters based on local posterior weights.
- **Somatic Scoring**: Evaluates the somatic probability of a variant using a Binomial Log-Likelihood Ratio (LLR) model.
  - Compares the somatic hypothesis ($H_1$, normal variant allele frequency (VAF) matches a sequencing error rate $e \approx 0.002$) against the germline hypothesis ($H_0$, normal VAF matches tumor VAF).
  - The LLR is Phred-scaled (`10 * log10(odds)`) to represent confidence. If a variant is highly present in the normal sample, the somatic score falls to 0.0.
- **Heuristic & Optional ML Filtering**: By default, scores candidate calls with a built-in heuristic based on variant allele frequency, read depth, mapping quality, and EM posterior gap (without penalizing low MAPQ on multi-mapped paralogs where `posterior_gap < 0.2`). Optionally supports loading custom pre-trained PyTorch (`DeepHistovarNet`) or scikit-learn models via `--model` and `--model-path`.
- **Histone Filtering**: Filters for genes following histone nomenclature (`H1-H4`, `HIST...`) during GTF parsing to limit the search space.
- **Reference & Hotspot Mapping**: Maps protein hotspots to genomic coordinates.
  - Dynamically parses GTF files to extract coding sequences (CDS) for histone genes.
  - Resolves legacy codon numbers to modern coordinates (which are +1 shifted due to the translation initiator Methionine).
  - Supports predefined cancer hotspots across H1, H2A, H2B, and H3 (e.g., H3 K27M/I, G34R/V/W, K36M, H2A E92K, H2B E76K/Q).

---

## Supported Genes and Hotspots

### Curated H3 Genes
Histovar explicitly maps and supports canonical H3 genes and their common aliases:
* **H3-3A** (`H3F3A`, `H3.3A`)
* **H3-3B** (`H3F3B`, `H3.3B`)
* **H3C1-H3C12** (e.g. `HIST1H3A` - `HIST1H3J`)
* **H3C13-H3C15** (`HIST2H3D`, `HIST2H3C`, `HIST2H3A`)
*(Other H1, H2A, H2B, and H4 histone genes matching the pattern `^(H[1-4][ABCDEF]?|HIST[1-5]H[1-4]|H[1-4]C[0-9])` are scanned dynamically when configured)*

### Predefined Hotspots
* **Histone H3**:
  - `R26` / `R27` (Arginine at position 27)
  - `K27` / `K28` (Lysine at position 28, variants: `M`, `I`)
  - `S28` / `S29` (Serine at position 29, variants: `C`)
  - `G34` / `G35` (Glycine at position 35, variants: `R`, `V`, `W`)
  - `K36` / `K37` (Lysine at position 37, variants: `M`)
* **Histone H2A**:
  - `E92` (Glutamate at position 92, variants: `K`)
  - `E121` (Glutamate at position 121, variants: `Q`)
* **Histone H2B**:
  - `G53` (Glycine at position 53, variants: `D`)
  - `E76` (Glutamate at position 76, variants: `K`, `Q`)
  - `E113` (Glutamate at position 113, variants: `K`)

---

## Installation

### Docker / Apptainer
```bash
# Docker
docker pull byoo/histovar:latest

# Apptainer
apptainer pull histovar.sif docker://byoo/histovar:latest
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

### Run with Docker
```bash
docker run --rm -v $(pwd):/data byoo/histovar:latest somatic \
    --tumor /data/tumor.bam \
    --gtf /data/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz \
    --out /data/results.tsv \
    --tsv \
    --genes "all"
```

### Run with Apptainer
```bash
apptainer run histovar.sif somatic \
    --tumor tumor.bam \
    --gtf MANE.GRCh38.v1.5.refseq_genomic.gtf.gz \
    --out results.tsv \
    --tsv \
    --genes "all"
```

### CLI Options
- `--tumor`: Path to the tumor BAM/CRAM alignment file. (Required)
- `--normal`: Path to the matched normal BAM/CRAM alignment file for somatic subtraction. (Optional)
- `--gtf`: Path to the reference gene annotation GTF file (can be `.gtf` or `.gtf.gz`). (Required)
- `--out`: Path to the output file (defaults to `calls.vcf`).
- `--tsv`: Flag to output in TSV format instead of VCF. (Default: False)
- `--genes`: Comma-separated list of genes to scan (e.g., `"H3C14,H3C15"`), or `"all"` to scan all histones. (Default: `"all"`)
- `--fasta`: Path to the reference genome FASTA file. (Required when using CRAM input files)
- `--min-reads`: Minimum number of alternate reads required supporting a variant to call it. (Default: 2)
- `--min-vaf`: Minimum variant allele frequency (VAF) required to call a variant. (Default: 0.01)
- `--pad`: Padding in base pairs around transcript boundaries to include mapping reads. (Default: 500)
- `--indels`: Flag to enable indel (insertions and deletions) discovery. (Default: False)
- `--model`: Scoring engine to use for filtering (`heuristic`, `dl`, `rf`, `boosting`, or `lr`). (Default: `heuristic`)
- `--model-path`: Path to a pre-trained model weights file (`.pt` or `.pkl`) when using `dl`, `rf`, `boosting`, or `lr`. (Optional)
- `--build`: Reference genome build assembly (`hg38` or `hg19`). (Default: `hg38`)

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
* **SOMATIC_SCORE**: Phred-scaled binomial log-likelihood ratio (LLR) comparing somatic ($H_1$: normal VAF $\approx 0.002$) and germline ($H_0$: normal VAF $\approx$ tumor VAF) hypotheses. If no normal sample is provided, falls back to `TUMOR_AF * TUMOR_DEPTH`.
* **EVENT / MATE_POS**: Links ambiguous mutations across paralogs.

---

## Training Custom ML Models

1. **Simulate Data**: Generate synthetic BAMs with true variants and artifacts.
2. **Train & Save Models**: Train a model with `MLScorer(model_type="dl").train(X, y)` and save weights via `scorer.save("model.pt")`.
3. **Deploy**: Pass `--model dl --model-path model.pt` (or `--model rf --model-path model.pkl`) to `histovar somatic`.

Benchmark performance:
```bash
histovar benchmark --truth data/simulated_truth.tsv --vcf calls.vcf
```
