# Histovar

Histovar is a specialized, high-performance somatic variant caller explicitly designed to identify critical cancer hotspots (e.g., K27M, G34R/V) in the highly duplicated histone H3 gene family. Standard variant callers frequently fail or generate false positives in these regions due to mapping ambiguities. Histovar employs targeted statistical modeling and deep learning to confidently call variants across paralogous regions in diseases like Diffuse Midline Glioma (DMG).

## Key Algorithms

- **Ultra-Fast BAM/CRAM Scanning**: Uses indexed `pysam.AlignmentFile` queries to perform a "fast-path" read-count check. This skips 99% of empty genomic coordinates in milliseconds without extracting complex pileups, dramatically speeding up execution.
- **Codon-Aware Discovery**: Dynamically resolves modern (+1) and legacy amino acid nomenclatures to their 3-base genomic codons. It thoroughly scans all 3 bases of the codon to ensure it never misses second- or third-base mutations (e.g., K28M).
- **Paralog-Aware Expectation-Maximization (EM)**: Uses an unsupervised EM algorithm to probabilistically assign multi-mapped reads to their true originating locus.
  - Parses BWA `XA` tags to find alternate mapping locations.
  - Uses a robust ±500bp padding logic around transcript boundaries to ensure reads mapping to UTRs or introns correctly link to the CDS.
  - Automatically identifies reads with `MAPQ < 5` that lack `XA` tags (e.g., from DRAGEN or Bowtie2 aligners) and distributes them across known identical paralog clusters (like H3C14/H3C15).
- **Deep Learning Pileup Artifact Filtering**: Differentiates true variants from systematic mapping artifacts using PyTorch (`DeepHistovarNet`) or Random Forest models based on MAPQ, strand bias, allele frequency, and the EM posterior gap.
- **Automated Histone Filtering**: During GTF parsing, Histovar automatically filters for genes following standard histone nomenclature (`H1-H4`, `HIST...`). This prevents irrelevant genomic regions from being scanned and cluttering the results, while still allowing explicit gene whitelisting via the CLI.

---

## Installation

### Docker (Recommended)
Histovar is available as a pre-built cross-platform (amd64/arm64) container.
```bash
docker pull byoo/histovar:latest
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

Run the somatic caller. By default, it uses the Deep Learning (`dl`) engine and outputs a VCF file.

### Basic Example (Docker)
Ensure you mount your local data directory using the `-v` flag:
```bash
docker run --rm -v $(pwd):/data byoo/histovar:latest somatic \
    --tumor /data/tumor.bam \
    --gtf /data/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz \
    --out /data/results.tsv \
    --tsv \
    --genes ""
```

### CLI Options
- `--tumor`: Path to the tumor BAM/CRAM file.
- `--normal`: (Optional) Path to the matched normal BAM/CRAM file for somatic subtraction.
- `--gtf`: Path to the reference GTF (Supports Ensembl `chr1` and RefSeq `NC_000` chromosome formats).
- `--out`: Path to the output VCF or TSV file.
- `--tsv`: Flag to output results as a tab-separated TSV instead of VCF format.
- `--genes`: Comma-separated list of genes to scan. **Default is `"H3C14,H3C15"` to speed up testing.** To scan the entire histone family (including `H3-3A` and `H3-3B`), set this to `--genes "all"`.
- `--fasta`: Path to the reference FASTA file. **Required if using CRAM files.**
- `--min-reads`: Minimum number of alternate reads required to call a variant (Default: 2). Use higher values to reduce sequencing error sensitivity.
- `--min-vaf`: Minimum Variant Allele Frequency (VAF) required to call a variant (Default: 0.01). Useful to filter out ultra-low frequency sequencing noise.
- `--indels`: Flag to enable discovery and scoring of insertions and deletions (experimental).
- `--model`: ML scoring engine to use (`dl` for Deep Learning, `rf` for Random Forest).
- `--build`: Genomic build (`hg38`, `hg19`).

---

## Interpreting the Output

When running Histovar in `--tsv` mode, the output contains several highly informative columns:

* **CHROM**: The chromosome where the variant is located.
* **POS**: The exact genomic coordinate (1-based) of the variant.
* **REF**: The reference sequence base.
* **ALT**: The alternate (mutated) sequence base.
* **FILTER**: The filter status of the variant (`PASS` if ML > 0.5, else `LowQual`).
* **AA**: The amino acid substitution expected from this mutation (e.g., `K28M` or `G35R`).
* **GENES**: The gene where the mutation was called. If the `POST` score is highly ambiguous, Histovar merges the paralogs (e.g., `H3C14|H3C15`).
* **POST**: The Expectation-Maximization (EM) algorithm's posterior confidence gap (0.0 to 1.0). 
  * `POST = 1.0`: The reads map uniquely and unambiguously to the specified gene.
  * `POST < 0.2`: The reads perfectly multi-map to multiple paralogs, making the true origin highly ambiguous.
* **ML**: The Machine Learning probability (0.0 to 1.0) predicting whether this is a genuine biological mutation or a sequencing/mapping artifact. This is a complex score that evaluates multiple geometric and alignment properties (Depth, VAF, MAPQ, Strand Bias, Position Bias, and EM Posterior).
* **TUMOR_AF**: The Variant Allele Frequency in the tumor sample.
* **NORMAL_AF**: The Variant Allele Frequency in the matched normal sample. Reports `NA` if no normal BAM is provided.
* **SOMATIC_SCORE**: A raw, deterministic heuristic mathematically calculated as `(TUMOR_AF - NORMAL_AF) * TUMOR_DEPTH`. Unlike `ML` which provides a confidence probability, `SOMATIC_SCORE` provides the raw magnitude of the tumor-normal frequency delta, strictly scaled by read depth.
* **EVENT / MATE_POS**: If an ambiguous mutation is detected across two paralogs, Histovar links them together. Both rows will share a deterministic `EVENT` ID (e.g., `ambig_K28_H3C14_H3C15`), and their `MATE_POS` columns will accurately point to the exact genomic coordinate of the other paralog.

---

## Training Custom ML Models
Histovar relies on supervised Machine Learning to filter artifacts. The EM algorithm is unsupervised and runs on the fly, but the ML models must be trained on a simulated ground-truth dataset.

1. **Simulate Data**: Generate synthetic BAMs containing true variants and known artifacts.
2. **Train Models**: Execute a training script that passes pileup features into the `MLScorer.train(X, y)` function.
3. **Deploy**: Save the resulting `.pt` or `.pkl` files and reference them in the pipeline.

Evaluate your variant calling performance against a truth set at any time:
```bash
histovar benchmark --truth data/simulated_truth.tsv --vcf calls.vcf
```
