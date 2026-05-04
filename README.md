# Histovar

Histovar is a specialized, paralog-aware somatic variant caller designed specifically for the duplicated histone H3 gene family. It accurately resolves ambiguous multi-mapping reads across highly homologous regions to call critical cancer hotspots (e.g., K27M, G34R/V, K36M) in diseases like Diffuse Midline Glioma.

## Key Features

- **Deep Learning Pileup Artifact Filtering**: Uses a PyTorch Multi-Layer Perceptron (`DeepHistovarNet`) to differentiate true variants from multi-mapping artifacts based on pileup signature matrices.
- **Paralog-Aware Expectation-Maximization**: Utilizes an unsupervised EM algorithm to probabilistically assign multi-mapping reads to their true originating locus (e.g., distinguishing H3C14 from H3C15).
- **BAM-native Parsing**: Scans for `XA` and `SA` tags to extract candidate multi-mapped loci directly from the aligner.
- **Classic Machine Learning Rescoring**: Embeds Random Forest/Boosting models as a fallback to filter systematic mapping artifacts based on features like strand bias, MAPQ, and posterior gap confidence.
- **Annotation-driven Hotspots**: Dynamically resolves coordinates for complex aliases (e.g., `H3F3A` -> `H3-3A`) and codon mappings across hg19, hg38, and T2T.

## Install

### Local Environment
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install torch  # Required for deep learning functionality
pip install -e .
```

### Docker
Alternatively, you can run Histovar entirely within a container to avoid environment dependencies:
```bash
docker build -t histovar:v2 .
```

## Training the ML Models

Histovar relies on supervised Machine Learning (PyTorch or Random Forest) to filter artifacts. Before using it on real samples, you must train the models using a simulated ground-truth dataset. 

*(Note: The EM algorithm for paralog assignment is an unsupervised statistical method and **does not** require training. It runs automatically on the fly before the ML step.)*

1. **Simulate Data**: Generate synthetic BAMs containing true variants and known artifacts.
2. **Train Models**: Execute a training script that passes pileup features into the `MLScorer.train(X, y)` function.
3. **Save/Load Weights**: Save the resulting `.pt` (PyTorch) or `.pkl` (Scikit-Learn) files to deploy the trained model in production.

## Run

### Local Execution

Run the somatic caller on a matched tumor-normal pair using the Deep Learning engine (default):
```bash
histovar somatic --tumor tumor.bam --normal normal.bam --model dl --out calls.vcf
```

#### New Options:
- `--tsv`: Output results in a tab-separated format (useful for downstream analysis in Excel/Pandas).
- `--indels`: Enable discovery and scoring of insertions and deletions (experimental).

Example with TSV and Indels:
```bash
histovar somatic --tumor tumor.bam --gtf ref.gtf --out results.tsv --indels
```

To run in **tumor-only mode**:
```bash
histovar somatic --tumor tumor.bam --model dl --out tumor_only.vcf
```

Alternatively, run using the classic Random Forest engine:
```bash
histovar somatic --tumor tumor.bam --normal normal.bam --model rf --out calls.vcf
```

Evaluate variant calling performance against a simulated truth set:
```bash
histovar benchmark
```

### Docker Execution
To run Histovar using the container, you must mount your local data directory into the container using the `-v` flag:

```bash
docker run --rm -v $(pwd):/data histovar:v2 somatic \
    --tumor /data/tumor.bam \
    --normal /data/normal.bam \
    --out /data/calls.vcf \
    --model dl
```
