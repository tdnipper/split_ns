# split_ns

CRISPR screen pipeline for UMI-aware sgRNA quantification and MAGeCK ranking.

## Workflow

1. Extract UMIs from raw reads (UMI-tools)
2. Merge paired-end reads (BBMerge)
3. Trim adapter sequences (Cutadapt)
4. Build Bowtie index from sgRNA library (Bowtie v1)
5. Align trimmed reads (Bowtie v1)
6. Sort & index BAM (SAMtools)
7. Deduplicate reads by UMI (UMI-tools dedup)
8. Count sgRNAs (MAGeCK count)
9. Rank genes by enrichment/depletion (MAGeCK test)
10. Aggregate QC reports (MultiQC)

## Samplesheet format

The input samplesheet is a CSV with the following columns:

| Column      | Required | Description                                                       |
|-------------|----------|-------------------------------------------------------------------|
| `sample`    | Yes      | Unique sample identifier (becomes the column header in the count table) |
| `fastq_1`   | Yes      | Path to R1 FASTQ file                                            |
| `fastq_2`   | No       | Path to R2 FASTQ file (leave empty for single-end)               |
| `condition` | Yes      | Experimental group label (e.g., `treatment`, `control`, `day14`) |

The `condition` column drives MAGeCK treatment/control grouping. All samples sharing the same condition value are grouped together. Use the `--treatment_condition` and `--control_condition` parameters to specify which condition values represent each group.

### Example samplesheet

```csv
sample,fastq_1,fastq_2,condition
ctrl_rep1,data/ctrl_rep1_R1.fastq.gz,data/ctrl_rep1_R2.fastq.gz,control
ctrl_rep2,data/ctrl_rep2_R1.fastq.gz,data/ctrl_rep2_R2.fastq.gz,control
treat_rep1,data/treat_rep1_R1.fastq.gz,data/treat_rep1_R2.fastq.gz,treatment
treat_rep2,data/treat_rep2_R1.fastq.gz,data/treat_rep2_R2.fastq.gz,treatment
```

## sgRNA library file

You only need to provide a single sgRNA library file via `--mageck_library`. The pipeline automatically converts it to FASTA format for Bowtie index building — no separate FASTA file is needed.

The library file must be a **tab-separated text file** with a header row and the following columns:

| Column   | Description                                      |
|----------|--------------------------------------------------|
| `sgRNA`  | Unique sgRNA identifier (e.g., `gene1_sg1`)      |
| `sequence` | sgRNA protospacer sequence (e.g., `ACGTACGTACGTACGTACGT`) |
| `gene`   | Target gene name (e.g., `gene1`)                 |

### Example library file

```
sgRNA	sequence	gene
TP53_sg1	ACGTACGTACGTACGTACGT	TP53
TP53_sg2	TGCATGCATGCATGCATGCA	TP53
BRCA1_sg1	AAGGCCTTAAGGCCTTAAGG	BRCA1
BRCA1_sg2	CCTTAAGGCCTTAAGGCCTT	BRCA1
NT_sg1	GACGATAGCGAGCTAGCTAG	non-targeting
```

This is the standard MAGeCK library format. The pipeline uses columns 1 and 2 (sgRNA name and sequence) to generate the FASTA reference for alignment, and passes the full file to MAGeCK count for sgRNA quantification.

## Usage

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --mageck_library /path/to/library.txt \
    --treatment_condition treatment \
    --control_condition control
```

### Key parameters

| Parameter                | Default        | Description                                              |
|--------------------------|----------------|----------------------------------------------------------|
| `--input`                | (required)     | Path to samplesheet CSV                                  |
| `--mageck_library`       | (required)     | Path to sgRNA library tab-separated `.txt` file          |
| `--treatment_condition`  | (required)     | Condition value for the treatment group                  |
| `--control_condition`    | (required)     | Condition value for the control group                    |
| `--mageck_control`       | `null`         | File listing control sgRNA IDs (one per line)            |
| `--outdir`               | `results`      | Output directory                                         |
| `--umi_pattern`          | `NNNNNNNNNN`   | UMI barcode pattern                                      |
| `--skip_fastqc`          | `false`        | Skip FastQC steps                                        |
| `--skip_trimming`        | `false`        | Skip adapter trimming                                    |
