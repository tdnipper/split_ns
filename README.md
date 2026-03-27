# split_ns

CRISPR screen pipeline for UMI-aware sgRNA quantification and MAGeCK ranking.

## Workflow

1. Extract UMIs from raw reads (UMI-tools)
2. Merge paired-end reads (BBMerge)
3. Trim flanking NS sequences (Cutadapt)
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

The `condition` column drives MAGeCK treatment/control grouping. Specify the control group with `--control_condition`. Every other unique condition value in the samplesheet is automatically treated as a separate treatment group and tested independently against the control via MAGeCK test.

### Example samplesheet

```csv
sample,fastq_1,fastq_2,condition
ctrl_rep1,data/ctrl_rep1_R1.fastq.gz,data/ctrl_rep1_R2.fastq.gz,control
ctrl_rep2,data/ctrl_rep2_R1.fastq.gz,data/ctrl_rep2_R2.fastq.gz,control
p1_rep1,data/p1_rep1_R1.fastq.gz,data/p1_rep1_R2.fastq.gz,p1
p2_rep1,data/p2_rep1_R1.fastq.gz,data/p2_rep1_R2.fastq.gz,p2
p3_rep1,data/p3_rep1_R1.fastq.gz,data/p3_rep1_R2.fastq.gz,p3
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
    --control_condition control
```

### Key parameters

| Parameter                | Default        | Description                                              |
|--------------------------|----------------|----------------------------------------------------------|
| `--input`                | (required)     | Path to samplesheet CSV                                  |
| `--mageck_library`       | (required)     | Path to sgRNA library tab-separated `.txt` file          |
| `--control_condition`    | (required)     | Condition value for the control group (all other conditions are tested against it) |
| `--mageck_control`       | `null`         | File listing control sgRNA IDs (one per line)            |
| `--outdir`               | `results`      | Output directory                                         |
| `--umi_pattern`          | `NNNNNNNNNN`   | UMI barcode pattern                                      |
| `--skip_fastqc`          | `false`        | Skip FastQC steps                                        |
| `--skip_trimming`        | `false`        | Skip adapter trimming                                    |

## UMI Deduplication

### How umi-tools dedup works

`umi_tools dedup` identifies PCR duplicates by grouping reads that share the same mapping position, then comparing their UMI sequences. Reads with identical or near-identical UMIs at the same position are collapsed to a single representative read, removing PCR-introduced copies while retaining true biological molecules. The tool processes a coordinate-sorted BAM file sequentially, maintaining an in-memory buffer of reads that are grouped and resolved according to the chosen method.

### Why sgRNA libraries are challenging for deduplication

In a standard sgRNA library, each guide sequence is its own short (~20 bp) reference contig, and every read that maps to a given guide aligns to position 0 of that contig — there is no coordinate diversity. This means all reads for any given guide form a single enormous bundle at one position. The more sophisticated deduplication methods (cluster, adjacency, directional) all require building a UMI similarity network across the entire bundle, which scales poorly when thousands or hundreds of thousands of reads pile up at a single position. In practice this can produce >60 GB memory usage and multi-hour runtimes on deeply sequenced screens.

### Available deduplication methods

umi-tools implements five methods for resolving which reads represent true molecules:

| Method | How it works | Accuracy | Speed |
|---|---|---|---|
| `unique` | Each distinct UMI sequence = one molecule. No network built, no error correction. | Lowest | Fastest |
| `percentile` | Drops UMIs whose count falls below 1% of the median at that position. No network built. | Low | Fast |
| `cluster` | Builds an undirected network connecting UMIs within 1 edit distance. Collapses each connected cluster to 1 molecule. | Moderate (undercounts) | Moderate |
| `adjacency` | Same network as cluster, but resolves clusters more carefully — allows multiple true molecules per cluster. | Good | Moderate |
| `directional` | Builds a directed network using both edit distance and read-count ratios to model PCR error biology. Most accurate. | Best | Slowest |

All three network-based methods (cluster, adjacency, directional) use a substring index for bundles with more than 25 unique UMIs to avoid a naive pairwise comparison of every UMI against every other, which substantially reduces computation for large bundles.

### Method used in this pipeline: `unique`

This pipeline runs `UMITOOLS_DEDUP` with `--method unique`. This means every distinct UMI sequence observed for a guide is counted as a separate original molecule, with no network construction and no UMI error correction.

**Why `unique` is the appropriate choice for this use case:**

The critical factor is the sgRNA reference structure. Because all reads for a given guide map to position 0 of a ~20 bp contig, every read forms a single massive bundle at one position. Network-based methods must construct a UMI similarity graph across the entire bundle — a step that becomes prohibitively expensive for deeply sequenced guides regardless of which network method is used. The `unique` method avoids this graph construction entirely, making it the most computationally tractable option for large screens.

The accuracy trade-off is acceptable for two reasons. First, any overcounting from sequencing errors in UMI bases is a systematic bias — it affects all guides and all samples proportionally. Since downstream analysis with MAGeCK performs between-sample normalisation and tests for relative enrichment or depletion between conditions, a uniform upward bias cancels out and does not affect guide rankings. Second, the 10 bp UMIs used in this pipeline limit the total number of artefactual UMI sequences that can arise from single-base sequencing errors, keeping the per-guide overcounting small relative to total read depth.

Note: the `--per-gene` and `--per-contig` flags are intentionally not used. These flags merge strand-specific read sub-bundles into a single larger bundle per guide, which increases the size of the UMI network and worsens the memory and runtime scaling problem. Without these flags, reads are grouped by strand within each contig, which keeps individual bundles smaller. For a well-designed sgRNA library where reads are expected to map predominantly to one strand, this has minimal impact on deduplication accuracy.

**When to consider a different method:**

If you need accurate absolute UMI counts — for example to estimate library complexity or PCR duplication rates — or if duplication rates vary substantially between samples, consider `--method directional`, the umi-tools developers' recommended default. It uses count-ratio logic to model PCR error biology and consistently outperforms all other methods in accuracy benchmarks. The trade-off is higher memory and runtime due to the large bundle sizes inherent in sgRNA data.

`cluster` is a pragmatic middle ground: it applies UMI error correction (unlike `unique`) while being cheaper than `directional`. Its undercounting bias — collapsing genuinely distinct UMIs that happen to be 1 edit distance apart into a single molecule — is, like the overcounting in `unique`, a systematic effect that largely cancels in differential comparisons.

Avoid `unique` if absolute per-guide molecule counts matter for your downstream analysis.

> Reference: Smith, Heger & Sudbery, *Genome Research* 2017.
