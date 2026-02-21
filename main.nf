#!/usr/bin/env nextflow

/*
========================================================================================
    CRISPR SCREEN PIPELINE — UMI-aware sgRNA quantification + MAGeCK ranking
========================================================================================
    Workflow:
        1. Extract UMIs from raw reads         (UMI-tools)
        2. Trim adapter sequences              (Trim Galore)
        3. Build Bowtie index from sgRNA lib   (Bowtie v1)
        4. Align trimmed reads                 (Bowtie v1)
        5. Sort & index BAM                    (SAMtools)
        6. Deduplicate reads by UMI            (UMI-tools dedup)
        7. Count sgRNAs & rank                 (MAGeCK count + test)
========================================================================================
*/

// ── Parameters ────────────────────────────────────────────────────────────────
// These are the inputs your pipeline needs. You can override any of these on the
// command line with --param_name value, e.g. --input samplesheet.csv

params.input          = null           // Path to samplesheet CSV (see format below)
params.sgrna_library  = null           // Path to sgRNA library FASTA
params.mageck_control = null           // File listing control sgRNA IDs (one per line)
params.mageck_treatment_id = null      // Sample label(s) for treatment in MAGeCK test e.g. "day14_rep1,day14_rep2"
params.mageck_control_id   = null      // Sample label(s) for control in MAGeCK test e.g. "day0_plasmid"
params.umi_pattern    = 'NNNNNNNNNN'     // UMI pattern: N = UMI base, X = non-UMI base
params.outdir         = 'results'
params.genome         = null           // Not used for sgRNA mapping, but kept for extensibility

// ── Input validation ──────────────────────────────────────────────────────────
// Nextflow will stop with a clear error if required params are missing

if (!params.input)               { error "Please provide a samplesheet with --input" }
if (!params.sgrna_library)       { error "Please provide an sgRNA library FASTA with --sgrna_library" }
if (!params.mageck_control)      { error "Please provide a control sgRNA list with --mageck_control" }
if (!params.mageck_treatment_id) { error "Please provide treatment sample ID(s) with --mageck_treatment_id" }
if (!params.mageck_control_id)   { error "Please provide control sample ID(s) with --mageck_control_id" }

// ── Import nf-core modules ────────────────────────────────────────────────────
// Each 'include' pulls in a self-contained process from the modules/ directory.
// After running `nf-core modules install <module>` these files will live under:
//   modules/nf-core/<tool>/<subtool>/main.nf

//include { UMITOOLS_EXTRACT } from './modules/nf-core/umitools/extract/main'
//include { TRIMGALORE } from './modules/nf-core/trimgalore/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'                                       
include { BOWTIE_BUILD } from './modules/nf-core/bowtie/build/main'
include { BOWTIE_ALIGN } from './modules/nf-core/bowtie/align/main'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { UMITOOLS_DEDUP } from './modules/nf-core/umitools/dedup/main'
include { MAGECK_COUNT } from './modules/nf-core/mageck/count/main'
include { MAGECK_TEST } from './modules/nf-core/mageck/test/main'

// ── Helper: parse samplesheet ─────────────────────────────────────────────────
/*
    Expected samplesheet format (CSV):
        sample,fastq_1,condition
        ctrl_rep1,/path/to/ctrl_rep1.fastq.gz,control
        treat_rep1,/path/to/treat_rep1.fastq.gz,treatment

    - 'sample'    : unique sample identifier
    - 'fastq_1'   : path to FASTQ (single-end; CRISPR screens are almost always SE)
    - 'condition' : used later by MAGeCK to define control vs treatment groups
*/
def parse_samplesheet(csv_file) {
    Channel
        .fromPath(csv_file)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id        : row.sample,
                condition : row.condition
            ]
            def fastq = file(row.fastq_1, checkIfExists: true)
            return [ meta, fastq ]
        }
}

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    // -- 1. Load reads from samplesheet ----------------------------------------
    // 'reads_ch' now holds tuples of: [ [id, condition], /path/to/file.fastq.gz ]
    reads_ch = parse_samplesheet(params.input)

    // -- 2. Extract UMIs --------------------------------------------------------
    // UMI-tools reads the UMI from the read header (default) or from a defined
    // pattern position. It moves the UMI sequence into the read name so that
    // downstream deduplication can group reads by their UMI later.
    //
    // The nf-core umitools/extract module expects: [ meta, [ fastq ] ]
    // We add an extra map() to wrap the single FASTQ in a list as required.
    // UMITOOLS_EXTRACT(
    //     reads_ch.map { meta, fq -> [ meta, [ fq ] ] },
    //     params.umi_pattern
    // )

    // -- 3. Trim adapters -------------------------------------------------------
    // Removes Illumina adapter sequences and low-quality bases from read ends.
    // For sgRNA screens this is especially important: reads must be trimmed
    // to the exact length of the sgRNA (typically 20 bp) for accurate alignment.
    // TRIMGALORE( UMITOOLS_EXTRACT.out.reads )

    // -- 2. Extract UMIs + Trim adapters --------------------------------------------------------
    // We can combine UMI extraction and adapter trimming into a single step using the nf-core
    // fastq_fastqc_umitools_trimgalore subworkflow. This runs UMI-tools extract and Trim Galore
    // together, and also produces FastQC reports for both raw and processed reads.
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
        reads_ch,
        params.umi_pattern
    )

    // -- 4. Build Bowtie index from sgRNA library FASTA -------------------------
    // Bowtie needs to pre-process the FASTA into an index before aligning.
    // We only need to build this ONCE regardless of how many samples we have,
    // so we pass the library as a plain file (not a channel of per-sample files).
    BOWTIE_BUILD(
        [ [id: 'sgrna_library'], file(params.sgrna_library) ]
    )

    // -- 5. Align reads to sgRNA library ----------------------------------------
    // Bowtie v1 is used here because sgRNA sequences are short (~20 bp).
    // We combine each sample's trimmed reads with the single shared index.
    // .combine() pairs every item in TRIMGALORE.out.reads with the bowtie index.
    BOWTIE_ALIGN(
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads.combine( BOWTIE_BUILD.out.index )
    )

    // -- 6. Sort BAM ------------------------------------------------------------
    // Alignments come out of Bowtie in the order they were processed (not
    // coordinate order). SAMtools sort reorders them by genomic/library position,
    // which is required for indexing and for UMI deduplication.
    SAMTOOLS_SORT( BOWTIE_ALIGN.out.bam )

    // -- 7. Index BAM -----------------------------------------------------------
    // Creates a .bai index file alongside the BAM. This lets tools rapidly
    // look up reads at any position without scanning the whole file.
    SAMTOOLS_INDEX( SAMTOOLS_SORT.out.bam )

    // -- 8. Deduplicate by UMI --------------------------------------------------
    // UMI-tools groups reads that:
    //   (a) map to the same genomic position, AND
    //   (b) share the same UMI sequence
    // and collapses them into a single representative read.
    // This removes PCR duplicates that could inflate counts for some sgRNAs.
    //
    // We join the sorted BAM with its index file so the tool can access both.
    bam_bai_ch = SAMTOOLS_SORT.out.bam
        .join( SAMTOOLS_INDEX.out.bai )

    UMITOOLS_DEDUP( bam_bai_ch )

    // -- 9. MAGeCK count --------------------------------------------------------
    // Because we did our OWN alignment (Bowtie) and deduplication (UMI-tools),
    // we pass the deduplicated BAMs directly to `mageck count --bam`.
    // This skips MAGeCK's internal alignment and uses our cleaner, UMI-deduped
    // read counts instead.
    //
    // The nf-core mageck/count module expects:
    //   input[0]: tuple val(meta), path(bam_files)   <- all BAMs bundled together
    //   input[1]: path(library)                       <- sgRNA library FASTA/CSV
    //
    // We collect all per-sample BAMs and bundle them under a single meta map,
    // because MAGeCK count runs once across ALL samples together.
    ch_bams_collected = UMITOOLS_DEDUP.out.bam
        .collect { meta, bam -> bam }
        .map { bams ->
            def combined_meta = [
                id            : 'all_samples',
                sample_labels : UMITOOLS_DEDUP.out.bam.collect { meta, bam -> meta.id }.join(',')
            ]
            [ combined_meta, bams ]
        }

    MAGECK_COUNT(
        ch_bams_collected,
        file(params.sgrna_library)
    )

    // -- 10. MAGeCK test (RRA ranking) -----------------------------------------
    // MAGeCK test uses the count table to rank sgRNAs and genes by depletion
    // or enrichment between conditions using the RRA (Robust Rank Aggregation)
    // algorithm. Genes whose sgRNAs all consistently drop out are ranked as
    // essential — that's the core readout of a KO screen.
    //
    // The nf-core mageck/test module expects:
    //   input[0]: tuple val(meta), path(count_table)
    //   input[1]: val(treatment)    <- sample label(s) for treatment condition
    //   input[2]: val(control)      <- sample label(s) for control condition
    //   input[3]: path(control_sgrna) <- optional: file of non-targeting controls
    //
    // IMPORTANT: Replace 'treatment_sample_name' and 'control_sample_name' below
    // with the actual sample IDs from your samplesheet (they must match the
    // column headers in the count table produced by MAGECK_COUNT).
    MAGECK_TEST(
        MAGECK_COUNT.out.count_table,
        params.mageck_treatment_id,   // e.g. "day14_rep1,day14_rep2"
        params.mageck_control_id,     // e.g. "day0_plasmid"
        file(params.mageck_control)   // non-targeting control sgRNA list (one per line)
    )
}

// ── Workflow completion summary ───────────────────────────────────────────────
workflow.onComplete {
    log.info """
    ========================================
    Pipeline complete!
    Status   : ${ workflow.success ? 'SUCCESS' : 'FAILED' }
    Results  : ${params.outdir}
    Duration : ${workflow.duration}
    ========================================
    """.stripIndent()
}
