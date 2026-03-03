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

// ── Import nf-core modules ────────────────────────────────────────────────────
// Each 'include' pulls in a self-contained process from the modules/ directory.
// After running `nf-core modules install <module>` these files will live under:
//   modules/nf-core/<tool>/<subtool>/main.nf

include { FASTQC } from './modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT } from './modules/nf-core/umitools/extract/main'
include { BBMAP_BBMERGE } from './modules/nf-core/bbmap/bbmerge/main'
include { BOWTIE_BUILD } from './modules/nf-core/bowtie/build/main'
include { BOWTIE_ALIGN } from './modules/nf-core/bowtie/align/main'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { UMITOOLS_DEDUP } from './modules/nf-core/umitools/dedup/main'
include { MAGECK_COUNT } from './modules/nf-core/mageck/count/main'
include { MAGECK_TEST } from './modules/nf-core/mageck/test/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'

// ── Helper: parse samplesheet ─────────────────────────────────────────────────
/*
    Expected samplesheet format (CSV):

    Single-end:
        sample,fastq_1,fastq_2,condition
        ctrl_rep1,/path/to/ctrl_rep1.fastq.gz,,control
        treat_rep1,/path/to/treat_rep1.fastq.gz,,treatment

    Paired-end:
        sample,fastq_1,fastq_2,condition
        ctrl_rep1,/path/to/ctrl_rep1_R1.fastq.gz,/path/to/ctrl_rep1_R2.fastq.gz,control
        treat_rep1,/path/to/treat_rep1_R1.fastq.gz,/path/to/treat_rep1_R2.fastq.gz,treatment

    - 'sample'    : unique sample identifier
    - 'fastq_1'   : path to R1 FASTQ (required)
    - 'fastq_2'   : path to R2 FASTQ (leave empty for single-end)
    - 'condition' : used later by MAGeCK to define control vs treatment groups

    Note: meta.single_end is set automatically based on whether fastq_2 is provided.
    For paired-end CRISPR screens, set --umi_discard_read 2 to discard R2 after
    UMI extraction so that only R1 (the sgRNA read) is passed to alignment.
*/
def parse_samplesheet(csv_file) {
    channel
        .fromPath(csv_file)
        .splitCsv(header: true)
        .map { row ->
            def single_end = !row.fastq_2 || row.fastq_2.trim() == ''
            def meta = [
                id         : row.sample,
                condition  : row.condition,
                single_end : single_end
            ]
            def reads = single_end
                ? [ file(row.fastq_1, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
            return [ meta, reads ]
        }
}

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    // ── Input validation ───────────────────────────────────────────────────────
    if (!params.input)               { error "Please provide a samplesheet with --input" }
    // if (!params.sgrna_library)       { error "Please provide an sgRNA library FASTA with --sgrna_library" }
    // if (!params.mageck_control)      { error "Please provide a control sgRNA list with --mageck_control" }
    // if (!params.mageck_treatment_id) { error "Please provide treatment sample ID(s) with --mageck_treatment_id" }
    // if (!params.mageck_control_id)   { error "Please provide control sample ID(s) with --mageck_control_id" }

    // -- 1. Load reads from samplesheet ----------------------------------------
    // 'reads_ch' holds tuples of:
    //   single-end: [ [id, condition, single_end:true],  [ /path/to/file.fastq.gz ] ]
    //   paired-end: [ [id, condition, single_end:false], [ /path/to/R1.fastq.gz, /path/to/R2.fastq.gz ] ]
    reads_ch = parse_samplesheet(params.input)

    // -- 2. FastQC on raw reads ------------------------------------------------
    FASTQC( reads_ch )
    fastqc_html = FASTQC.out.html
    fastqc_zip  = FASTQC.out.zip
    fastqc_versions = FASTQC.out.versions_fastqc

    // -- 2. Extract UMIs from raw reads ----------------------------------------
    // UMI-tools extracts UMIs based on the specified pattern and appends them to
    // the read headers.
    umi_reads = channel.empty()
    UMITOOLS_EXTRACT( reads_ch )
    umi_reads = UMITOOLS_EXTRACT.out.reads
    umi_log = UMITOOLS_EXTRACT.out.log

    // -- 3. BBMerge ------------------------------------------------
    // Merge paired end reads to single read.
    merged_reads = channel.empty()
    BBMAP_BBMERGE( umi_reads, params.interleave )
    merged_reads = BBMAP_BBMERGE.out.merged
    bbmerge_log = BBMAP_BBMERGE.out.log

    // -- 4. Trim adapter sequences -----------------------------------------------
    // Trim Galore removes adapter sequences and low-quality bases from the reads.
    // Trim to custom adapter sequences to remove non sgRNA sequence from amplicon.


    // -- 4. Build Bowtie index from sgRNA library FASTA -------------------------
    // Bowtie needs to pre-process the FASTA into an index before aligning.
    // We only need to build this ONCE regardless of how many samples we have,
    // so we pass the library as a plain file (not a channel of per-sample files).
    // BOWTIE_BUILD( 
    //     [ [id: 'sgrna_library'], file(params.sgrna_library) ]
    // )

    // -- 5. Align reads to sgRNA library ----------------------------------------
    // Bowtie v1 is used here because sgRNA sequences are short (~20 bp).
    // We combine each sample's trimmed reads with the single shared index.
    // .combine() pairs every item in TRIMGALORE.out.reads with the bowtie index.
    // BOWTIE_ALIGN(
    //     UMITOOLS_EXTRACT.out.reads.combine( BOWTIE_BUILD.out.index )
    // )

    // -- 6. Sort BAM ------------------------------------------------------------
    // Alignments come out of Bowtie in the order they were processed (not
    // coordinate order). SAMtools sort reorders them by genomic/library position,
    // which is required for indexing and for UMI deduplication.
    // SAMTOOLS_SORT( BOWTIE_ALIGN.out.bam )

    // -- 7. Index BAM -----------------------------------------------------------
    // Creates a .bai index file alongside the BAM. This lets tools rapidly
    // look up reads at any position without scanning the whole file.
    // SAMTOOLS_INDEX( SAMTOOLS_SORT.out.bam )

    // -- 8. Deduplicate by UMI --------------------------------------------------
    // UMI-tools groups reads that:
    //   (a) map to the same genomic position, AND
    //   (b) share the same UMI sequence
    // and collapses them into a single representative read.
    // This removes PCR duplicates that could inflate counts for some sgRNAs.
    //
    // We join the sorted BAM with its index file so the tool can access both.
    // bam_bai_ch = SAMTOOLS_SORT.out.bam
    //     .join( SAMTOOLS_INDEX.out.bai )

    // UMITOOLS_DEDUP( bam_bai_ch, params.dedup_stats )

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
    // ch_bams_collected = UMITOOLS_DEDUP.out.bam
    //     .collect { meta, bam -> bam }
    //     .map { bams ->
    //         def combined_meta = [
    //             id            : 'all_samples',
    //             sample_labels : UMITOOLS_DEDUP.out.bam.collect { meta, bam -> meta.id }.join(',')
    //         ]
    //         [ combined_meta, bams ]
    //     }

    // MAGECK_COUNT(
    //     ch_bams_collected,
    //     file(params.sgrna_library)
    // )

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
    // MAGECK_TEST(
    //     MAGECK_COUNT.out.count_table,
    //     params.mageck_treatment_id,   // e.g. "day14_rep1,day14_rep2"
    //     params.mageck_control_id,     // e.g. "day0_plasmid"
    //     file(params.mageck_control)   // non-targeting control sgRNA list (one per line)
    // )

    // FINAL -- MultiQC
    // MultiQC aggregates all the QC reports (FastQC, BBMerge stats, UMI-tools logs)
    // into a single interactive HTML report. This gives a nice overview of data
    // quality and processing metrics across all samples.
    MULTIQC(
        [ fastqc_html, bbmerge_log, umi_log ],
        params.outdir
    )

}

// ── Workflow completion summary ────────────────────────────────────────────
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
