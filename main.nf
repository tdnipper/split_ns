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
include { CUTADAPT } from './modules/nf-core/cutadapt/main'
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'
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
    if (!params.sgrna_library)       { error "Please provide an sgRNA library FASTA with --sgrna_library" }
    // if (!params.mageck_control)      { error "Please provide a control sgRNA list with --mageck_control" }
    // if (!params.mageck_treatment_id) { error "Please provide treatment sample ID(s) with --mageck_treatment_id" }
    // if (!params.mageck_control_id)   { error "Please provide control sample ID(s) with --mageck_control_id" }

    ch_multiqc_files = channel.empty()

    // -- 1. Load reads from samplesheet ----------------------------------------
    // 'ch_reads' holds tuples of:
    //   single-end: [ [id, condition, single_end:true],  [ /path/to/file.fastq.gz ] ]
    //   paired-end: [ [id, condition, single_end:false], [ /path/to/R1.fastq.gz, /path/to/R2.fastq.gz ] ]
    ch_reads = parse_samplesheet(params.input)

    // -- 2. FastQC on raw reads ------------------------------------------------
    FASTQC( ch_reads )
    fastqc_html = FASTQC.out.html
    fastqc_zip  = FASTQC.out.zip
    ch_multiqc_files = ch_multiqc_files.mix( fastqc_html, fastqc_zip )

    // -- 2. Extract UMIs from raw reads ----------------------------------------
    // UMI-tools extracts UMIs based on the specified pattern and appends them to
    // the read headers.
    umi_reads = channel.empty()
    UMITOOLS_EXTRACT( ch_reads )
    umi_reads = UMITOOLS_EXTRACT.out.reads
    umi_log = UMITOOLS_EXTRACT.out.log
    ch_multiqc_files = ch_multiqc_files.mix( umi_log )

    // -- 3. BBMerge ------------------------------------------------
    // Merge paired end reads to single read.
    merged_reads = channel.empty()
    BBMAP_BBMERGE( umi_reads, params.interleave )
    merged_reads = BBMAP_BBMERGE.out.merged
        .map { meta, reads -> [ meta + [single_end: true], reads ] }
    ch_bbmerge_log = BBMAP_BBMERGE.out.log
    ch_multiqc_files = ch_multiqc_files.mix( ch_bbmerge_log )

    // -- 4. Trim adapter sequences -----------------------------------------------
    // Trim Galore removes adapter sequences and low-quality bases from the reads.
    // Trim to custom adapter sequences to remove non sgRNA sequence from amplicon.
    // trimmed_reads = channel.empty()
    // TRIMGALORE( merged_reads )
    // trimmed_reads = TRIMGALORE.out.reads
    // trim_log = TRIMGALORE.out.log
    // ch_multiqc_files = ch_multiqc_files.mix( trim_log )
    
    // -- 4. Trim adapter sequences -----------------------------------------------
    // Cutadapt removes adapter sequences and low-quality bases from the reads.
    // Trim to custom adapter sequences to remove non sgRNA sequence from amplicon.
    ch_trimmed_reads = channel.empty()
    CUTADAPT( merged_reads )
    ch_trimmed_reads = CUTADAPT.out.reads
    ch_cutadapt_log = CUTADAPT.out.log
    ch_multiqc_files = ch_multiqc_files.mix( ch_cutadapt_log )

    //-- 4. Build Bowtie index from sgRNA library FASTA -------------------------
    ch_library = Channel
        .fromPath(params.sgrna_library, checkIfExists: true)
        .map { lib -> [ [id: 'sgrna_library'], lib ] }

    BOWTIE_BUILD(ch_library)

    // -- 5. Align reads to sgRNA library ----------------------------------------
    // Bowtie v1 is used here because sgRNA sequences are short (~20 bp).
    // We combine each sample's trimmed reads with the single shared index.
    // .combine() pairs every item in TRIMGALORE.out.reads with the bowtie index.
    ch_aligned_reads = channel.empty()
    BOWTIE_ALIGN(
        ch_trimmed_reads, BOWTIE_BUILD.out.index, true
    )
    ch_aligned_reads = BOWTIE_ALIGN.out.bam

    // -- 6. Sort BAM ------------------------------------------------------------
    // Alignments come out of Bowtie in the order they were processed (not
    // coordinate order). SAMtools sort reorders them by genomic/library position,
    // which is required for indexing and for UMI deduplication.
    ch_sorted_indexed_bams = channel.empty()
    SAMTOOLS_SORT( ch_aligned_reads, ch_library, params.index_format )
    ch_sorted_indexed_bams = SAMTOOLS_SORT.out.bam
    
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
    ch_name_replacements = ch_reads
        .map { meta, reads ->
            def paired = !meta.single_end
            def suffixes = paired ? ['_1', '_2'] : ['']
            def mappings = []

            def fastq1_simplename = file(reads[0][0]).simpleName
            if (fastq1_simplename != meta.id) {
                mappings << [fastq1_simplename, "${meta.id}${suffixes[0]}"]
                if (paired) {
                    mappings << [file(reads[0][1]).simpleName, "${meta.id}${suffixes[1]}"]
                }
            }

            return mappings.collect { mapping -> mapping.join('\t') }
        }
        .flatten()
        .collectFile(name: 'name_replacement.txt', newLine: true)
        .ifEmpty([])

    ch_multiqc_config        = channel.fromPath("$projectDir/modules/nf-core/multiqc/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_logo          = params.multiqc_logo   ? channel.fromPath(params.multiqc_logo)   : channel.empty()
    ch_multiqc_sample_names = params.multiqc_sample_names ? channel.fromPath(params.multiqc_sample_names) : channel.value([])
    // MULTIQC(
    //     channel.of([[id: 'multiqc']])
    //         .combine( ch_multiqc_files.collect() )
    //         .combine( ch_multiqc_config.toList() )
    //         .combine( ch_multiqc_logo.toList() )
    //         .combine( ch_name_replacements )
    //         .combine( ch_multiqc_sample_names )
    // )
    // ch_multiqc_report = MULTIQC.out.report

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
