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

include { FQ_LINT } from './modules/nf-core/fq/lint/main'
include { FASTQC as FASTQC_RAW     } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from './modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT } from './modules/nf-core/umitools/extract/main'
include { BBMAP_BBMERGE } from './modules/nf-core/bbmap/bbmerge/main'
include { CUTADAPT } from './modules/nf-core/cutadapt/main'
include { BOWTIE_BUILD } from './modules/nf-core/bowtie/build/main'
include { BOWTIE_ALIGN } from './modules/nf-core/bowtie/align/main'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { UMITOOLS_DEDUP } from './modules/nf-core/umitools/dedup/main'
include { MAGECK_COUNT } from './modules/nf-core/mageck/count/main'
include { MAGECK_TEST } from './modules/nf-core/mageck/test/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { LIBRARY_TO_FASTA } from './modules/local/library_to_fasta/main'

// ── Helper: parse samplesheet ─────────────────────────────────────────────────
/*
    Expected samplesheet format (CSV):

        sample,fastq_1,fastq_2,condition
        ctrl_rep1,/path/to/ctrl_rep1_R1.fastq.gz,/path/to/ctrl_rep1_R2.fastq.gz,control
        treat_rep1,/path/to/treat_rep1_R1.fastq.gz,/path/to/treat_rep1_R2.fastq.gz,treatment

    - 'sample'    : unique sample identifier
    - 'fastq_1'   : path to R1 FASTQ (required)
    - 'fastq_2'   : path to R2 FASTQ (required — this pipeline requires paired-end reads)
    - 'condition' : used later by MAGeCK to define control vs treatment groups
*/
def parse_samplesheet(csv_file) {
    channel
        .fromPath(csv_file)
        .splitCsv(header: true)
        .map { row ->
            if (!row.fastq_2 || row.fastq_2.trim() == '') {
                error "Sample '${row.sample}' is single-end. This pipeline requires paired-end reads (fastq_2 must be provided) because BBMerge is used to merge R1 and R2."
            }
            def meta = [
                id         : row.sample,
                condition  : row.condition,
                single_end : false
            ]
            def reads = [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
            return [ meta, reads ]
        }
}

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    // ── Input validation ───────────────────────────────────────────────────────
    if (!params.input)               { error "Please provide a samplesheet with --input" }
    if (!params.mageck_library)       { error "Please provide an sgRNA library .txt or .csv file with --mageck_library" }
    if (!params.control_condition)   { error "Please provide a control condition with --control_condition" }

    ch_multiqc_files = channel.empty()

    // -- 1. Load reads from samplesheet ----------------------------------------
    // 'ch_reads' holds tuples of:
    //   single-end: [ [id, condition, single_end:true],  [ /path/to/file.fastq.gz ] ]
    //   paired-end: [ [id, condition, single_end:false], [ /path/to/R1.fastq.gz, /path/to/R2.fastq.gz ] ]
    ch_reads = parse_samplesheet(params.input)

    // -- 2. Fastq Lint on paired-end reads ------------------------------------------------
    ch_reads_paired = ch_reads.filter { meta, reads -> !meta.single_end }
    FQ_LINT( ch_reads_paired )
    ch_multiqc_files = ch_multiqc_files.mix( FQ_LINT.out.lint.collect { _meta, lint -> lint } )
    ch_versions_fq = FQ_LINT.out.versions_fq

    // -- 2. FastQC on raw reads ------------------------------------------------
    FASTQC_RAW( ch_reads )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_RAW.out.html.collect { _meta, html -> html }, FASTQC_RAW.out.zip.collect { _meta, zip -> zip } )

    // -- Gate: block downstream until FQ_LINT + FASTQC_RAW complete -----------
    // If FQ_LINT fails the pipeline terminates before any further processing.
    // FQ_LINT only runs on paired-end reads; .ifEmpty(true) lets the gate pass
    // when there are no paired-end reads.
    ch_qc_done = FQ_LINT.out.lint
        .collect()
        .ifEmpty(true)
        .mix(
            FASTQC_RAW.out.zip.collect()
        )
        .collect()
        .map { true }

    ch_reads_gated = ch_reads.combine( ch_qc_done )
        .map { meta, reads, _gate -> [ meta, reads ] }

    // -- 3. Extract UMIs from raw reads ----------------------------------------
    // UMI-tools extracts UMIs based on the specified pattern and appends them to
    // the read headers.
    umi_reads = channel.empty()
    UMITOOLS_EXTRACT( ch_reads_gated )
    umi_reads = UMITOOLS_EXTRACT.out.reads
    ch_multiqc_files = ch_multiqc_files.mix( UMITOOLS_EXTRACT.out.log.collect { _meta, log -> log } )

    // -- 4. BBMerge ------------------------------------------------
    // Merge paired end reads to single read.
    merged_reads = channel.empty()
    BBMAP_BBMERGE( umi_reads, params.interleave )
    merged_reads = BBMAP_BBMERGE.out.merged
        .map { meta, reads -> [ meta + [single_end: true], reads ] }
    ch_multiqc_files = ch_multiqc_files.mix( BBMAP_BBMERGE.out.log.collect() )
    
    // -- 5. Trim adapter sequences -----------------------------------------------
    // Cutadapt removes adapter sequences and low-quality bases from the reads.
    // Trim to custom adapter sequences to remove non sgRNA sequence from amplicon.
    ch_trimmed_reads = channel.empty()
    CUTADAPT( merged_reads )
    ch_trimmed_reads = CUTADAPT.out.reads
    ch_multiqc_files = ch_multiqc_files.mix( CUTADAPT.out.log.collect { _meta, log -> log } )

    // -- 5b. FastQC on trimmed reads -------------------------------------------
    FASTQC_TRIMMED( ch_trimmed_reads )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_TRIMMED.out.html.collect { _meta, html -> html }, FASTQC_TRIMMED.out.zip.collect { _meta, zip -> zip } )

    //-- 6. Convert library TXT to FASTA, then build Bowtie index ---------------
    ch_mageck_library = Channel
        .fromPath(params.mageck_library, checkIfExists: true)
        .map { lib -> [ [id: 'sgrna_library'], lib ] }

    LIBRARY_TO_FASTA(ch_mageck_library)
    ch_library = LIBRARY_TO_FASTA.out.fasta

    BOWTIE_BUILD(ch_library)

    // -- 7. Align reads to sgRNA library ----------------------------------------
    // Bowtie v1 is used here because sgRNA sequences are short (~20 bp).
    // We combine each sample's trimmed reads with the single shared index.
    // .combine() pairs every item in CUTADAPT.out.reads with the bowtie index.
    ch_aligned_reads = channel.empty()
    BOWTIE_ALIGN(
        ch_trimmed_reads, BOWTIE_BUILD.out.index.first(), true
    )
    ch_aligned_reads = BOWTIE_ALIGN.out.bam

    // -- 8. Sort BAM ------------------------------------------------------------
    // Alignments come out of Bowtie in the order they were processed (not
    // coordinate order). SAMtools sort reorders them by genomic/library position,
    // which is required for indexing and for UMI deduplication.
    ch_sorted_indexed_bams = channel.empty()
    SAMTOOLS_SORT( ch_aligned_reads, ch_library.first(), params.index_format )
    ch_sorted_indexed_bams = SAMTOOLS_SORT.out.bam .join( SAMTOOLS_SORT.out.bai )
    
    // -- 9. Deduplicate by UMI --------------------------------------------------
    // UMI-tools groups reads that:
    //   (a) map to the same genomic position, AND
    //   (b) share the same UMI sequence
    // and collapses them into a single representative read.
    // This removes PCR duplicates that could inflate counts for some sgRNAs.
    ch_dedup_bams = channel.empty()
    UMITOOLS_DEDUP( ch_sorted_indexed_bams, params.dedup_stats )
    ch_dedup_bams = UMITOOLS_DEDUP.out.bam
    ch_multiqc_files = ch_multiqc_files.mix( UMITOOLS_DEDUP.out.log.collect { _meta, log -> log } )
    
    // -- 10. MAGeCK count --------------------------------------------------------
    // Because we did our OWN alignment (Bowtie) and deduplication (UMI-tools),
    // we pass the deduplicated BAMs directly to `mageck count --bam`.
    // This skips MAGeCK's internal alignment and uses our cleaner, UMI-deduped
    // read counts instead.
    //
    // We collect all per-sample BAMs and bundle them under a single meta map,
    // because MAGeCK count runs once across ALL samples together.
    ch_mageck_counts = channel.empty()
    ch_mageck_count_input = ch_dedup_bams
        .collect(flat: false)
        .map { tuples ->
            def ids  = tuples.collect { it[0].id }
            def bams = tuples.collect { it[1] }
            def meta = [id: 'all_samples', sample_labels: ids.join(',')]
            [meta, bams]
        }

    MAGECK_COUNT(
        ch_mageck_count_input,
        file(params.mageck_library)
    )
    ch_mageck_counts = MAGECK_COUNT.out.count

    // -- 11. MAGeCK test (RRA ranking) -----------------------------------------
    // MAGeCK test uses the count table to rank sgRNAs and genes by depletion
    // or enrichment between conditions using the RRA (Robust Rank Aggregation)
    // algorithm. Genes whose sgRNAs all consistently drop out are ranked as
    // essential — that's the core readout of a KO screen.
    //
    // Treatment/control sample IDs are derived automatically from the
    // 'condition' column in the samplesheet. Every unique condition that
    // is not --control_condition is treated as a separate treatment group
    // and tested independently against the control.
    //   --control_condition "day0"
    //   --mageck_control /path/to/control_sgrnas.txt  (optional)

    // Derive control sample IDs
    control_ids = ch_reads
        .map { meta, _reads -> meta }
        .filter { it.condition == params.control_condition }
        .ifEmpty { error "No samples matched control_condition '${params.control_condition}'. Check that this value matches entries in the 'condition' column of your samplesheet." }
        .map { it.id }
        .collect()
        .map { it.join(',') }

    // Derive treatment groups: every unique condition that isn't the control
    // produces one MAGECK_TEST run
    ch_treatment_groups = ch_reads
        .map { meta, _reads -> meta }
        .filter { it.condition != params.control_condition }
        .ifEmpty { error "No non-control conditions found in samplesheet. All samples have condition '${params.control_condition}'." }
        .map { [it.condition, it.id] }
        .groupTuple()
        .map { condition, ids -> [condition, ids.join(',')] }

    // Build MAGECK_TEST input: one run per treatment condition
    ch_mageck_test_input = ch_mageck_counts
        .combine(ch_treatment_groups)
        .combine(control_ids)
        .map { meta, counts, condition, t_ids, c_ids ->
            def new_meta = meta + [
                id: condition,
                treatment_ids: t_ids,
                control_ids: c_ids
            ]
            [new_meta, counts]
        }

    MAGECK_TEST(
        ch_mageck_test_input
    )

    // FINAL -- MultiQC
    // MultiQC aggregates all the QC reports (FastQC, BBMerge stats, UMI-tools logs)
    // into a single interactive HTML report. This gives a nice overview of data
    // quality and processing metrics across all samples.
    ch_name_replacements = ch_reads
        .map { meta, reads ->
            def paired = !meta.single_end
            def suffixes = paired ? ['_1', '_2'] : ['']
            def mappings = []

            def fastq1_simplename = reads[0].simpleName
            if (fastq1_simplename != meta.id) {
                mappings << [fastq1_simplename, "${meta.id}${suffixes[0]}"]
                if (paired) {
                    mappings << [reads[1].simpleName, "${meta.id}${suffixes[1]}"]
                }
            }

            return mappings.collect { mapping -> mapping.join('\t') }
        }
        .flatten()
        .collectFile(name: 'name_replacement.txt', newLine: true)
        .ifEmpty([])

    ch_multiqc_config        = channel.fromPath("$projectDir/modules/nf-core/multiqc/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config      ? channel.fromPath(params.multiqc_config)      : channel.empty()
    ch_multiqc_logo          = params.multiqc_logo        ? channel.fromPath(params.multiqc_logo)        : channel.empty()
    ch_multiqc_sample_names  = params.multiqc_sample_names ? channel.fromPath(params.multiqc_sample_names) : channel.value([])
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        ch_name_replacements,
        ch_multiqc_sample_names
    )
    ch_multiqc_report = MULTIQC.out.report

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
