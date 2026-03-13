process LIBRARY_TO_FASTA {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(library_txt)

    output:
    tuple val(meta), path("${prefix}.fa"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F'\\t' 'NR>1 { print ">"\$1; print \$2 }' ${library_txt} > ${prefix}.fa
    """
}
