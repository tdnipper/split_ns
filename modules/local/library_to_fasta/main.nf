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
    sed 's/\\r\$//' ${library_txt} \\
    | awk -F'\\t' '
        BEGIN { OFS="\\t" }
        {
            # If no tabs detected, try commas then whitespace
            if (NF < 2) {
                if (index(\$0, ",") > 0) {
                    n = split(\$0, a, ",")
                } else {
                    n = split(\$0, a, " +")
                }
                \$0 = a[1]
                for (i = 2; i <= n; i++) \$0 = \$0 OFS a[i]
            }
            # Trim leading/trailing whitespace from each field
            for (i = 1; i <= NF; i++) {
                gsub(/^[[:space:]]+|[[:space:]]+\$/, "", \$i)
            }
            # Skip header and blank lines, emit FASTA
            if (NR > 1 && \$1 != "") {
                print ">"\$1
                print \$2
            }
        }
    ' > ${prefix}.fa

    # Validate that output is non-empty
    if [ ! -s ${prefix}.fa ]; then
        echo "ERROR: No sequences extracted from library file '${library_txt}'." >&2
        echo "Check that it is tab- or comma-separated with columns: sgRNA, sequence, gene" >&2
        exit 1
    fi
    """
}
