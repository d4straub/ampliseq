process COMBINE_TABLE {
    tag "${meta.classifier}_${meta.database}"
    label 'process_low'

    conda "bioconda::bioconductor-biostrings=2.58.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-biostrings:2.58.0--r40h037d062_0' :
        'biocontainers/bioconductor-biostrings:2.58.0--r40h037d062_0' }"

    input:
    tuple val(meta),path(tax),path(table),path(seq)

    output:
    tuple val(meta), path("rel-table-ASV_with-*-tax.tsv")  , emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def outfile = "rel-table-ASV_with-${meta.classifier}-${meta.database}-tax.tsv"
    """
    combine_table.r ${table} ${seq} ${tax}
    mv combined_ASV_table.tsv ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
