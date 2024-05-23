process BENCHMARKING_SEQUENCES {
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    path(observed)
    path(expected)

    output:
    path("*.svg")                            , emit: svgs
    path("*.png")                            , emit: pngs
    path("*benchmarking_stats_mean.tsv")     , emit: mean_tsv
    path("*benchmarking_stats_long.tsv")     , emit: long_tsv
    path("*benchmarking_sequences_exact.log"), emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fbeta = task.ext.fbeta ?: 2
    def prefix = task.ext.prefix ?: ""
    def id_header = task.ext.id_header ?: "sequence"
    """
    benchmarking_sequences_exact.r \\
        "$expected" \\
        "$observed" \\
        "$fbeta" \\
        "$prefix" \\
        "$id_header" \\
        >${prefix}benchmarking_sequences_exact.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
    END_VERSIONS
    """
}
