process BENCHMARKING_TAXONOMY {
    tag "${meta.classifier}_${meta.database}"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(observed), path(expected)

    output:
    path("*benchmarking_taxonomy_persequence.tsv")           , emit: long_tsv
    path("*benchmarking_taxonomy_persequence_summary.tsv")   , emit: summary_tsv
    path("*benchmarking_taxonomy_persequence_stats_long.tsv"), emit: stats_long
    path("*benchmarking_taxonomy_persequence_stats_mean.tsv"), emit: stats_mean
    path("*benchmarking_taxonomy_persequence.log")           , emit: log
    path("*.svg")                                            , emit: svgs
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.classifier}_${meta.database}"
    def args = task.ext.args ?: [
        '"sequence"',
        "${meta.database}" == "silva" ? '" "' : '""',
        "${meta.classifier}" == "qiime2" ? '"__"' : '""',
        '"2"'
        ].join(' ').trim()
        //'"sequence" "" "" 2' //$id_header $rm_after_string $rm_before_string $fbeta
    """
    benchmarking_taxonomy_persequence.r \\
        "$expected" \\
        "$observed" \\
        "$prefix" \\
        $args \\
        >${prefix}benchmarking_taxonomy_persequence.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
    END_VERSIONS
    """
}
