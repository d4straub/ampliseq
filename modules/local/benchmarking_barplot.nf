process BENCHMARKING_BARPLOT {
    tag "${meta.classifier}_${meta.database}"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(observed)
    path(expected)

    output:
    path("*benchmarking_barplot_stats_long.tsv"), emit: long_tsv
    path("*benchmarking_barplot_stats_mean.tsv"), emit: mean_tsv
    path("*.png")                               , emit: pngs
    path("*.svg")                               , emit: svgs
    path("*benchmarking_barplot.log")           , emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.classifier}_${meta.database}_"
    def taxaseparator = task.ext.taxaseparator ?: ";"
    def mergedtaxasep = task.ext.mergedtaxasep ?: '\\|' //this is used for Sidle
    def fbeta = task.ext.fbeta ?: "2"
    def taxlevel = task.ext.taxlevel ?: "7"
    """
    benchmarking_barplot.r \\
        "$expected" \\
        "level-${taxlevel}.csv" \\
        "$taxaseparator" \\
        "$mergedtaxasep" \\
        "$fbeta" \\
        "$prefix" \\
        >${prefix}benchmarking_barplot.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
    END_VERSIONS
    """
}
