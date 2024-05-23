process BENCHMARKING_TAXONOMY_PLOT {
    label 'process_single'

    conda "conda-forge::r-base=4.2.3 conda-forge::r-rmarkdown=2.22 conda-forge::r-tidyverse=2.0.0 conda-forge::r-knitr=1.43 conda-forge::r-dt=0.28 conda-forge::r-dtplyr=1.3.1 conda-forge::r-formattable=0.2.1 conda-forge::r-purrr=1.0.1 conda-forge::r-vegan=2.6_4 conda-forge::r-optparse=1.7.3 conda-forge::r-ggplot2=3.4.2 conda-forge::r-dplyr=1.1.2 conda-forge::r-data.table=1.14.8 conda-forge::r-patchwork=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    path(benchmarking_taxonomy_persequence_summary)
    path(benchmarking_taxonomy_persequence_stats_long)

    output:
    path("*_summary.tsv")   , emit: summary
    path("*_stats_long.tsv"), emit: stats_long
    path("*.log")           , emit: logs
    path("*.svg")           , emit: svgs
    path("*.png")           , emit: pngs
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "Fbeta,F1,recall,precision"
    """
    cat *benchmarking_taxonomy_persequence_summary.tsv > concatenated_benchmarking_taxonomy_persequence_summary.tsv
    benchmarking_taxonomy_plot_summary.r concatenated_benchmarking_taxonomy_persequence_summary.tsv >concatenated_benchmarking_taxonomy_persequence_summary.log

    cat *benchmarking_taxonomy_persequence_stats_long.tsv > concatenated_benchmarking_taxonomy_persequence_stats_long.tsv
    benchmarking_taxonomy_plot_stats.r concatenated_benchmarking_taxonomy_persequence_stats_long.tsv $args >concatenated_benchmarking_taxonomy_persequence_stats_long.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
        ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | sed 's/.*‘//' | sed 's/’.*//')
    END_VERSIONS
    """
}
