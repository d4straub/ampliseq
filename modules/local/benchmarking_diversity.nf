process BENCHMARKING_DIVERSITY {
    tag "${md5sum_version}"
    label 'process_single'

    conda "conda-forge::r-base=4.2.3 conda-forge::r-rmarkdown=2.22 conda-forge::r-tidyverse=2.0.0 conda-forge::r-knitr=1.43 conda-forge::r-dt=0.28 conda-forge::r-dtplyr=1.3.1 conda-forge::r-formattable=0.2.1 conda-forge::r-purrr=1.0.1 conda-forge::r-vegan=2.6_4 conda-forge::r-optparse=1.7.3 conda-forge::r-ggplot2=3.4.2 conda-forge::r-dplyr=1.1.2 conda-forge::r-data.table=1.14.8 conda-forge::r-patchwork=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    path(distances)
    val(md5sum_version)

    output:
    path("*.svg")           , emit: svgs
    path("*.png")           , emit: pngs
    path("*.tsv")           , emit: tsv
    path("*.log")           , emit: log
    path("*.md5sum_version"), emit: md5sum_version
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fbeta = task.ext.fbeta ?: 2
    def prefix = task.ext.prefix ?: ""
    def id_header = task.ext.id_header ?: "sequence"
    """
    benchmarking_diversity.r weighted_unifrac_distance.txt "$md5sum_version" >weighted_unifrac_distance.log
    benchmarking_diversity.r jaccard_distance.txt "$md5sum_version" >jaccard_distance.log
    benchmarking_diversity.r bray_curtis_distance.txt "$md5sum_version" >bray_curtis_distance.log
    benchmarking_diversity.r unweighted_unifrac_distance.txt "$md5sum_version" >unweighted_unifrac_distance.log

    echo "md5sum_version $md5sum_version" > "${md5sum_version}.md5sum_version"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
        ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | sed 's/.*‘//' | sed 's/’.*//')
    END_VERSIONS
    """
}
