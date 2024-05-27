process BENCHMARKING_DIVERSITY {
    tag "${md5sum_version}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(table)
    path(tree)
    val(md5sum_version)


    output:
    path("*_distance.txt")  , emit: distance_txt
    path("*.svg")           , emit: svgs
    path("*.png")           , emit: pngs
    path("*.tsv")           , emit: tsv
    path("*.log")           , emit: log
    path("*.md5sum_version"), emit: md5sum_version
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def metric = task.ext.metric ?: "weighted_unifrac" // options are: 'generalized_unifrac', 'unweighted_unifrac', 'weighted_normalized_unifrac', 'weighted_unifrac'; described in https://forum.qiime2.org/t/alpha-and-beta-diversity-explanations-and-commands/2282
    def args = task.ext.args ?: "" // options are: https://docs.qiime2.org/2023.7/plugins/available/diversity/beta-phylogenetic/
    """
    # FIX: detecting a viable GPU on your system, but the GPU is unavailable for compute, causing UniFrac to fail.
    # COMMENT: might be fixed in version after QIIME2 2023.5
    export UNIFRAC_USE_GPU=N

    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime diversity beta-phylogenetic \\
        --i-phylogeny "${tree}" \\
        --i-table "${table}" \\
        --p-metric "${metric}" \\
        $args \\
        --p-threads ${task.cpus} \\
        --o-distance-matrix ${metric}_distance-matrix.qza

    qiime tools export \\
        --input-path ${metric}_distance-matrix.qza \\
        --output-path ${metric}_matrix
    mv ${metric}_matrix/distance-matrix.tsv ${metric}_distance.txt

    benchmarking_diversity.r ${metric}_distance.txt "${md5sum_version}" >${metric}_distance.log

    echo "md5sum_version $md5sum_version" > "${md5sum_version}.md5sum_version"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
        ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | sed 's/.*‘//' | sed 's/’.*//')
    END_VERSIONS
    """
}
