process QIIME2_DIVERSITY_CORE {
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(metadata)
    path(table)
    path(tree)
    path(stats)
    val(mindepth)

    output:
    path("diversity_core/*_pcoa_results.qza")   , emit: pcoa
    path("diversity_core/*_vector.qza")         , emit: vector
    path("diversity_core/*_distance_matrix.qza"), emit: distance
    path("*_distance.txt")                      , emit: distance_txt
    path "versions.yml"                         , emit: versions
    path("*rarefaction.txt")                    , emit: depth

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def stats_cmd = stats ? "mindepth=\$(count_table_minmax_reads.py $stats minimum 2>&1)" : "mindepth=$mindepth"
    def metadata_cmd = task.ext.metadata_script ? "${task.ext.metadata_script} $metadata" : "cp $metadata ${metadata}.final.tsv"
    """
    # FIX: detecting a viable GPU on your system, but the GPU is unavailable for compute, causing UniFrac to fail.
    # COMMENT: might be fixed in version after QIIME2 2023.5
    export UNIFRAC_USE_GPU=N

    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # metadata manipulation
    $metadata_cmd

    $stats_cmd
    if [ \"\$mindepth\" -lt \"$mindepth\" ]; then mindepth=$mindepth; fi

    # report the rarefaction depth and return warning, if needed
    if [ \"\$mindepth\" -lt \"1000\" ]; then
        echo \$mindepth >\"WARNING The sampling depth of \$mindepth seems too small for rarefaction.txt\"
    elif [ \"\$mindepth\" -lt \"5000\" ]; then
        echo \$mindepth >\"WARNING The sampling depth of \$mindepth is very small for rarefaction.txt\"
    elif [ \"\$mindepth\" -lt \"10000\" ]; then
        echo \$mindepth >\"WARNING The sampling depth of \$mindepth is quite small for rarefaction.txt\"
    else
        echo \$mindepth >\"Use the sampling depth of \$mindepth for rarefaction.txt\"
    fi

    qiime diversity core-metrics-phylogenetic \\
        --m-metadata-file "${metadata}.final.tsv" \\
        --i-phylogeny ${tree} \\
        --i-table ${table} \\
        --p-sampling-depth \$mindepth \\
        --output-dir diversity_core \\
        --p-n-jobs-or-threads ${task.cpus} \\
        --verbose

    qiime tools export \\
        --input-path diversity_core/weighted_unifrac_distance_matrix.qza \\
        --output-path weighted_unifrac_distance_matrix
    mv weighted_unifrac_distance_matrix/distance-matrix.tsv weighted_unifrac_distance.txt
    qiime tools export \\
        --input-path diversity_core/jaccard_distance_matrix.qza \\
        --output-path jaccard_distance_matrix
    mv jaccard_distance_matrix/distance-matrix.tsv jaccard_distance.txt
    qiime tools export \\
        --input-path diversity_core/bray_curtis_distance_matrix.qza \\
        --output-path bray_curtis_distance_matrix
    mv bray_curtis_distance_matrix/distance-matrix.tsv bray_curtis_distance.txt
    qiime tools export \\
        --input-path diversity_core/unweighted_unifrac_distance_matrix.qza \\
        --output-path unweighted_unifrac_distance_matrix
    mv unweighted_unifrac_distance_matrix/distance-matrix.tsv unweighted_unifrac_distance.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
