process QIIME2_MERGE_TABLES {
    tag "${file1} + ${file2}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(file1, stageAs: 'one/*')
    path(file2, stageAs: 'two/*')

    output:
    path("table_merged.qza"), emit: qza
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: "error_on_overlapping_sample"
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime feature-table merge \\
        --i-tables $file1 \\
        --i-tables $file2 \\
        --o-merged-table table_merged.qza \\
        --p-overlap-method $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
