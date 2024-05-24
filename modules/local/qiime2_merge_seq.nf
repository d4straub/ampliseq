process QIIME2_MERGE_SEQ {
    tag "${file1} + ${file2}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(file1, stageAs: 'one/*')
    path(file2, stageAs: 'two/*')

    output:
    path("seq_merged.qza"), emit: qza
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime feature-table merge-seqs \\
        --i-data $file1 \\
        --i-data $file2 \\
        --o-merged-data seq_merged.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
