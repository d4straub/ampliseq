process QIIME2_SEQFILTERTABLE {
    tag "${repseq}-filter-by-${table}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-ubuntu-2025.4-conda.yml"
    container "qiime2/amplicon:2025.4"

    input:
    path(table)
    path(repseq)

    output:
    path("filtered-sequences.qza"), emit: qza
    path "versions.yml"           , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime feature-table filter-seqs \\
        --i-data $repseq \\
        --i-table $table \\
        --o-filtered-data filtered-sequences.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
