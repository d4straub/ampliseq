process QIIME2_INSEQ {
    tag "${seq}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-ubuntu-2025.4-conda.yml"
    container "qiime2/amplicon:2025.4"

    input:
    path(seq)

    output:
    path("rep-seqs.qza"), emit: qza
    path "versions.yml", emit: versions

    script:
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime tools import \\
        --input-path "$seq" \\
        --type 'FeatureData[Sequence]' \\
        --output-path rep-seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
