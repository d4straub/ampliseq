process QIIME2_INTREE {
    tag "${meta.id}:${meta.model}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-ubuntu-2025.4-conda.yml"
    container "qiime2/amplicon:2025.4"

    input:
    tuple val(meta), path(tree)

    output:
    path("tree.qza")   , emit: qza
    path "versions.yml", emit: versions

    script:
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime tools import \\
        --type 'Phylogeny[Rooted]' \\
        --input-path $tree \\
        --output-path tree.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
