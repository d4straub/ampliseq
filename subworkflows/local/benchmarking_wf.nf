include { BENCHMARKING_SEQUENCES     } from '../../modules/local/benchmarking_sequences'
include { BENCHMARKING_TAXONOMY      } from '../../modules/local/benchmarking_taxonomy'
include { BENCHMARKING_TAXONOMY_PLOT } from '../../modules/local/benchmarking_taxonomy_plot'
include { BENCHMARKING_BARPLOT       } from '../../modules/local/benchmarking_barplot'
include { QIIME2_INASV as QIIME2_INASV_BM } from '../../modules/local/qiime2_inasv'
include { QIIME2_MERGE_TABLES        } from '../../modules/local/qiime2_merge_tables'
include { QIIME2_MERGE_SEQ           } from '../../modules/local/qiime2_merge_seq'
include { QIIME2_INSEQ as QIIME2_INSEQ_BM } from '../../modules/local/qiime2_inseq'
include { QIIME2_TREE as QIIME2_TREE_BM   } from '../../modules/local/qiime2_tree'
include { BENCHMARKING_DIVERSITY     } from '../../modules/local/benchmarking_diversity'

workflow BENCHMARKING_WF {
    take:
    val_md5sum_version
    ch_asv
    ch_expected_sequences
    ch_tax
    ch_expected_taxonomies
    ch_barplot
    ch_expected_barplot
    params_benchmarking_diversity
    ch_observed_table_qza
    ch_observed_seq_qza

    main:
    ch_benchmarking_versions = Channel.empty()

    BENCHMARKING_SEQUENCES ( ch_asv, ch_expected_sequences, val_md5sum_version )
    ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_SEQUENCES.out.versions)

    BENCHMARKING_TAXONOMY ( ch_tax.combine(ch_expected_taxonomies) )
    ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_TAXONOMY.out.versions)

    BENCHMARKING_TAXONOMY_PLOT ( BENCHMARKING_TAXONOMY.out.summary_tsv.collect(), BENCHMARKING_TAXONOMY.out.stats_long.collect() )
    ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_TAXONOMY_PLOT.out.versions)

    BENCHMARKING_BARPLOT ( ch_barplot, ch_expected_barplot )
    ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_BARPLOT.out.versions)

    if ( params_benchmarking_diversity ) {
        expected = params_benchmarking_diversity.tokenize(',')
        ch_expected_table = Channel.fromPath( expected[0], checkIfExists: true)
        ch_expected_seq = Channel.fromPath( expected[1], checkIfExists: true)
        QIIME2_INASV_BM ( ch_expected_table )
        QIIME2_MERGE_TABLES ( QIIME2_INASV_BM.out.qza, ch_observed_table_qza )
        ch_benchmarking_versions = ch_benchmarking_versions.mix(QIIME2_MERGE_TABLES.out.versions)
        QIIME2_INSEQ_BM ( ch_expected_seq )
        QIIME2_MERGE_SEQ ( QIIME2_INSEQ_BM.out.qza, ch_observed_seq_qza )
        ch_benchmarking_versions = ch_benchmarking_versions.mix(QIIME2_MERGE_SEQ.out.versions)
        QIIME2_TREE_BM ( QIIME2_MERGE_SEQ.out.qza )
        BENCHMARKING_DIVERSITY ( QIIME2_MERGE_TABLES.out.qza, QIIME2_TREE_BM.out.qza, val_md5sum_version )
        ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_DIVERSITY.out.versions)
    }

    emit:
    versions       = ch_benchmarking_versions
}
