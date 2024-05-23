include { BENCHMARKING_SEQUENCES     } from '../../modules/local/benchmarking_sequences'
include { BENCHMARKING_TAXONOMY      } from '../../modules/local/benchmarking_taxonomy'
include { BENCHMARKING_TAXONOMY_PLOT } from '../../modules/local/benchmarking_taxonomy_plot'
include { BENCHMARKING_BARPLOT       } from '../../modules/local/benchmarking_barplot'

workflow BENCHMARKING_WF {
    take:
    ch_asv
    ch_expected_sequences
    ch_tax
    ch_expected_taxonomies
    ch_barplot
    ch_expected_barplot

    main:
    ch_benchmarking_versions = Channel.empty()

    BENCHMARKING_SEQUENCES ( ch_asv, ch_expected_sequences )
    ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_SEQUENCES.out.versions)

    BENCHMARKING_TAXONOMY ( ch_tax.combine(ch_expected_taxonomies) )
    ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_TAXONOMY.out.versions)

    BENCHMARKING_TAXONOMY_PLOT ( BENCHMARKING_TAXONOMY.out.summary_tsv.collect(), BENCHMARKING_TAXONOMY.out.stats_long.collect() )
    ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_TAXONOMY_PLOT.out.versions)

    BENCHMARKING_BARPLOT ( ch_barplot, ch_expected_barplot )
    ch_benchmarking_versions = ch_benchmarking_versions.mix(BENCHMARKING_BARPLOT.out.versions)

    emit:
    versions       = ch_benchmarking_versions
}
