/*
 * Export filtered tables from QIIME2
 */

include { QIIME2_EXPORT_ABSOLUTE } from '../../modules/local/qiime2_export_absolute'
include { QIIME2_EXPORT_RELASV   } from '../../modules/local/qiime2_export_relasv'
include { QIIME2_EXPORT_RELTAX   } from '../../modules/local/qiime2_export_reltax'
include { COMBINE_TABLE          } from '../../modules/local/combine_table'

workflow QIIME2_EXPORT {
    take:
    ch_asv
    ch_seq
    ch_tax
    ch_tax_tsv
    tax_agglom_min
    tax_agglom_max

    main:
    ch_versions_qiime2_export = Channel.empty()

    //export_filtered_dada_output (optional)
    QIIME2_EXPORT_ABSOLUTE ( ch_asv, ch_seq, ch_tax, tax_agglom_min, tax_agglom_max )
    ch_versions_qiime2_export = ch_versions_qiime2_export.mix(QIIME2_EXPORT_ABSOLUTE.out.versions)

    //RelativeAbundanceASV (optional)
    QIIME2_EXPORT_RELASV ( ch_asv )
    ch_versions_qiime2_export = ch_versions_qiime2_export.mix(QIIME2_EXPORT_RELASV.out.versions)

    //RelativeAbundanceReducedTaxa (optional)
    QIIME2_EXPORT_RELTAX ( ch_asv, ch_tax, tax_agglom_min, tax_agglom_max )
    ch_versions_qiime2_export = ch_versions_qiime2_export.mix(QIIME2_EXPORT_RELTAX.out.versions)

    //combine_table.r (optional), similar to DADA2_table.tsv but with additionally taxonomy merged
    COMBINE_TABLE ( ch_tax_tsv.combine(QIIME2_EXPORT_RELASV.out.tsv).combine(QIIME2_EXPORT_ABSOLUTE.out.fasta) )
    ch_versions_qiime2_export = ch_versions_qiime2_export.mix(COMBINE_TABLE.out.versions)

    emit:
    comb_asv  = COMBINE_TABLE.out.tsv
    abs_fasta = QIIME2_EXPORT_ABSOLUTE.out.fasta
    abs_tsv   = QIIME2_EXPORT_ABSOLUTE.out.tsv
    rel_tsv   = QIIME2_EXPORT_RELASV.out.tsv
    versions  = ch_versions_qiime2_export
}
