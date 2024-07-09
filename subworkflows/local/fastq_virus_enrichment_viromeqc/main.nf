//
// Estimate viral enrichment in reads
//

include { VIROMEQC_INSTALL  } from '../../../modules/local/viromeqc/install/main'
include { VIROMEQC_VIROMEQC } from '../../../modules/local/viromeqc/viromeqc/main'

workflow FASTQ_VIRUS_ENRICHMENT_VIROMEQC {
    take:
    fastq_gz                // [ [ meta.id ] , [ 1.fastq.gz, 2.fastq.gz ]  , reads (mandatory)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Install ViromeQC index
    //
    ch_viromeqc_index = VIROMEQC_INSTALL ( ).viromeqc_index
    ch_versions = ch_versions.mix(VIROMEQC_INSTALL.out.versions.first())

    //
    // MODULE: Estimate viral enrichment
    //
    ch_enrichment_tsv = VIROMEQC_VIROMEQC ( fastq_gz, ch_viromeqc_index ).enrichment
    ch_versions = ch_versions.mix(VIROMEQC_VIROMEQC.out.versions.first())

    emit:
    enrichment_tsv  = ch_enrichment_tsv // [ [ meta.id ] , enrichment.tsv ]  , viral enrichment estimates
    versions        = ch_versions       // [ versions.yml ]
}
