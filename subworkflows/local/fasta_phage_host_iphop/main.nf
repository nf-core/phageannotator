//
// Predict phage hosts using iPHoP
//

include { IPHOP_DOWNLOAD    } from '../../../modules/nf-core/iphop/download/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'
include { IPHOP_PREDICT     } from '../../../modules/nf-core/iphop/predict/main'

workflow FASTA_PHAGE_HOST_IPHOP {
    take:
    virus_fasta_gz  // [ [ meta ], fasta.gz ]   , virus sequences (mandatory)
    iphop_db        // [ checkv_db ]            , iPHoP database directory (optional)

    main:
    ch_versions = Channel.empty()

    // if iphop_db exists, skip IPHOP_DOWNLOAD
    if ( iphop_db ){
        ch_iphop_db = iphop_db
    } else {
        //
        // MODULE: download IPHoP database
        //
        ch_iphop_db = IPHOP_DOWNLOAD( ).iphop_db
        ch_versions = ch_versions.mix(IPHOP_DOWNLOAD.out.versions.first())
    }

    //
    // MODULE: Gunzip fasta for IPHoP
    //
    ch_viruses_fasta = GUNZIP ( virus_fasta_gz ).gunzip
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())

    //
    // MODULE: Predict virus host
    //
    IPHOP_PREDICT ( virus_fasta_gz, ch_iphop_db )
    ch_host_predictions_tsv  = IPHOP_PREDICT.out.iphop_genus
    ch_versions = ch_versions.mix(IPHOP_PREDICT.out.versions.first())

    emit:
    host_predictions_tsv    = ch_host_predictions_tsv  // [ [ meta ], genus_predictions.tsv ]   , TSV file containing host genus predictions
    versions                = ch_versions               // [ versions.yml ]
}
