//
// Predict phage hosts using iPHoP
//

include { IPHOP_DOWNLOAD    } from '../../../modules/nf-core/iphop/download/main'
include { UNTAR             } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR2   } from '../../../modules/nf-core/untar/main'
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
        // if partial run is specified the minimal, incomplete database will be downloaded
        if ( !params.iphop_partial_test ) {
            //
            // MODULE: download standard IPHoP database
            //
            ch_iphop_db = IPHOP_DOWNLOAD( ).iphop_db
            ch_versions = ch_versions.mix(IPHOP_DOWNLOAD.out.versions)
        } else {
            //
            // MODULE: unpack minimal iphop db
            //
            ch_iphop_db = UNTAR ( [ [ id:'iphop_minimal_db' ], file(params.test_data['modules_nfcore']['iphop_test_db_tar_gz'], checkIfExists: true) ] ).untar.map { it[1] }
            ch_versions = ch_versions.mix(UNTAR.out.versions)
        }
    }

    //
    // MODULE: Gunzip fasta for IPHoP
    //
    ch_viruses_fasta = GUNZIP ( virus_fasta_gz ).gunzip
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())

    // create input for partial iphop run
    ch_iphop_partial_input = Channel.empty()
    if ( params.iphop_partial_test ) {
        //
        // MODULE: unpack data for partial iphop run
        //
        ch_iphop_partial_input = UNTAR2 ( [ [ id:'iphop_partial_input' ], file(params.test_data['modules_nfcore']['iphop_data_tar_gz'], checkIfExists: true) ] ).untar.map { it[1] }
        ch_versions = ch_versions.mix(UNTAR2.out.versions)
    }

    //
    // MODULE: Predict virus host
    //
    IPHOP_PREDICT ( ch_viruses_fasta, ch_iphop_db, ch_iphop_partial_input )
    ch_host_predictions_tsv  = IPHOP_PREDICT.out.iphop_genus
    ch_versions = ch_versions.mix(IPHOP_PREDICT.out.versions.first())

    emit:
    host_predictions_tsv    = ch_host_predictions_tsv  // [ [ meta ], genus_predictions.tsv ]   , TSV file containing host genus predictions
    versions                = ch_versions               // [ versions.yml ]
}
