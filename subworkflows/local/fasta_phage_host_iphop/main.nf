//
// Predict phage hosts using iPHoP
//

include { IPHOP_DOWNLOAD    } from '../../../modules/nf-core/iphop/download/main'
include { UNTAR             } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR2   } from '../../../modules/nf-core/untar/main'
include { IPHOP_PREDICT     } from '../../../modules/nf-core/iphop/predict/main'

workflow FASTA_PHAGE_HOST_IPHOP {
    take:
    virus_fasta     // [ [ meta ], fasta ]  , virus sequences (mandatory)
    iphop_db        // [ checkv_db ]        , iPHoP database directory (optional)

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
            ch_iphop_db = UNTAR ( [ [ id:'iphop_minimal_db' ], file("https://github.com/nf-core/test-datasets/raw/phageannotator/modules/nfcore/iphop/download/Test_db_rw.tar.gz", checkIfExists: true) ] ).untar.map { it[1] }
            ch_versions = ch_versions.mix(UNTAR.out.versions)
        }
    }

    // create input for partial iphop run
    ch_iphop_partial_input = []
    if ( params.iphop_partial_test ) {
        //
        // MODULE: unpack data for partial iphop run
        //
        ch_iphop_partial_input = UNTAR2 ( [ [ id:'iphop_partial_input' ], file("https://github.com/nf-core/test-datasets/raw/phageannotator/modules/nfcore/iphop/predict/iPHoP_data.tar.gz", checkIfExists: true) ] ).untar.map { it[1] }
        ch_versions = ch_versions.mix(UNTAR2.out.versions)
    }

    //
    // MODULE: Predict virus host
    //
    ch_host_predictions_tsv = IPHOP_PREDICT ( virus_fasta, ch_iphop_db, ch_iphop_partial_input ).iphop_genus
    ch_versions = ch_versions.mix(IPHOP_PREDICT.out.versions.first())

    emit:
    host_predictions_tsv    = ch_host_predictions_tsv   // [ [ meta ], genus_predictions.tsv ]   , TSV file containing host genus predictions
    versions                = ch_versions               // [ versions.yml ]
}
