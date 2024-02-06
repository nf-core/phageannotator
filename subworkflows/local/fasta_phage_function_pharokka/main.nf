//
// Functionally annotate phage sequences using Pharokka
//

include { PHAROKKA_INSTALLDATABASES } from '../../../modules/nf-core/pharokka/installdatabases/main'
include { PHAROKKA_PHAROKKA         } from '../../../modules/nf-core/pharokka/pharokka/main'

workflow FASTA_PHAGE_FUNCTION_PHAROKKA {
    take:
    virus_fasta     // [ [ meta ], fasta ]  , virus sequences (mandatory)
    pharokka_db     // [ pharokka_db ]      , Pharokka database directory (optional)

    main:
    ch_versions = Channel.empty()

    // if pharokka_db exists, skip database download
    if ( pharokka_db ){
        ch_pharokka_db = pharokka_db
    } else {
        //
        // MODULE: download standard pharokka database
        //
        ch_pharokka_db = PHAROKKA_INSTALLDATABASES( ).pharokka_db
        ch_versions = ch_versions.mix(PHAROKKA_INSTALLDATABASES.out.versions)
    }

    //
    // MODULE: Functionally annotate phages
    //
    ch_pharokka_gbk_gz = PHAROKKA_PHAROKKA ( virus_fasta, ch_pharokka_db ).gbk
    pharokka_final_output_tsv = PHAROKKA_PHAROKKA.out.cds_final_merged_output
    ch_versions = ch_versions.mix(PHAROKKA_PHAROKKA.out.versions.first())

    emit:
    versions                    = ch_versions               // [ versions.yml ]
    pharokka_gbk_gz             = ch_pharokka_gbk_gz        // [ [ meta ], gbk.gz ]
    pharokka_final_output_tsv   = pharokka_final_output_tsv // [ [ meta ], tsv ]
}
