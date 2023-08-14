//
// De novo identification of viral sequences in assemblies
//

include { GENOMAD_DOWNLOAD   } from '../../../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND   } from '../../../modules/nf-core/genomad/endtoend/main'

workflow DE_NOVO_VIRUS_IDENTIFICATION {
    take:
    assemblies   // [ [ meta] , fasta    ], input assemblies (mandatory)

    main:
    ch_versions = Channel.empty()

    if ( params.genomad_db ){
        ch_genomad_db = file(params.genomad_db, checkIfExists: true)
    } else {
        ch_genomad_db = GENOMAD_DOWNLOAD ( ).genomad_db
        ch_versions = ch_versions.mix( GENOMAD_DOWNLOAD.out.versions )
    }

    GENOMAD_ENDTOEND ( assemblies, ch_genomad_db )
    ch_versions = ch_versions.mix( GENOMAD_ENDTOEND.out.versions )

    emit:
    identified_viruses = GENOMAD_ENDTOEND.out.virus_fasta   // [ [ meta ], fasta ]
    versions = ch_versions                                  // [ versions.yml ]

}
