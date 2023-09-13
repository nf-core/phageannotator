//
// Classify and annotate sequences with geNomad
//

include { GENOMAD_DOWNLOAD   } from '../../../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND   } from '../../../modules/nf-core/genomad/endtoend/main'
include { AWK as AWK_GENOMAD } from '../../../modules/local/awk/main'

workflow FASTA_VIRUS_CLASSIFICATION_GENOMAD {
    take:
    fasta_gz    // [ [ meta] , fasta.gz ], assemblies/genomes (mandatory)

    main:
    ch_versions = Channel.empty()

    // if genomad_db != null, skip GENOMAD_DOWNLOAD
    if ( params.genomad_db ){
        ch_genomad_db = file(params.genomad_db, checkIfExists: true)
    }
    // if genomad_db == null, run GENOMAD_DOWNLOAD
    else {
        //
        // MODULE: download geNomad database
        //
        ch_genomad_db = GENOMAD_DOWNLOAD ( ).genomad_db
        ch_versions = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions.first())
    }

    //
    // MODULE: Classify/annotate viral sequences
    //
    ch_viruses_fasta_gz = GENOMAD_ENDTOEND ( fasta_gz, ch_genomad_db ).virus_fasta
    ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions.first())

    //
    // MODULE: Combine geNomad summaries across samples
    //
    ch_virus_summaries_tsv = AWK_GENOMAD ( GENOMAD_ENDTOEND.out.virus_summary.map { [ [ id:'all_samples' ], it[1] ] } ).file_out
    ch_versions = ch_versions.mix(AWK_GENOMAD.out.versions.first())

    emit:
    viruses_fasta_gz = ch_viruses_fasta_gz  // [ [ meta ], fasta.gz ]   , FASTA file containing viral sequences
    versions = ch_versions                  // [ versions.yml ]

}
