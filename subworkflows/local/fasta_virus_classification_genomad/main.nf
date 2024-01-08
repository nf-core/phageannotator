//
// Classify and annotate sequences with geNomad
//

include { GENOMAD_DOWNLOAD   } from '../../../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND   } from '../../../modules/nf-core/genomad/endtoend/main'    // TODO: Update nf-core module to gzip output files

workflow FASTA_VIRUS_CLASSIFICATION_GENOMAD {
    take:
    fasta_gz    // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)
    genomad_db  // [ genomad_db ]           , genomad database directory (optional)

    main:
    ch_versions = Channel.empty()

    // if genomad_db exists, skip GENOMAD_DOWNLOAD
    if ( genomad_db ){
        ch_genomad_db = genomad_db
    } else {
        //
        // MODULE: download geNomad database
        //
        ch_genomad_db = GENOMAD_DOWNLOAD( ).genomad_db
        ch_versions = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions.first())
    }

    //
    // MODULE: Classify/annotate viral sequences
    //
    ch_viruses_fna_gz = GENOMAD_ENDTOEND ( fasta_gz, ch_genomad_db ).virus_fasta
    ch_virus_summaries_tsv = GENOMAD_ENDTOEND.out.virus_summary
    ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions.first())

    emit:
    genomad_db          = ch_genomad_db             // [ genomad_db/ ]                  , directory containing genomad_db files
    viruses_fna_gz      = ch_viruses_fna_gz         // [ [ meta ], fna.gz ]             , FASTA file containing viral sequences
    virus_summaries_tsv = ch_virus_summaries_tsv    // [ [ meta ], virus_summary.tsv ]  , TSV file containing virus information
    versions            = ch_versions               // [ versions.yml ]

}
