//
// Compare sequences by performing an all-v-all BLAST
//
include { BLAST_MAKEBLASTDB } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN      } from '../../../modules/nf-core/blast/blastn/main'

workflow FASTA_ALL_V_ALL_BLAST {

    take:
    fasta_gz    // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Make BLASTN database
    //
    ch_blast_db = BLAST_MAKEBLASTDB ( fasta_gz ).db
    ch_versions = ch_versions.mix( BLAST_MAKEBLASTDB.out.versions )

    //
    // MODULE: Perform BLAST
    //
    ch_blast_txt = BLAST_BLASTN ( fasta_gz , ch_blast_db ).txt
    ch_versions = ch_versions = ch_versions.mix( BLAST_MAKEBLASTDB.out.versions )

    emit:
    blast_txt   = ch_blast_txt  // [ [ meta ], blast_output.tsv ]   , TSV file containing BLAST results
    versions    = ch_versions   // [ versions.yml ]

}
