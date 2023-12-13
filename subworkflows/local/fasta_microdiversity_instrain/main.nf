//
// Assess virus microdiversity with instrain
//

include { GUNZIP as GUNZIP_FASTA    } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PROTEINS } from '../../../modules/nf-core/gunzip/main'
include { INSTRAIN_PROFILE          } from '../../../modules/nf-core/instrain/profile/main'
include { INSTRAIN_COMPARE          } from '../../../modules/nf-core/instrain/compare/main'


workflow FASTA_MICRODIVERSITY_INSTRAIN {
    take:
    bam             // [ [ meta ], bam ]        , BAM files from reads aligned to FASTA file (mandatory)
    fasta_gz        // [ [ meta ], fasta.gz ]   , FASTA file used in read alignment (mandatory)
    proteins_fna_gz // [ [ meta ], fna.gz ]     , FASTA file for protein-coding genes (optional)
    instrain_stb    // [ [ meta ], stb.tsv ]    , TSV file with two columns for associationg scaffolds to bins (optional)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: gunzip FASTA file
    //
    fasta = GUNZIP_FASTA ( fasta_gz ).gunzip

    //
    // MODULE: gunzip proteins fna file
    //
    proteins_fna = GUNZIP_PROTEINS ( proteins_fna_gz ).gunzip
    ch_versions = ch_versions.mix(GUNZIP_PROTEINS.out.versions)

    // remove meta information for instrain profile inputs
    ch_stb_file_tsv_nometa = instrain_stb.map { it[1] }.first()
    ch_fasta_nometa = fasta.map { it[1] }.first()
    ch_proteins_fna_nometa = proteins_fna.map { it[1] }.first()

    //
    // MODULE: Profile microdiveristy within each sample
    //
    ch_instrain_profiles = INSTRAIN_PROFILE ( bam, ch_fasta_nometa, ch_proteins_fna_nometa, ch_stb_file_tsv_nometa ).profile
    ch_instrain_gene_tsv = INSTRAIN_PROFILE.out.gene_info
    ch_versions = ch_versions.mix(INSTRAIN_PROFILE.out.versions)

    // combine bams and profiles across samples
    ch_bam_profiles_combined = bam.join(ch_instrain_profiles)
    ch_instrain_compare_input = ch_bam_profiles_combined.map { [ [ id:'all_samples' ], it[1], it[2] ] }.groupTuple( sort: 'deep' )

    //
    // MODULE: Compare microdiversity across samples (within groups OR across all samples)
    //
    ch_instrain_compare = INSTRAIN_COMPARE ( ch_instrain_compare_input, ch_stb_file_tsv_nometa ).compare
    ch_instrain_comparisons_tsv = INSTRAIN_COMPARE.out.comparisons_table
    ch_instrain_pooled_snvs_tsv = INSTRAIN_COMPARE.out.pooled_snv

    ch_versions = ch_versions.mix(INSTRAIN_COMPARE.out.versions)

    emit:
    gene_info_tsv   = ch_instrain_gene_tsv  // [ [ meta ], gene_info.tsv ]
    versions        = ch_versions           // [ versions.yml ]
}
