//
// Assess virus microdiversity with instrain
//

include { INSTRAIN_PROFILE          } from '../../../modules/nf-core/instrain/profile/main'
include { INSTRAIN_COMPARE          } from '../../../modules/nf-core/instrain/compare/main'


workflow FASTA_MICRODIVERSITY_INSTRAIN {
    take:
    bam             // [ [ meta ], bam ]        , BAM files from reads aligned to FASTA file (mandatory)
    genome_fasta    // [ [ meta ], fasta ]      , FASTA file used in read alignment (mandatory)
    proteins_fna    // [ [ meta ], fna ]        , FASTA file for protein-coding genes (optional)
    instrain_stb    // [ [ meta ], stb.tsv ]    , TSV file with two columns for associationg scaffolds to bins (optional)

    main:
    ch_versions = Channel.empty()

    // remove meta information for instrain profile inputs
    ch_genome_fasta_nometa = genome_fasta.map { it[1] }.first()
    if ( proteins_fna ){
        ch_proteins_fna_nometa = proteins_fna.map { it[1] }.first()
    } else {
        ch_proteins_fna_nometa = []
    }
    if ( instrain_stb ){
        ch_stb_file_tsv_nometa = instrain_stb.map { it[1] }.first()
    } else {
        ch_stb_file_tsv_nometa = []
    }

    //
    // MODULE: Profile microdiveristy within each sample
    //
    ch_instrain_profiles = INSTRAIN_PROFILE ( bam, ch_genome_fasta_nometa, ch_proteins_fna_nometa, ch_stb_file_tsv_nometa ).profile
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
