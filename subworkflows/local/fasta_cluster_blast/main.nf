//
// Compare sequences by performing an all-v-all BLAST
//
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { BLAST_MAKEBLASTDB         } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN              } from '../../../modules/nf-core/blast/blastn/main'
include { ANICLUSTER_ANICALC        } from '../../../modules/local/anicluster/anicalc/main'
include { ANICLUSTER_ANICLUST       } from '../../../modules/local/anicluster/aniclust/main'
include { ANICLUSTER_EXTRACTREPS    } from '../../../modules/local/anicluster/extractreps/main'

workflow FASTA_CLUSTER_BLAST {

    take:
    fasta_gz    // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)
    min_ani     // val [ 0 - 100 ]          , minimum ANI for alignment to be counted
    min_qcov    // val [ 0 - 100 ]          , minimum query cover for clustering
    min_tcov    // val [ 0 - 100 ]          , minimum test cover for clustering

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

    //
    // MODULE: Calculate average nucleotide identity (ANI) and alignment fraction (AF) based on BLAST
    //
    ch_ani_tsv = ANICLUSTER_ANICALC ( ch_blast_txt ).ani
    ch_versions = ch_versions.mix( ANICLUSTER_ANICALC.out.versions )

    // create input for ANICLUSTER_ANICALC
    ch_aniclust_input = fasta_gz.join( ch_ani_tsv )

    //
    // MODULE: Cluster virus sequences based on ANI and AF
    //
    ch_clusters_tsv = ANICLUSTER_ANICLUST ( ch_aniclust_input.map { [ it[0], it[1] ] }, ch_aniclust_input.map { [ it[0], it[2] ] }, min_ani, min_qcov, min_tcov ).clusters
    ch_versions = ch_versions.mix( ANICLUSTER_ANICLUST.out.versions )

    // create input for extracting cluster representatives
    ch_extractreps_input = fasta_gz.join( ch_clusters_tsv )

    //
    // MODULE: Extract cluster representatives
    //
    ch_anicluster_reps_fasta_gz = ANICLUSTER_EXTRACTREPS ( ch_extractreps_input ).representatives
    ch_versions = ch_versions.mix( ANICLUSTER_EXTRACTREPS.out.versions )

    emit:
    cluster_reps_fasta_gz   = ch_anicluster_reps_fasta_gz   // [ [ meta ], cluster_reps.fasta.gz ]  , FASTA file cluster representatives
    clusters_tsv            = ch_clusters_tsv               // [ [ meta ], clusters.tsv          ]  , TSV file containing cluster membership
    versions                = ch_versions                   // [ versions.yml ]

}
