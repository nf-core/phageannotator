//
// Predict provirus activity with Propagate
//
include { PROPAGATE_EXTRACTPROVIRUSES                           } from '../../../modules/local/propagate/extractproviruses/main'
include { CAT_CAT as CAT_FASTA                                  } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_COORDS                                 } from '../../modules/nf-core/cat/cat/main'
include { FASTA_CLUSTER_BLAST as FASTA_GROUP_DEREPLICATE_BLAST  } from '../fasta_cluster_blast/main'
include { PROPAGATE_DEREPCOORDINATES                            } from '../../../modules/local/propagate/derepcoordinates/main'
include { PROPAGATE_PROPAGATE                                   } from '../../../modules/local/propagate/propagate/main'

workflow FASTQ_FASTA_PROVIRUS_ACTIVITY_PROPAGATE {
    take:
    fastq_gz                // [ [ meta ], reads_1.fastq.gz, reads_1.fastq.gz ] , reads (mandatory)
    fasta_gz                // [ [ meta ], assembly.fasta.gz ]                  , assemblies (mandatory)
    ch_virus_summaries_tsv  // [ [ meta ], virus_summary.tsv ]                  , genomad virus summary (mandatory)
    ch_contamination_tsv    // [ [ meta ], contamination.tsv ]                  , checkV contamination summary (mandatory)
    min_ani                 // val [ 0 - 100 ]                                  , minimum ANI for blast hits to be counted
    min_qcov                // val [ 0 - 100 ]                                  , minimum query coverage when clustering
    min_tcov                // val [ 0 - 100 ]                                  , minimum test coverage when clustering

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Identify integrated proviruses (and optionally assign to clusters)
    //
    ch_provirus_scaffolds_fasta_gz  = PROPAGATE_EXTRACTPROVIRUSES ( fasta_gz, ch_virus_summaries_tsv, ch_contamination_tsv ).provirus_scaffolds
    ch_provirus_coords_tsv          = PROPAGATE_EXTRACTPROVIRUSES.out.provirus_coords
    ch_versions                     = ch_versions.mix(PROPAGATE_EXTRACTPROVIRUSES.out.versions)

    // organize proviral contigs by group
    ch_grouped_proviruses_fasta_gz = fasta_gz
        .map {
            meta, fasta ->
                def group = meta.group

                return [ group, fasta ]
        }
        .groupTuple( sort: 'deep' )

    //
    // MODULE: Concatenate all proviral contigs from the same group
    //
    ch_combined_proviruses_fasta_gz = CAT_CAT ( ch_grouped_proviruses_fasta_gz ).file_out
    ch_versions                     = ch_versions.mix(CAT_CAT.out.versions)

    //
    // SUBWORKFLOW: Cluster (Dereplicate) contigs containing integrated proviruses within groups
    //
    ch_derep_scaffolds_fasta_gz = FASTA_GROUP_DEREPLICATE_BLAST ( ch_combined_proviruses_fasta_gz, min_ani, min_qcov, min_tcov ).cluster_reps_fasta_gz
    ch_derep_clusters_tsv       = FASTA_GROUP_DEREPLICATE_BLAST.out.clusters_tsv
    ch_versions                 = ch_versions.mix( FASTA_GROUP_DEREPLICATE_BLAST.out.versions )

    // collect coords file by group into one coords file
    ch_grouped_coords_tsv = ch_provirus_coords_tsv
        .map {
            meta, coords ->
                return [ meta.group, coords ]
        }
        .groupTuple( sort:'deep')

    //
    // MODULE: Combine coords files within group into one coords file
    //
    ch_combined_coords_tsv  = CAT_COORDS ( ch_grouped_coords_tsv ).file_out
    ch_versions             = ch_versions.mix(CAT_CAT.out.versions)

    // combine coords and cluster files by group
    ch_provirus_coords_input = ch_combined_coords_tsv.combine ( ch_derep_clusters_tsv )

    //
    // MODULE: Create a coordinates file for dereplicated proviruses
    //
    ch_derep_provirus_coords_tsv    = PROPAGATE_DEREPCOORDINATES ( ch_provirus_coords_input.map { it[0], it[1] }, ch_provirus_coords_input.map { it[0], it[2] } ).coords
    ch_versions                     = ch_versions.mix(PROPAGATE_DEREPCOORDINATES.out.versions)

    // organize reads by group and combine with dereplicated contigs & dereplicated coords
    ch_propagate_input = fastq_gz
        .map {
            meta, fastq ->
                def id      = meta.id
                def group   = meta.group

                return [ group, id, fastq ]
        }
        .combine ( ch_dereplicated_proviruses_fasta_gz, by:0 )
        .combine ( ch_derep_provirus_coords_tsv, by: 0 )
        .map {
            group, id, fastq, fasta, coords ->
                def meta = [:]

                meta.id     = id
                meta.group  = group

                return [ meta, fastq, fasta, coords ]
        }

    //
    // MODULE: Predict provirus activity
    //
    ch_propagate_results_tsv    = PROPAGATE_PROPAGATE (
        ch_propagate_input.map { it[0], it[1] },
        ch_propagate_input.map { it[0], it[2] },
        ch_propagate_input.map { it[0], it[3] }
        ).propagate_results
    ch_versions                 = ch_versions.mix(PROPAGATE_PROPAGATE.out.versions)

    emit:
    propagate_results_tsv   = ch_propagate_results_tsv  // [ [ meta ], propagate_results.tsv ]  , TSV file containing provirus activity predictions
    versions                = ch_versions               // [ versions.yml ]
}
