//
// Predict provirus activity with Propagate
//
include { PROPAGATE_EXTRACTPROVIRUSES                           } from '../../../modules/local/propagate/extractproviruses/main'
include { CAT_CAT as CAT_FASTA                                  } from '../../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_COORDS                                 } from '../../../modules/nf-core/cat/cat/main'
include { GUNZIP                                                } from '../../../modules/nf-core/gunzip/main'
include { FASTA_CLUSTER_BLAST as FASTA_GROUP_DEREPLICATE_BLAST  } from '../fasta_cluster_blast/main'
include { PROPAGATE_DEREPCOORDINATES                            } from '../../../modules/local/propagate/derepcoordinates/main'
include { PROPAGATE_PROPAGATE                                   } from '../../../modules/local/propagate/propagate/main'

workflow FASTQ_FASTA_PROVIRUS_ACTIVITY_PROPAGATE {
    take:
    fastq_gz                // [ [ meta ], reads_1.fastq.gz, reads_1.fastq.gz ] , reads (mandatory)
    fasta_gz                // [ [ meta ], assembly.fasta.gz ]                  , assemblies (mandatory)
    ch_virus_summaries_tsv  // [ [ meta ], virus_summary.tsv ]                  , genomad virus summary (mandatory)
    ch_contamination_tsv    // [ [ meta ], contamination.tsv ]                  , checkV contamination summary (mandatory)
    min_ani                 // val [ 0 - 100 ]                                  , minimum ANI for dereplication
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

    // combine provirus assemblies by group
    ch_grouped_proviruses_fasta_gz = ch_provirus_scaffolds_fasta_gz
        .map {
            meta, fasta ->
                def meta_new = [:]
                meta_new.id = meta.group
                return [ meta_new, fasta ]
        }
        .groupTuple( sort: 'deep' )

    //
    // MODULE: Concatenate provirus assemblies within groups
    //
    ch_combined_proviruses_fasta_gz = CAT_FASTA ( ch_grouped_proviruses_fasta_gz ).file_out
    ch_versions                     = ch_versions.mix(CAT_FASTA.out.versions)

    //
    // SUBWORKFLOW: Dereplicate provirus-containing assemblies within groups
    //
    ch_derep_scaffolds_fasta_gz = FASTA_GROUP_DEREPLICATE_BLAST (
        ch_combined_proviruses_fasta_gz,
        min_ani,
        min_qcov,
        min_tcov
    ).cluster_reps_fasta_gz
    ch_derep_clusters_tsv       = FASTA_GROUP_DEREPLICATE_BLAST.out.clusters_tsv
    ch_versions                 = ch_versions.mix ( FASTA_GROUP_DEREPLICATE_BLAST.out.versions )

    //
    // MODULE: Gunzip assemblies for input into propagate
    //
    ch_derep_scaffolds_fasta    = GUNZIP ( ch_derep_scaffolds_fasta_gz ).gunzip
    ch_versions                 = ch_versions.mix ( GUNZIP.out.versions )

    // Combine coords files within groups
    ch_grouped_coords_tsv = ch_provirus_coords_tsv
        .map {
            meta, coords ->
                def meta_new = [:]
                meta_new.id = meta.group
                return [ meta_new, coords ]
        }
        .groupTuple( sort:'deep')

    //
    // MODULE: Combine coords files within group into one coords file
    //
    ch_combined_coords_tsv  = CAT_COORDS ( ch_grouped_coords_tsv ).file_out
    ch_versions             = ch_versions.mix ( CAT_COORDS.out.versions )

    // combine coords and cluster files by group
    ch_derep_coords_input = ch_combined_coords_tsv.join ( ch_derep_clusters_tsv )
        .multiMap {
            meta, coords, clusters ->
                coords: [ meta, coords ]
                clusters: [ meta, clusters ]
        }

    //
    // MODULE: Identify provirus coordinates in dereplicated assemblies
    //
    ch_derep_provirus_coords_tsv    = PROPAGATE_DEREPCOORDINATES ( ch_derep_coords_input.coords, ch_derep_coords_input.clusters ).derep_coords
    ch_versions                     = ch_versions.mix ( PROPAGATE_DEREPCOORDINATES.out.versions )

    // join reads with assemblies and provirus coords by group
    ch_propagate_input = fastq_gz
        .map {
            meta, fastq ->
                def meta_new = [:]

                meta_new.id = meta.group

                return [ meta_new, meta.id, fastq ]
        }
        .join ( ch_derep_scaffolds_fasta, by:0 )
        .join ( ch_derep_provirus_coords_tsv, by: 0 )
        .multiMap {
            meta, id, fastq, fasta, coords ->
                def meta_new = [:]

                meta_new.id     = id
                meta_new.group  = meta.id

                reads: [ meta_new, fastq ]
                assemblies: [ meta_new, fasta ]
                coords: [ meta_new, coords ]
        }

    //
    // MODULE: Predict provirus activity
    //
    ch_propagate_results_tsv    = PROPAGATE_PROPAGATE (
        ch_propagate_input.reads,
        ch_propagate_input.assemblies,
        ch_propagate_input.coords
        ).results
    ch_versions                 = ch_versions.mix(PROPAGATE_PROPAGATE.out.versions)

    emit:
    propagate_results_tsv   = ch_propagate_results_tsv  // [ [ meta ], propagate_results.tsv ]  , TSV file containing provirus activity predictions
    versions                = ch_versions               // [ versions.yml ]
}
