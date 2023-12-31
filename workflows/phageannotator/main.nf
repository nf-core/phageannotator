/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: Update nf-core modules

//
// MODULES: Local modules
//
include { SEQKIT_SEQ                                } from '../../modules/local/seqkit/seq/main'                                    // TODO: Add to nf-core
include { AWK as AWK_GENOMAD                        } from '../../modules/local/awk/main'                                           // TODO: Add to nf-core
include { APPENDSCREENHITS                          } from '../../modules/local/appendscreenhits/main'
include { AWK as AWK_CHECKV                         } from '../../modules/local/awk/main'                                           // TODO: Add to nf-core
include { QUALITYFILTERVIRUSES                      } from '../../modules/local/qualityfilterviruses/main'
include { ANICLUSTER_ANICALC                        } from '../../modules/local/anicluster/anicalc/main'
include { ANICLUSTER_ANICLUST                       } from '../../modules/local/anicluster/aniclust/main'
include { ANICLUSTER_EXTRACTREPS                    } from '../../modules/local/anicluster/extractreps/main'
include { COVERM_CONTIG                             } from '../../modules/local/coverm/contig/main'                                 // TODO: Add to nf-core
include { PRODIGAL_PRODIGALGV                       } from '../../modules/local/prodigal/prodigalgv/main'                           // TODO: Add to nf-core
include { INSTRAIN_STB                              } from '../../modules/local/instrain/stb/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH                                } from '../../subworkflows/local/fastq_fasta_reference_containment_mash/main'   // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_VIRUS_CLASSIFICATION_GENOMAD                                    } from '../../subworkflows/local/fasta_virus_classification_genomad/main'       // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_VIRUS_QUALITY_CHECKV                                            } from '../../subworkflows/local/fasta_virus_quality_checkv/main'               // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_ALL_V_ALL_BLAST                                                 } from '../../subworkflows/local/fasta_all_v_all_blast/main'
include { FASTA_PHAGE_HOST_IPHOP                                                } from '../../subworkflows/local/fasta_phage_host_iphop/main'                   // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_MICRODIVERSITY_INSTRAIN                                         } from '../../subworkflows/local/fasta_microdiversity_instrain/main'            // TODO: Add to nf-core; Add nf-tests to nf-core modules


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_CAT as CAT_MASHSCREEN             } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_VIRUSES                } from '../../modules/nf-core/cat/cat/main'
include { BOWTIE2_BUILD                         } from '../../modules/nf-core/bowtie2/build/main'
include { GENOMAD_ENDTOEND as GENOMAD_TAXONOMY  } from '../../modules/nf-core/genomad/endtoend/main'
include { GUNZIP                                } from '../../modules/nf-core/gunzip/main'
include { BACPHLIP                              } from '../../modules/nf-core/bacphlip/main'

//
// SUBWORKFLOW: Installed directory from nf-core/subworkflows
//
include { FASTQ_ALIGN_BOWTIE2           } from '../../subworkflows/nf-core/fastq_align_bowtie2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PHAGEANNOTATOR {

    take:
    fastq_gz    // [ [ meta ], reads.fastq.gz ]     , reads (mandatory)
    fasta_gz    // [ [ meta ], assembly.fasta.gz ]  , assemblies/genomes (mandatory)

    main:
    ch_versions = Channel.empty()


    /*----------------------------------------------------------------------------
        Filter input assemblies
    ------------------------------------------------------------------------------*/
    //
    // MODULE: Filter assemblies by length
    //
    ch_filtered_input_fasta_gz = SEQKIT_SEQ ( fasta_gz ).fastx
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions.first())


    /*----------------------------------------------------------------------------
        OPTIONAL: Identify reference virus genomes contained in reads
    ------------------------------------------------------------------------------*/
    // if skip_reference_containment == false, run subworkflow
    if ( !params.skip_reference_containment ) {
        // if reference based identification requested, a reference FASTA file must be included
        if ( !params.reference_virus_fasta ) {
            error "[nf-core/phageannotator] ERROR: reference containment requested, but no --reference_virus_fasta provided"
        }

        // create channel from params.reference_virus_fasta
        ch_reference_virus_fasta_gz = Channel.of([ [ id:'reference_viruses' ], file( params.reference_virus_fasta, checkIfExists:true ) ])
        ch_reference_virus_fasta_gz

        // create channel from params.reference_virus_sketch
        if ( !params.reference_virus_sketch ){
            ch_reference_virus_sketch_msh = null
        } else {
            ch_reference_virus_sketch_msh = [ [ id:'reference_viruses' ], file( params.reference_virus_sketch, checkIfExists:true ) ]
        }

        //
        // SUBWORKFLOW: Identify contained reference genomes
        //
        ch_containment_results_tsv = FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH ( fastq_gz, ch_filtered_input_fasta_gz, ch_reference_virus_fasta_gz.first(), ch_reference_virus_sketch_msh ).mash_screen_results
        ch_versions = ch_versions.mix(FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH.out.versions.first())

        // join mash screen and assembly fasta by meta.id
        ch_append_screen_hits_input = ch_containment_results_tsv.join( ch_filtered_input_fasta_gz, by:0 )

        //
        // MODULE: Append screen hits to assemblies
        //
        ch_assembly_w_references_fasta_gz = APPENDSCREENHITS ( ch_append_screen_hits_input, ch_reference_virus_fasta_gz.first() ).assembly_w_screen_hits
        ch_versions = ch_versions.mix(APPENDSCREENHITS.out.versions.first())

        //
        // MODULE: Combine mash screen outputs across samples
        //
        ch_combined_mash_screen_tsv = CAT_MASHSCREEN( ch_containment_results_tsv.map{ [ [ id:'all_samples' ], it[1] ] }.groupTuple( sort: 'deep' ) ).file_out
        ch_versions = ch_versions.mix(CAT_MASHSCREEN.out.versions.first())
    } else {
        // if skip_reference_containment == true, skip subworkflow and use input assemblies
        ch_assembly_w_references_fasta_gz = ch_filtered_input_fasta_gz
        ch_combined_mash_screen_tsv = null
    }


    /*----------------------------------------------------------------------------
        Classify/annotate viral sequences
    ------------------------------------------------------------------------------*/
    // create channel from params.genomad_db
    if ( !params.genomad_db ){
        ch_genomad_db = null
    } else {
        ch_genomad_db = file( params.genomad_db, checkIfExists:true )
    }

    //
    // SUBWORKFLOW: Classify and annotate sequences
    //
    ch_viruses_fna_gz = FASTA_VIRUS_CLASSIFICATION_GENOMAD ( ch_assembly_w_references_fasta_gz, ch_genomad_db ).viruses_fna_gz
    ch_genomad_db_dir = FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.genomad_db
    ch_versions = ch_versions.mix(FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.versions.first())

    // create a channel for combining geNomad virus summaries (sorted so output is the same for tests)
    ch_awk_genomad_input = FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.virus_summaries_tsv
                            .map { [ [ id:'all_samples' ], it[1] ] }
                            .groupTuple(sort: 'deep')

    //
    // MODULE: Combine geNomad summaries across samples
    //
    ch_combined_virus_summaries_tsv = AWK_GENOMAD ( ch_awk_genomad_input ).file_out
    ch_versions = ch_versions.mix(AWK_GENOMAD.out.versions.first())


    /*----------------------------------------------------------------------------
        Assess virus quality and filter
    ------------------------------------------------------------------------------*/
    // create channel from params.checkv_db
    if ( !params.checkv_db ){
        ch_checkv_db = null
    } else {
        ch_checkv_db = file( params.checkv_db, checkIfExists:true )
    }

    //
    // SUBWORKFLOW: Assess virus quality
    //
    FASTA_VIRUS_QUALITY_CHECKV ( ch_viruses_fna_gz, ch_checkv_db )
    ch_versions = ch_versions.mix(FASTA_VIRUS_QUALITY_CHECKV.out.versions.first())

    // create a channel for combining Checkv quality summaries (sorted so output is the same for tests)
    ch_awk_checkv_input = FASTA_VIRUS_QUALITY_CHECKV.out.quality_summary_tsv
                            .map { [ [ id:'all_samples' ], it[1] ] }
                            .groupTuple( sort: 'deep' )

    //
    // MODULE: Combine quality summaries across samples
    //
    ch_combined_quality_summaries_tsv = AWK_CHECKV ( ch_awk_checkv_input ).file_out
    ch_versions = ch_versions.mix(AWK_CHECKV.out.versions.first())

    // create channel for input into QUALITY_FILTER_VIRUSES
    ch_quality_filter_viruses_input1 = FASTA_VIRUS_QUALITY_CHECKV.out.viruses_fna_gz.join(FASTA_VIRUS_QUALITY_CHECKV.out.proviruses_fna_gz)
    ch_quality_filter_viruses_input2 = ch_quality_filter_viruses_input1.join(FASTA_VIRUS_QUALITY_CHECKV.out.quality_summary_tsv)

    //
    // MODULE: Quality filter viruses
    //
    ch_filtered_viruses_fna_gz = QUALITYFILTERVIRUSES ( ch_quality_filter_viruses_input2 ).filtered_viruses
    ch_versions = ch_versions.mix(QUALITYFILTERVIRUSES.out.versions.first())


    /*----------------------------------------------------------------------------
        Cluster viruses using all-v-all BLAST approach
    ------------------------------------------------------------------------------*/
    // create a channel for combining filtered viruses (sorted so output is the same for tests)
    ch_cat_viruses_input = ch_filtered_viruses_fna_gz
                            .map { [ [ id:'all_samples' ], it[1] ] }
                            .groupTuple( sort: 'deep' )

    //
    // MODULE: Concatenate all quality filtered viruses into one file
    //
    ch_filtered_viruses_combined_fna_gz = CAT_VIRUSES ( ch_cat_viruses_input ).file_out

    //
    // SUBWORKFLOW: Perform all-v-all BLAST
    //
    ch_blast_txt = FASTA_ALL_V_ALL_BLAST ( ch_filtered_viruses_combined_fna_gz ).blast_txt
    ch_versions = ch_versions.mix( FASTA_ALL_V_ALL_BLAST.out.versions )

    //
    // MODULE: Calculate average nucleotide identity (ANI) and alignment fraction (AF) based on BLAST
    //
    ch_ani_tsv = ANICLUSTER_ANICALC ( ch_blast_txt ).ani
    ch_versions = ch_versions.mix( ANICLUSTER_ANICALC.out.versions )

    // create input for ANICLUSTER_ANICALC
    ch_aniclust_input = ch_filtered_viruses_combined_fna_gz.join( ch_ani_tsv )

    //
    // MODULE: Cluster virus sequences based on ANI and AF
    //
    ch_clusters_tsv = ANICLUSTER_ANICLUST ( ch_aniclust_input ).clusters
    ch_versions = ch_versions.mix( ANICLUSTER_ANICLUST.out.versions )

    // create input for extracting cluster representatives
    ch_extractreps_input = ch_filtered_viruses_combined_fna_gz.join( ch_clusters_tsv )

    //
    // MODULE: Extract cluster representatives
    //
    ch_anicluster_reps_fasta_gz = ANICLUSTER_EXTRACTREPS ( ch_extractreps_input ).representatives
    ch_versions = ch_versions.mix( ANICLUSTER_EXTRACTREPS.out.versions )


    /*----------------------------------------------------------------------------
        Align reads to viruses
    ------------------------------------------------------------------------------*/
    //
    // MODULE: Make bowtie2 index
    //
    ch_anicluster_reps_bt2 = BOWTIE2_BUILD ( ch_anicluster_reps_fasta_gz ).index.first()
    ch_versions = ch_versions.mix( BOWTIE2_BUILD.out.versions )

    //
    // SUBWORKFLOW: Align reads to bowtie2 index
    //
    ch_cluster_rep_alignment_bam = FASTQ_ALIGN_BOWTIE2 ( fastq_gz, ch_anicluster_reps_bt2, false, false, ch_anicluster_reps_fasta_gz ).bam
    ch_versions = ch_versions.mix( FASTQ_ALIGN_BOWTIE2.out.versions )

    //
    // MODULE: Calculate abundance metrics from BAM file
    //
    ch_combined_bams = ch_cluster_rep_alignment_bam.map { [ [ id:'all_samples' ], it[1] ] }.groupTuple( sort: 'deep' )
    ch_alignment_results_tsv = COVERM_CONTIG ( ch_combined_bams ).alignment_results
    ch_versions = ch_versions.mix( COVERM_CONTIG.out.versions )


    /*----------------------------------------------------------------------------
        Predict phage hosts
    ------------------------------------------------------------------------------*/
    // create channel from params.checkv_db
    if ( !params.iphop_db ){
        ch_iphop_db = null
    } else {
        ch_iphop_db = file( params.iphop_db, checkIfExists:true )
    }

    //
    // SUBWORKFLOW: Download database and predict phage hosts
    //
    ch_host_predictions_tsv = FASTA_PHAGE_HOST_IPHOP ( ch_anicluster_reps_fasta_gz, ch_iphop_db ).host_predictions_tsv
    ch_versions = ch_versions.mix( FASTA_PHAGE_HOST_IPHOP.out.versions )


    /*----------------------------------------------------------------------------
        Assign viral taxonomy
    ------------------------------------------------------------------------------*/
    //
    // SUBWORKFLOW: Assign taxonomy using ICTV taxa specific marker genes
    //
    ch_marker_taxonomy_tsv = GENOMAD_TAXONOMY ( ch_anicluster_reps_fasta_gz, ch_genomad_db_dir ).taxonomy
    ch_versions = ch_versions.mix( GENOMAD_TAXONOMY.out.versions )


    /*----------------------------------------------------------------------------
        Predict virus lifestyle
    ------------------------------------------------------------------------------*/
    // gunzip fasta for input into bacphlip
    ch_anicluster_reps_fasta = GUNZIP ( ch_anicluster_reps_fasta_gz ).gunzip
    ch_versions = ch_versions.mix( GUNZIP.out.versions )

    //
    // MODULE: Predict phage lifestyle using lysogeny specific genes
    //
    ch_bacphlip_lifestyle_tsv = BACPHLIP ( ch_anicluster_reps_fasta ).bacphlip_results
    ch_versions = ch_versions.mix( BACPHLIP.out.versions )


    /*----------------------------------------------------------------------------
        Identify protein-coding regions
    ------------------------------------------------------------------------------*/
    ch_prodigalgv_proteins_fna_gz = PRODIGAL_PRODIGALGV ( ch_anicluster_reps_fasta ).fna
    ch_prodigalgv_proteins_faa_gz = PRODIGAL_PRODIGALGV.out.faa
    ch_versions = ch_versions.mix( PRODIGAL_PRODIGALGV.out.versions )


    /*----------------------------------------------------------------------------
        Analyze virus microdiversity
    ------------------------------------------------------------------------------*/
    //
    // MODULE: Generate instrain scaffold to bin file
    //
    ch_stb_file_tsv  = INSTRAIN_STB ( ch_anicluster_reps_fasta ).stb
    ch_versions = ch_versions.mix(INSTRAIN_STB.out.versions)

    //
    // SUBWORKFLOW: Assess virus microdiversity within and across samples
    //
    ch_gene_info_tsv = FASTA_MICRODIVERSITY_INSTRAIN ( ch_cluster_rep_alignment_bam, ch_anicluster_reps_fasta_gz, ch_prodigalgv_proteins_fna_gz, ch_stb_file_tsv ).gene_info_tsv
    ch_versions = ch_versions = ch_versions.mix(FASTA_MICRODIVERSITY_INSTRAIN.out.versions)


    emit:
    reference_containment_tsv   = ch_combined_mash_screen_tsv
    virus_classification_tsv    = ch_combined_virus_summaries_tsv
    virus_quality_tsv           = ch_combined_quality_summaries_tsv
    filtered_viruses_fna_gz     = ch_filtered_viruses_fna_gz
    anicluster_reps_fna_gz      = ch_anicluster_reps_fasta_gz
    // alignment_results_tsv       = ch_alignment_results_tsv // Inconsistent hash
    host_predictions_tsv        = ch_host_predictions_tsv
    marker_taxonomy_tsv         = ch_marker_taxonomy_tsv
    // bacphlip_lifestyle_tsv      = ch_bacphlip_lifestyle_tsv // Inconsistent hash
    prodigalgv_proteins_faa_gz  = ch_prodigalgv_proteins_faa_gz
    instrain_gene_info          = ch_gene_info_tsv
    versions                    = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
