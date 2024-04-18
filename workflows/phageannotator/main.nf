/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: Update nf-core modules versions

//
// MODULES: Local modules
//

include { SEQKIT_SEQ                                } from '../../modules/local/seqkit/seq/main'                                    // TODO: Add to nf-core
include { APPENDSCREENHITS                          } from '../../modules/local/appendscreenhits/main'
include { EXTRACTVIRALASSEMBLIES                    } from '../../modules/local/extractviralassemblies/main'
include { QUALITYFILTERVIRUSES                      } from '../../modules/local/qualityfilterviruses/main'
include { ANICLUSTER_ANICALC                        } from '../../modules/local/anicluster/anicalc/main'
include { ANICLUSTER_ANICLUST                       } from '../../modules/local/anicluster/aniclust/main'
include { ANICLUSTER_EXTRACTREPS                    } from '../../modules/local/anicluster/extractreps/main'
include { COVERM_CONTIG                             } from '../../modules/local/coverm/contig/main'                                 // TODO: Add to nf-core
include { INSTRAIN_STB                              } from '../../modules/local/instrain/stb/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FASTQ_VIRUS_ENRICHMENT_VIROMEQC           } from '../../subworkflows/local/fastq_virus_enrichment_viromeqc/main'
include { FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH    } from '../../subworkflows/local/fastq_fasta_reference_containment_mash/main'   // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_VIRUS_CLASSIFICATION_GENOMAD        } from '../../subworkflows/local/fasta_virus_classification_genomad/main'       // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTQ_FASTA_CONTIG_EXTENSION_COBRA        } from '../../subworkflows/local/fastq_fasta_contig_extension_cobra/main'       // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_VIRUS_QUALITY_CHECKV                } from '../../subworkflows/local/fasta_virus_quality_checkv/main'               // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_ALL_V_ALL_BLAST                     } from '../../subworkflows/local/fasta_all_v_all_blast/main'
include { FASTA_PHAGE_HOST_IPHOP                    } from '../../subworkflows/local/fasta_phage_host_iphop/main'                   // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_PHAGE_FUNCTION_PHAROKKA             } from '../../subworkflows/local/fasta_phage_function_pharokka/main'
include { FASTA_MICRODIVERSITY_INSTRAIN             } from '../../subworkflows/local/fasta_microdiversity_instrain/main'            // TODO: Add to nf-core; Add nf-tests to nf-core modules


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_CAT as CAT_VIRUSES                } from '../../modules/nf-core/cat/cat/main'
include { BOWTIE2_BUILD                         } from '../../modules/nf-core/bowtie2/build/main'
include { GENOMAD_ENDTOEND as GENOMAD_TAXONOMY  } from '../../modules/nf-core/genomad/endtoend/main'
include { GUNZIP as GUNZIP_CLUSTER_REPS         } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_VIRUS_PROTEINS       } from '../../modules/nf-core/gunzip/main'
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

workflow PHAGEANNOTATOR {

    take:
    fastq_gz    // [ [ meta ], reads.fastq.gz ]     , reads (mandatory)
    fasta_gz    // [ [ meta ], assembly.fasta.gz ]  , assemblies/genomes (mandatory)

    main:
    ch_versions         = Channel.empty()


    /*----------------------------------------------------------------------------
        Estimate viral enrichment in reads
    ------------------------------------------------------------------------------*/
    // if run_viromeqc == true, run subworkflow
    if ( params.run_viromeqc ) {
        ch_virus_enrichment_tsv = FASTQ_VIRUS_ENRICHMENT_VIROMEQC ( fastq_gz ).enrichment_tsv
        ch_versions = ch_versions.mix(FASTQ_VIRUS_ENRICHMENT_VIROMEQC.out.versions)
    } else {
        // if run_viromeqc == false, skip subworkflow
        ch_virus_enrichment_tsv = Channel.empty()
    }


    //
    // MODULE: Filter assemblies by length
    //
    ch_filtered_input_fasta_gz = SEQKIT_SEQ ( fasta_gz ).fastx
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)


    /*----------------------------------------------------------------------------
        Reference virus identification
    ------------------------------------------------------------------------------*/
    // if run_reference_containment == true, run subworkflow
    if ( params.run_reference_containment ) {
        // if reference based identification requested, a reference FASTA file must be included
        if ( !params.reference_virus_fasta ) {
            error "[nf-core/phageannotator] ERROR: reference containment requested, but no --reference_virus_fasta provided"
        }

        // create channel from params.reference_virus_fasta
        ch_reference_virus_fasta_gz = Channel.value([ [ id:'reference_viruses' ], file( params.reference_virus_fasta, checkIfExists:true ) ])

        // create channel from params.reference_virus_sketch
        if ( !params.reference_virus_sketch ){
            ch_reference_virus_sketch_msh = null
        } else {
            ch_reference_virus_sketch_msh = Channel.value( [ [ id:'reference_viruses' ], file( params.reference_virus_sketch, checkIfExists:true ) ] )
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
    } else {
        // if run_reference_containment == false, skip subworkflow and use input assemblies
        ch_assembly_w_references_fasta_gz = ch_filtered_input_fasta_gz
    }


    /*----------------------------------------------------------------------------
        De novo virus classification
    ------------------------------------------------------------------------------*/
    // if skip_genomad == false OR run_genomad_taxonomy == true, prepare genomad db
    if ( !params.skip_genomad || params.run_genomad_taxonomy ) {
        //
        // De-novo virus classification using assemblies
        //
        // create channel from params.genomad_db
        if ( !params.genomad_db ){
            ch_genomad_db = null
        } else {
            ch_genomad_db = Channel.value( file( params.genomad_db, checkIfExists:true ) )
        }

        // if skip_genomad == false run genomad subworkflow
        //
        // SUBWORKFLOW: Classify and annotate sequences
        //
        if ( !params.skip_genomad ){
            ch_viruses_fna_gz = FASTA_VIRUS_CLASSIFICATION_GENOMAD ( ch_assembly_w_references_fasta_gz, ch_genomad_db ).viruses_fna_gz
            ch_genomad_db_dir = FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.genomad_db
            ch_versions = ch_versions.mix(FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.versions.first())
            ch_virus_summaries_tsv = FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.virus_summaries_tsv
        } else {
            ch_genomad_db_dir = FASTA_VIRUS_CLASSIFICATION_GENOMAD ( Channel.empty(), ch_genomad_db ).genomad_db
            ch_versions = ch_versions.mix(FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.versions.first())
            ch_viruses_fna_gz = ch_assembly_w_references_fasta_gz
            ch_virus_summaries_tsv = Channel.empty()
        }
    } else {
        // if skip_genomad == true use assemblies with references
        ch_viruses_fna_gz = ch_assembly_w_references_fasta_gz
        ch_virus_summaries_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Extend viral contigs
    ------------------------------------------------------------------------------*/
    if ( params.run_cobra ) {
        //
        // MODULE: Create a TSV file containing viral contig names (from assemblies)
        //
        ch_viral_assemblies_tsv = EXTRACTVIRALASSEMBLIES ( ch_viruses_fna_gz ).viral_assemblies
        ch_versions = ch_versions.mix(EXTRACTVIRALASSEMBLIES.out.versions)

        //
        // SUBWORKFLOW: Extend assembled contigs
        //
        ch_extended_viruses_fasta_gz = FASTQ_FASTA_CONTIG_EXTENSION_COBRA (
            fastq_gz,
            fasta_gz,
            ch_viral_assemblies_tsv,
            params.cobra_assembler,
            params.cobra_mink,
            params.cobra_maxk
        ).extended_fasta
        ch_virus_extension_summary_tsv = FASTQ_FASTA_CONTIG_EXTENSION_COBRA.out.cobra_summary_tsv
        ch_versions = ch_versions.mix( FASTQ_FASTA_CONTIG_EXTENSION_COBRA.out.versions )
    } else {
        ch_extended_viruses_fasta_gz = ch_viruses_fna_gz
        ch_virus_extension_summary_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Assess virus quality and filter
    ------------------------------------------------------------------------------*/
    // if skip_checkv == false, run subworkflow
    if ( !params.skip_checkv ) {
        // create channel from params.checkv_db
        if ( !params.checkv_db ){
            ch_checkv_db = null
        } else {
            ch_checkv_db = Channel.value( file( params.checkv_db, checkIfExists:true ) )
        }

        //
        // SUBWORKFLOW: Assess virus quality
        //
        ch_quality_summary_tsv = FASTA_VIRUS_QUALITY_CHECKV ( ch_viruses_fna_gz, ch_checkv_db ).quality_summary_tsv
        ch_versions = ch_versions.mix(FASTA_VIRUS_QUALITY_CHECKV.out.versions.first())

        // create channel for input into QUALITY_FILTER_VIRUSES
        ch_quality_filter_viruses_input1 = FASTA_VIRUS_QUALITY_CHECKV.out.viruses_fna_gz.join(FASTA_VIRUS_QUALITY_CHECKV.out.proviruses_fna_gz)
        ch_quality_filter_viruses_input2 = ch_quality_filter_viruses_input1.join(FASTA_VIRUS_QUALITY_CHECKV.out.quality_summary_tsv)

        //
        // MODULE: Quality filter viruses
        //
        ch_filtered_viruses_fna_gz = QUALITYFILTERVIRUSES ( ch_quality_filter_viruses_input2 ).filtered_viruses
        ch_versions = ch_versions.mix(QUALITYFILTERVIRUSES.out.versions.first())
    } else {
        // if skip_checkv == false, use non-quality filtered viruses
        ch_filtered_viruses_fna_gz = ch_viruses_fna_gz
        ch_quality_summary_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Cluster viruses using all-v-all BLAST approach
    ------------------------------------------------------------------------------*/
    // if skip_virus_clustering == false, run subworkflow
    if ( !params.skip_virus_clustering  ) {
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
    } else {
        // if skip_virus_clustering == true, use unclustered viruses
        ch_anicluster_reps_fasta_gz = ch_filtered_viruses_fna_gz
        ch_clusters_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Align reads to viruses
    ------------------------------------------------------------------------------*/
    // if skip_read_alignment == false OR run_instrain == true, run subworkflow
    if ( !params.skip_read_alignment || params.run_instrain ) {
        //
        // MODULE: Make bowtie2 index
        //
        ch_anicluster_reps_bt2 = BOWTIE2_BUILD ( ch_anicluster_reps_fasta_gz ).index.first()
        ch_versions = ch_versions.mix( BOWTIE2_BUILD.out.versions )

        //
        // SUBWORKFLOW: Align reads to bowtie2 index
        //
        ch_cluster_rep_alignment_bam = FASTQ_ALIGN_BOWTIE2 ( fastq_gz, ch_anicluster_reps_bt2, false, false, ch_anicluster_reps_fasta_gz.first() ).bam
        ch_versions = ch_versions.mix( FASTQ_ALIGN_BOWTIE2.out.versions )

        //
        // MODULE: Calculate abundance metrics from BAM file
        //
        ch_combined_bams = ch_cluster_rep_alignment_bam.map { [ [ id:'all_samples' ], it[1] ] }.groupTuple( sort: 'deep' )
        ch_alignment_results_tsv = COVERM_CONTIG ( ch_combined_bams ).alignment_results
        ch_versions = ch_versions.mix( COVERM_CONTIG.out.versions )
    } else {
        // if skip_read_alignment == true AND run_instrain == false, skip subworkflow
        ch_alignment_results_tsv = Channel.empty()
        if ( !params.skip_instrain ) {
            error "[nf-core/phageannotator] ERROR: skip_read_alignment = true but skip_instrain = false; read alignment must take place for inStrain to run"
        }
    }


    /*----------------------------------------------------------------------------
        Assign viral taxonomy
    ------------------------------------------------------------------------------*/
    // if run_genomad_taxonomy == true run subworkflow
    if ( params.run_genomad_taxonomy ) {
        //
        // SUBWORKFLOW: Assign taxonomy using ICTV taxa specific marker genes
        //
        ch_marker_taxonomy_tsv = GENOMAD_TAXONOMY ( ch_anicluster_reps_fasta_gz, ch_genomad_db_dir ).taxonomy
        ch_versions = ch_versions.mix( GENOMAD_TAXONOMY.out.versions )
    } else {
        // if run_genomad_taxonomy == false, skip subworkflow
        ch_marker_taxonomy_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Predict phage hosts
    ------------------------------------------------------------------------------*/
    // gunzip fasta for input into iphop, bacphlip, pharokka,
    ch_anicluster_reps_fasta = GUNZIP_CLUSTER_REPS ( ch_anicluster_reps_fasta_gz ).gunzip
    ch_versions = ch_versions.mix( GUNZIP_CLUSTER_REPS.out.versions )

    // if run_iphop == true, run subworkflow
    if ( params.run_iphop ){
        // create channel from params.checkv_db
        if ( !params.iphop_db ){
            ch_iphop_db = null
        } else {
            ch_iphop_db = file( params.iphop_db, checkIfExists:true )
        }

        //
        // SUBWORKFLOW: Download database and predict phage hosts
        //
        ch_host_predictions_tsv = FASTA_PHAGE_HOST_IPHOP ( ch_anicluster_reps_fasta, ch_iphop_db ).host_predictions_tsv
        ch_versions = ch_versions.mix( FASTA_PHAGE_HOST_IPHOP.out.versions )
    } else {
        // if run_iphop == false, skip subworkflow
        ch_host_predictions_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Predict virus lifestyle
    ------------------------------------------------------------------------------*/
    // if run_bacphlip == true, run subworkflow
    if ( params.run_bacphlip ) {
        //
        // MODULE: Predict phage lifestyle using lysogeny specific genes
        //
        ch_bacphlip_lifestyle_tsv = BACPHLIP ( ch_anicluster_reps_fasta ).bacphlip_results
        ch_versions = ch_versions.mix( BACPHLIP.out.versions )

        // TODO: Add ability to automatically add -meta flag to pharokka when multiple sequences are in input fasta
    } else {
        // if run_bacphlip == false, skip subworkflow
        ch_bacphlip_lifestyle_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Phage functional annotation
    ------------------------------------------------------------------------------*/
    // if run_pharokka == true OR run_instrain == true, run subworkflow
    if ( params.run_pharokka || params.run_instrain ) {
        // create channel from params.pharokka_db
        if ( !params.pharokka_db ){
            ch_pharokka_db = null
        } else {
            ch_pharokka_db = Channel.value( file( params.pharokka_db, checkIfExists:true ) )
        }

        //
        // SUBWORKFLOW: Functionally annotate phage sequences
        //
        ch_pharokka_gbk_gz = FASTA_PHAGE_FUNCTION_PHAROKKA ( ch_anicluster_reps_fasta, ch_pharokka_db ).pharokka_gbk_gz
        ch_pharokka_output_tsv = FASTA_PHAGE_FUNCTION_PHAROKKA.out.pharokka_final_output_tsv
        ch_versions = ch_versions.mix( FASTA_PHAGE_FUNCTION_PHAROKKA.out.versions )

        // TODO: Add ability to automatically add -meta flag to pharokka when multiple sequences are in input fasta

        // gunzip proteins for input into instrain
        ch_pharokka_gbk = GUNZIP_VIRUS_PROTEINS ( ch_pharokka_gbk_gz ).gunzip
        ch_versions = ch_versions.mix( GUNZIP_VIRUS_PROTEINS.out.versions )

        // add gene identifier to gbk for inStrain
        ch_pharokka_gbk_mod = ch_pharokka_gbk
        .map { meta, gbk ->
            def gbk_mod = file("${workDir}/${meta.id}_mod.gbk")

            gbk.withReader { source ->
                gbk_mod.withWriter { target ->
                    String line
                    while( line=source.readLine() ) {
                        target << line.replaceAll('/ID=','/gene=') << '\n'
                    }
                }
            }
            return [ meta, gbk_mod]
        }
    } else {
        // if run_pharokka == false AND run_instrain == false, skip subworkflow
        ch_pharokka_gbk_mod = []
        ch_pharokka_output_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Analyze phage microdiversity
    ------------------------------------------------------------------------------*/
    // if run_instrain == true, run subworkflow
    if ( params.run_instrain ) {
        //
        // MODULE: Generate instrain scaffold to bin file
        //
        ch_stb_file_tsv  = INSTRAIN_STB ( ch_anicluster_reps_fasta ).stb
        ch_versions = ch_versions.mix(INSTRAIN_STB.out.versions)

        //
        // SUBWORKFLOW: Assess virus microdiversity within and across samples
        //
        ch_gene_info_tsv = FASTA_MICRODIVERSITY_INSTRAIN ( ch_cluster_rep_alignment_bam, ch_anicluster_reps_fasta, ch_pharokka_gbk_mod, ch_stb_file_tsv ).gene_info_tsv
        ch_versions = ch_versions = ch_versions.mix(FASTA_MICRODIVERSITY_INSTRAIN.out.versions)
    } else {
        // if run_instrain == false, skip subworkflow
        ch_gene_info_tsv = Channel.empty()
    }


    emit:
    virus_enrichment_tsv        = ch_virus_enrichment_tsv
    virus_classification_tsv    = ch_virus_summaries_tsv
    // virus_extension_tsv         = ch_virus_extension_summary_tsv // Don't want to update all other workflow snapshots
    virus_quality_tsv           = ch_quality_summary_tsv
    filtered_viruses_fna_gz     = ch_filtered_viruses_fna_gz
    anicluster_reps_fna_gz      = ch_anicluster_reps_fasta_gz
    // alignment_results_tsv       = ch_alignment_results_tsv   // Inconsistent hash
    host_predictions_tsv        = ch_host_predictions_tsv
    marker_taxonomy_tsv         = ch_marker_taxonomy_tsv
    // bacphlip_lifestyle_tsv      = ch_bacphlip_lifestyle_tsv  // Inconsistent hash
    pharokka_output_tsv         = ch_pharokka_output_tsv
    instrain_gene_info          = ch_gene_info_tsv
    versions                    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
