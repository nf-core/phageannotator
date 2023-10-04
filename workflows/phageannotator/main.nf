/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { SEQKIT_SEQ                                } from '../../modules/local/seqkit/seq/main'                                    // TODO: Add to nf-core
include { AWK as AWK_GENOMAD                        } from '../../modules/local/awk/main'                                           // TODO: Add to nf-core
include { APPENDSCREENHITS                          } from '../../modules/local/appendscreenhits/main'
include { AWK as AWK_CHECKV                         } from '../../modules/local/awk/main'                                           // TODO: Add to nf-core
include { QUALITYFILTERVIRUSES                      } from '../../modules/local/qualityfilterviruses/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH    } from '../../subworkflows/local/fastq_fasta_reference_containment_mash/main'   // TODO: Add to nf-core
include { FASTA_VIRUS_CLASSIFICATION_GENOMAD        } from '../../subworkflows/local/fasta_virus_classification_genomad/main'       // TODO: Add to nf-core
include { FASTA_VIRUS_QUALITY_CHECKV                } from '../../subworkflows/local/fasta_virus_quality_checkv/main'               // TODO: Add to nf-core

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_CAT as CAT_MASHSCREEN     } from '../../modules/nf-core/cat/cat/main'

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
    fasta_gz    // [ [ meta ], assembly.fastq.gz ]  , assemblies/genomes (mandatory)

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
        ch_reference_virus_fasta_gz = [ [ id:'reference_viruses' ], file( params.reference_virus_fasta, checkIfExists:true ) ]

        // create channel from params.reference_virus_sketch
        if ( !params.reference_virus_sketch ){
            ch_reference_virus_sketch_msh = null
        } else {
            ch_reference_virus_sketch_msh = [ [ id:'reference_viruses' ], file( params.reference_virus_sketch, checkIfExists:true ) ]
        }

        //
        // SUBWORKFLOW: Identify contained reference genomes
        //
        ch_containment_results_tsv = FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH ( fastq_gz, ch_filtered_input_fasta_gz, ch_reference_virus_fasta_gz, ch_reference_virus_sketch_msh ).mash_screen_results
        ch_versions = ch_versions.mix(FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH.out.versions.first())

        // join mash screen and assembly fasta by meta.id
        ch_append_screen_hits_input = ch_containment_results_tsv.join( ch_filtered_input_fasta_gz, by:0 )

        //
        // MODULE: Append screen hits to assemblies
        //
        ch_assembly_w_references_fasta_gz = APPENDSCREENHITS ( ch_append_screen_hits_input, ch_reference_virus_fasta_gz ).assembly_w_screen_hits
        ch_versions = ch_versions.mix(APPENDSCREENHITS.out.versions.first())

        //
        // MODULE: Combine mash screen outputs across samples
        //
        ch_combined_mash_screen_tsv = CAT_MASHSCREEN( ch_containment_results_tsv.map{ [ [ id:'all_samples' ], it[1] ] }.groupTuple() ).file_out
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
    ch_versions = ch_versions.mix(FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.versions.first())

    // create channel for genomad virus summary files
    ch_virus_summaries_tsv = FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.virus_summaries_tsv

    // create a channel for combining geNomad virus summaries
    ch_awk_genomad_input = ch_virus_summaries_tsv.map { [ [ id:'all_samples' ], it[1] ] }.groupTuple()

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
    FASTA_VIRUS_QUALITY_CHECKV ( ch_viruses_fna_gz, ch_genomad_db )
    ch_versions = ch_versions.mix(FASTA_VIRUS_QUALITY_CHECKV.out.versions.first())

    // create a channel for quality summaries
    ch_quality_summaries_tsv = FASTA_VIRUS_QUALITY_CHECKV.out.quality_summary_tsv.map { [ [ id:'all_samples' ], it[1] ] }.groupTuple()

    //
    // MODULE: Combine quality summaries across samples
    //
    ch_combined_quality_summaries_tsv = AWK_CHECKV ( ch_quality_summaries_tsv ).file_out
    ch_versions = ch_versions.mix(AWK_CHECKV.out.versions.first())

    // create channel for input into QUALITY_FILTER_VIRUSES
    ch_quality_filter_viruses_input1 = FASTA_VIRUS_QUALITY_CHECKV.out.viruses_fna_gz.join(FASTA_VIRUS_QUALITY_CHECKV.out.proviruses_fna_gz)
    ch_quality_filter_viruses_input2 = ch_quality_filter_viruses_input1.join(FASTA_VIRUS_QUALITY_CHECKV.out.quality_summary_tsv)

    //
    // MODULE: Quality filter viruses
    //
    ch_filtered_viruses_fna_gz = QUALITYFILTERVIRUSES ( ch_quality_filter_viruses_input2 ).filtered_viruses
    ch_versions = ch_versions.mix(QUALITYFILTERVIRUSES.out.versions.first())


    emit:
    filtered_viruses_fna_gz     = ch_filtered_viruses_fna_gz
    reference_containment_tsv   = ch_combined_mash_screen_tsv
    virus_classification_tsv    = ch_combined_virus_summaries_tsv
    virus_quality_tsv           = ch_combined_quality_summaries_tsv
    versions                    = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
