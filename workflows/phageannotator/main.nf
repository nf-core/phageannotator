/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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
include { APPEND_SCREEN_HITS                        } from '../../modules/local/append_screen_hits/main'
include { AWK as AWK_CHECKV                         } from '../../modules/local/awk/main'                                           // TODO: Add to nf-core
include { QUALITY_FILTER_VIRUSES                    } from '../../modules/local/quality_filter_viruses/main'

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
include { FASTQC                        } from '../../modules/nf-core/fastqc/main'
include { CAT_CAT as CAT_MASH_SCREEN    } from '../../modules/nf-core/cat/cat/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                       } from '../../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PHAGEANNOTATOR {

    ch_versions = Channel.empty()


    /*----------------------------------------------------------------------------
        Read in samplesheet
    ------------------------------------------------------------------------------*/
    // read inputs using nf-validation plugin
    Channel
        .fromSamplesheet("input")
        .multiMap { meta, fastq_1, fastq_2, fasta ->
            fastq_gz: [ meta, [ fastq_1, fastq_2 ] ]
            fasta_gz: [ meta, [ fasta ] ]
        }
        .set { ch_input }


    /*----------------------------------------------------------------------------
        Analyze input reads
    ------------------------------------------------------------------------------*/
    //
    // MODULE: Analyze reads
    //
    FASTQC ( ch_input.fastq_gz )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    /*----------------------------------------------------------------------------
        Filter input assemblies
    ------------------------------------------------------------------------------*/
    //
    // MODULE: Filter assemblies by length
    //
    ch_filtered_input_fasta_gz = SEQKIT_SEQ ( ch_input.fasta_gz ).fastx
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
        ch_containment_results_tsv = FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH ( ch_input.fastq_gz, ch_filtered_input_fasta_gz, ch_reference_virus_fasta_gz, ch_reference_virus_sketch_msh ).mash_screen_results
        ch_versions = ch_versions.mix(FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH.out.versions.first())

        // join mash screen and assembly fasta by meta.id
        ch_append_screen_hits_input = ch_containment_results_tsv.join( ch_filtered_input_fasta_gz, by:0 )

        //
        // MODULE: Append screen hits to assemblies
        //
        ch_assembly_w_references_fasta_gz = APPEND_SCREEN_HITS ( ch_append_screen_hits_input, ch_reference_virus_fasta_gz ).assembly_w_screen_hits
        ch_versions = ch_versions.mix(APPEND_SCREEN_HITS.out.versions.first())

        //
        // MODULE: Combine mash screen outputs across samples
        //
        ch_combined_mash_screen_tsv = CAT_MASH_SCREEN( ch_containment_results_tsv.map{ [ [ id:'all_samples' ], it[1] ] }.groupTuple() ).file_out
        ch_versions = ch_versions.mix(CAT_MASH_SCREEN.out.versions.first())
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
    ch_viruses_fasta_gz = FASTA_VIRUS_CLASSIFICATION_GENOMAD ( ch_assembly_w_references_fasta_gz, ch_genomad_db ).viruses_fasta_gz
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
    FASTA_VIRUS_QUALITY_CHECKV ( ch_viruses_fasta_gz, ch_genomad_db )
    ch_versions = ch_versions.mix(FASTA_VIRUS_QUALITY_CHECKV.out.versions.first())

    // create a channel for quality summaries
    ch_quality_summary_tsv = FASTA_VIRUS_QUALITY_CHECKV.out.quality_summary_tsv

    //
    // MODULE: Combine quality summaries across samples
    //
    ch_combined_quality_summaries_tsv = AWK_CHECKV ( ch_quality_summary_tsv ).file_out
    ch_versions = ch_versions.mix(AWK_CHECKV.out.versions.first())

    // create channel for input into QUALITY_FILTER_VIRUSES
    ch_quality_filter_viruses_input = FASTA_VIRUS_QUALITY_CHECKV.out.viruses.join(FASTA_VIRUS_QUALITY_CHECKV.out.proviruses).join(ch_quality_summary_tsv)

    //
    // MODULE: Quality filter viruses
    //
    ch_filtered_viruses_fasta_gz = QUALITY_FILTER_VIRUSES ( ch_quality_filter_viruses_input ).filtered_viruses
    ch_versions = ch_versions.mix(QUALITY_FILTER_VIRUSES.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPhageannotator.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPhageannotator.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
