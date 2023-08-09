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
include { APPEND_SCREEN_HITS           } from '../../modules/local/append_screen_hits/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC                        } from '../../modules/nf-core/fastqc/main'
include { GENOMAD_DOWNLOAD              } from '../../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND              } from '../../modules/nf-core/genomad/endtoend/main'
include { MASH_SKETCH                   } from '../../modules/nf-core/mash/sketch/main'
include { MASH_SCREEN                   } from '../../modules/nf-core/mash/screen/main'
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

    Channel
        .fromSamplesheet("input")
        .multiMap { meta, fastq_1, fastq_2, fasta ->
            fastq: [ meta, [ fastq_1, fastq_2 ] ]
            fasta: [ meta, [fasta] ]
        }
        .set { ch_input }

    //
    // MODULE: Analyze reads with FASTQC
    //
    FASTQC ( ch_input.fastq )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //----------------------------------------------------------------------
    //  Reference-based virus identification
    //----------------------------------------------------------------------
    //
    // MODULE: Create mash sketch of viral genomes
    //
    if ( !params.skip_reference_based_id && !params.reference_id_fasta ) {
        error "[nf-core/phageannotator] ERROR: reference-based identification requested, but no --reference_id_fasta provided"
    }
    if ( !params.skip_reference_based_id && params.reference_id_fasta ) {
        ch_mash_sketch = MASH_SKETCH ( [ [ id: 'reference_fasta' ], file(params.reference_id_fasta, checkIfExists: true) ] ).mash
        ch_versions = ch_versions.mix(MASH_SKETCH.out.versions.first())

        //
        // MODULE: Screen reads for contained genomes
        //
        ch_mash_screen = MASH_SCREEN ( ch_input.fastq, ch_mash_sketch ).screen
        ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

        //
        // MODULE: Append screen hits to input contigs
        //
        ch_assembly_w_screen_hits = APPEND_SCREEN_HITS ( params.reference_id_fasta, ch_mash_screen, ch_input.fasta ).fasta_w_screen_hits
        ch_versions = ch_versions.mix(APPEND_SCREEN_HITS.out.versions.first())
    } else {
        ch_assembly_w_screen_hits = ch_input.fasta
    }

    //----------------------------------------------------------------------
    //  De novo virus identification
    //----------------------------------------------------------------------
    //
    // MODULE: Download genomad's database
    //
    ch_genomad_db = GENOMAD_DOWNLOAD ( ).genomad_db
    ch_versions = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions.first())

    //
    // MODULE: Identify and annotate viruses with geNomad
    //
    GENOMAD_ENDTOEND ( ch_assembly_w_screen_hits, ch_genomad_db )
    ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions.first())



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
