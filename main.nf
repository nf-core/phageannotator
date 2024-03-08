#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/phageannotator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/phageannotator
    Website: https://nf-co.re/phageannotator
    Slack  : https://nfcore.slack.com/channels/phageannotator
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// PLUGIN: Installed from nextflow plugins
//
include { paramsSummaryMap                      } from 'plugin/nf-validation'

//
// MODULE
//
include { MULTIQC                               } from './modules/nf-core/multiqc/main'
include { FASTQC                                } from './modules/nf-core/fastqc/main'

//
// SUBWORKFLOW
//
include { PIPELINE_INITIALISATION       } from './subworkflows/local/utils_nfcore_phageannotator_pipeline'
include { paramsSummaryMultiqc          } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from './subworkflows/local/utils_nfcore_phageannotator_pipeline'
include { PIPELINE_COMPLETION           } from './subworkflows/local/utils_nfcore_phageannotator_pipeline'

//
// WORKFLOW
//
include { PHAGEANNOTATOR            } from './workflows/phageannotator/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PHAGEANNOTATOR {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_fastq_gz = samplesheet.map { meta, fastq, fasta -> return [ meta, fastq ] }
    ch_fasta_gz = samplesheet.map { meta, fastq, fasta -> return [ meta, fasta ] }

    //
    // WORKFLOW: Run pipeline
    //
    PHAGEANNOTATOR (
        ch_fastq_gz,
        ch_fasta_gz
    )

    emit:
    versions = PHAGEANNOTATOR.out.versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        PIPELINE_INITIALISATION.out.samplesheet.map { meta, fastq, fasta -> return [ meta, fastq ] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PHAGEANNOTATOR (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(NFCORE_PHAGEANNOTATOR.out.versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    ch_multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        ch_multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
