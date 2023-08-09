//
// De novo identification of viral sequences in assemblies
//

include { MASH_SKETCH           } from '../../modules/nf-core/mash/sketch/main'
include { MASH_SCREEN           } from '../../modules/nf-core/mash/screen/main'
include { APPEND_SCREEN_HITS    } from '../../modules/local/append_screen_hits/main'

workflow REFERENCE_BASED_IDENTIFICATION {
    take:
    reads       // [ [ meta] , fastq    ], input reads (mandatory)
    assemblies  // [ [ meta] , fasta    ], input assemblies (mandatory)

    main:
    ch_versions = Channel.empty()

    // use reference sketch if provided, create one if not
    if ( params.reference_id_sketch && !params.skip_reference_based_id ) {
        ch_mash_sketch = [ [ id: 'reference' ], file(params.reference_id_sketch, checkIfExists: true) ]
    } else {
        //
        // MODULE: Create mash sketch of viral genomes
        //
        ch_mash_sketch = MASH_SKETCH ( [ [ id: 'reference' ], file(params.reference_id_fasta, checkIfExists: true) ] ).mash
        ch_versions = ch_versions.mix(MASH_SKETCH.out.versions.first())
    }

    //
    // MODULE: Screen reads for contained genomes
    //
    ch_mash_screen = MASH_SCREEN ( reads, ch_mash_sketch ).screen
    ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

    //
    // MODULE: Append screen hits to input contigs
    //
    ch_assemblies_w_screen_hits = APPEND_SCREEN_HITS ( params.reference_id_fasta, ch_mash_screen, assemblies ).fasta_w_screen_hits
    ch_versions = ch_versions.mix(APPEND_SCREEN_HITS.out.versions.first())


    emit:
    assemblies_w_screen_hits    = ch_assemblies_w_screen_hits   // [ [ meta] , fasta    ]
    versions                    = ch_versions                   // [ versions.yml       ]

}
