//
// De novo identification of viral sequences in assemblies
//

include { MASH_SKETCH as MASH_SKETCH_ASSEMBLIES } from '../../../modules/nf-core/mash/sketch/main'
include { MASH_SKETCH as MASH_SKETCH_REFERENCES } from '../../../modules/nf-core/mash/sketch/main'
include { MASH_PASTE                            } from '../../../modules/local/mash/paste/main'
include { MASH_SCREEN                           } from '../../../modules/nf-core/mash/screen/main'
include { APPEND_SCREEN_HITS                    } from '../../../modules/local/append_screen_hits/main'
include { CAT_CAT as CAT_MASH_SCREEN            } from '../../../modules/nf-core/cat/cat/main'

workflow FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH {
    take:
    fastq_gz    // [ [ meta] , [ 1.fastq.gz, 2.fastq.gz ]   , input reads (mandatory)
    fasta_gz    // [ [ meta] , fasta.gz ]                   , input assemblies (mandatory)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Create a sketch of assemblies
    //
    ch_assembly_sketch_msh = MASH_SKETCH_ASSEMBLIES ( fasta_gz ).mash
    ch_versions = ch_versions.mix(MASH_SKETCH_ASSEMBLIES.out.versions.first())

    // if reference based identification requested, a reference FASTA file must be included
    if ( !params.skip_reference_containment && !params.reference_virus_fasta ) {
        error "[nf-core/phageannotator] ERROR: reference containment requested, but no --reference_virus_fasta provided"
    }

    // if provided, use reference sketch. If not, create one
    if ( params.reference_virus_sketch ) {
        ch_reference_sketch_msh = [ [ id: 'reference' ], file(params.reference_virus_sketch, checkIfExists: true) ]
    } else {
        //
        // MODULE: Create sketch of reference viruses
        //
        ch_reference_sketch_msh = MASH_SKETCH_REFERENCES ( [ [ id: 'reference' ], file(params.reference_virus_fasta, checkIfExists: true) ] ).mash
        ch_versions = ch_versions.mix(MASH_SKETCH_REFERENCES.out.versions.first())
    }

    //
    // MODULE: Combine assembly and reference sketches
    //
    ch_combined_sketch_msh = MASH_PASTE ( ch_assembly_sketch_msh, ch_reference_sketch_msh ).msh
    ch_versions = ch_versions.mix(MASH_PASTE.out.versions.first())

    // join reads and combined sketch within samples
    ch_mash_screen_input = fastq_gz.join(ch_combined_sketch_msh, by:0 )
    //
    // MODULE: Identify contained genomes
    //
    ch_mash_screen_tsv = MASH_SCREEN ( ch_mash_screen_input ).screen
    ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

    // join mash screen and assembly fasta within samples
    ch_append_screen_hits_input = ch_mash_screen_tsv.join( fasta_gz, by:0 )
    // create channel for reference virus fasta
    ch_reference_virus_fasta_gz = [ [ id: 'reference' ], file(params.reference_virus_fasta, checkIfExists: true) ]
    //
    // MODULE: Append screen hits to assemblies
    //
    ch_assembly_w_references_fasta_gz = APPEND_SCREEN_HITS ( ch_append_screen_hits_input, ch_reference_virus_fasta_gz ).assembly_w_screen_hits
    ch_versions = ch_versions.mix(APPEND_SCREEN_HITS.out.versions.first())

    //
    // MODULE: Combine mash screen outputs across samples
    //
    ch_combined_mash_screen_tsv = CAT_MASH_SCREEN( ch_mash_screen_tsv.map{ [ [ id:'all_samples' ], it[1] ] }.groupTuple() ).file_out
    ch_versions = ch_versions.mix(CAT_MASH_SCREEN.out.versions.first())

    emit:
    assembly_w_references_fasta_gz  = ch_assembly_w_references_fasta_gz // [ [ meta] , fasta.gz ] , concatenated assemblies and contained references
    versions                        = ch_versions                       // [ versions.yml ]
}
