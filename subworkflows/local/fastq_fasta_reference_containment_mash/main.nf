//
// Identify sequences contained in readset
//

include { MASH_SKETCH as MASH_SKETCH_ASSEMBLIES } from '../../../modules/nf-core/mash/sketch/main'      // TODO: Update nf-core module to remove optional -r argument
include { MASH_SKETCH as MASH_SKETCH_REFERENCES } from '../../../modules/nf-core/mash/sketch/main'      // TODO: Update nf-core module to remove optional -r argument
include { MASH_PASTE                            } from '../../../modules/local/mash/paste/main'         // TODO: Add module to nf-core
include { MASH_SCREEN                           } from '../../../modules/nf-core/mash/screen/main'      // TODO: Update nf-core module to add meta to sketch

workflow FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH {
    take:
    fastq_gz                // [ [ meta.id ] , [ 1.fastq.gz, 2.fastq.gz ]  , reads (mandatory)
    assembly_fasta_gz       // [ [ meta.id ] , fasta.gz ]                  , assemblies (optional)
    reference_fasta_gz      // [ [ meta.id ] , fasta.gz ]                  , reference sequences (mandatory)
    reference_sketch_msh    // [ [ meta.id ] , fasta.msh ]                 , reference_sequences_sketch (optional)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Create a sketch of assemblies
    //
    ch_assembly_sketch_msh = MASH_SKETCH_ASSEMBLIES ( assembly_fasta_gz ).mash
    ch_versions = ch_versions.mix(MASH_SKETCH_ASSEMBLIES.out.versions.first())

    // if provided, use reference sketch. If not, create one
    if ( reference_sketch_msh ) {
        ch_reference_sketch_msh = reference_sketch_msh
    } else {
        //
        // MODULE: Create sketch of reference sequences
        //
        ch_reference_sketch_msh = MASH_SKETCH_REFERENCES ( reference_fasta_gz ).mash
        ch_versions = ch_versions.mix(MASH_SKETCH_REFERENCES.out.versions.first())
    }

    //
    // MODULE: Combine assembly and reference sketches
    //
    ch_combined_sketch_msh = MASH_PASTE ( ch_assembly_sketch_msh, ch_reference_sketch_msh.collect() ).msh
    ch_versions = ch_versions.mix(MASH_PASTE.out.versions.first())

    // join reads and combined sketch by meta.id
    ch_mash_screen_input = fastq_gz.join( ch_combined_sketch_msh, by:0 )

    //
    // MODULE: Identify contained genomes
    //
    ch_mash_screen_tsv = MASH_SCREEN ( ch_mash_screen_input.map { [ it[0], it[1] ] }, ch_mash_screen_input.map { [ it[0], it[2] ] } ).screen
    ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

    emit:
    mash_screen_results = ch_mash_screen_tsv    // [ [ meta.id ] , fasta.gz ]   , concatenated assemblies and contained references
    versions            = ch_versions           // [ versions.yml ]
}
