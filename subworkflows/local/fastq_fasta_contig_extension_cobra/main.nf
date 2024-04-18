//
// Extend contig ends using COBRA
//

include { BOWTIE2_BUILD                 } from '../../../modules/nf-core/bowtie2/build/main'
include { FASTQ_ALIGN_BOWTIE2           } from '../../nf-core/fastq_align_bowtie2/main'
include { GUNZIP                        } from '../../../modules/nf-core/gunzip/main'
include { COVERM_CONTIG as COVERM_COBRA } from '../../../modules/local/coverm/contig/main'
include { COBRAMETA                     } from '../../../modules/nf-core/cobrameta/main'

workflow FASTQ_FASTA_CONTIG_EXTENSION_COBRA {
    take:
    fastq_gz            // [ [ meta ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] , read files (mandatory)
    fasta_gz            // [ [ meta ], fasta.gz ]                               , assemblies/scaffolds (mandatory)
    viral_contigs_tsv   // [ [ meta ], virus_contigs.tsv ]                      , TSV file containing viral contig names
    cobra_assembler     // val:string {metaspades|megahit|idba}                 , string specifying which assembler was used
    mink                // val:int [ 0-999 ]                                    , integer specifying minimum kmer size used by assembler
    maxk                // val:int [ 0-999 ]                                    , integer specifying maximum kmer size used by assembler

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Make bowtie2 index
    //
    ch_fasta_bt2 = BOWTIE2_BUILD ( fasta_gz ).index.first()
    ch_versions = ch_versions.mix( BOWTIE2_BUILD.out.versions )

    //
    // SUBWORKFLOW: Align reads to bowtie2 index
    //
    ch_fasta_alignment_bam = FASTQ_ALIGN_BOWTIE2 ( fastq_gz, ch_fasta_bt2, false, false, fasta_gz ).bam
    ch_versions = ch_versions.mix( FASTQ_ALIGN_BOWTIE2.out.versions )

    //
    // MODULE: Calculate abundance metrics from BAM file
    //
    ch_coverage_tsv = COVERM_COBRA ( ch_fasta_alignment_bam ).alignment_results
    ch_versions = ch_versions.mix( COVERM_COBRA.out.versions )

    // Remove header from coverage file
    ch_coverage_mod = ch_coverage_tsv
    .map { meta, tsv ->
        def tsv_mod = file("${workDir}/${meta.id}_coverage_mod.tsv")

        tsv.withReader { source ->
            tsv_mod.withWriter { target ->
                String line
                while( line=source.readLine() ) {
                    if ( !line.startsWith('Contig' )) {
                        target << line << '\n'
                    }
                }
            }
        }

        return [ meta, tsv_mod ]
    }

    //
    // MODULE: Gunzip FASTA file
    //
    ch_fasta = GUNZIP ( fasta_gz ).gunzip
    ch_versions = ch_versions.mix( GUNZIP.out.versions )

    // prepare input for cobra
    ch_cobra_input = ch_fasta
        .join(ch_coverage_mod)
        .join(viral_contigs_tsv)
        .join(ch_fasta_alignment_bam)
        .multiMap { it ->
            fasta: [ it[0], it[1] ]
            coverage: [ it[0], it[2] ]
            query: [ it[0], it[3] ]
            bam: [ it[0], it[4] ]
        }


    //
    // MODULE: Extend contigs using COBRA
    //
    ch_extended_fasta_gz = COBRAMETA (
        ch_cobra_input.fasta,
        ch_cobra_input.coverage,
        ch_cobra_input.query,
        ch_cobra_input.bam,
        cobra_assembler,
        mink,
        maxk
    ).all_cobra_assemblies
    ch_versions = ch_versions.mix( COBRAMETA.out.versions )

    emit:
    extended_fasta      = ch_extended_fasta_gz  // [ [ meta ], extended_contigs.fna.gz ]    , FASTA file containing extended contigs
    cobra_summary_tsv   = COBRAMETA.out.joining_summary
    versions            = ch_versions.unique()  // [ versions.yml ]
}
