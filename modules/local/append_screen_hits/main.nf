process APPEND_SCREEN_HITS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::biopython=1.78 conda-forge::pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' :
        'quay.io/biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' }"

    input:
    path reference_fasta
    tuple val(meta), path(mash_screen)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta_w_mash_hits.fna.gz") , emit: fasta_w_mash_hits
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    append_screen_hits.py \\
        --reference_fasta $reference_fasta \\
        --mash_screen_results $mash_screen \\
        --assembly_fasta $fasta \\
        --prefix $prefix \\
        --output ${prefix}.fasta_w_mash_hits.fna

    gzip ${prefix}.fasta_w_mash_hits.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(echo \$(biopython_version.py 2>&1))
        pandas: \$(echo \$(pandas_version.py 2>&1))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
