process APPENDSCREENHITS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' :
        'biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' }"

    input:
    tuple val(meta), path(mash_screen),  path(assembly_fasta)
    tuple val(meta2), path(reference_fasta)

    output:
    tuple val(meta), path("*.fasta_w_screen_hits.fna.gz")   , emit: assembly_w_screen_hits
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    append_screen_hits.py \\
        --reference_fasta $reference_fasta \\
        --mash_screen_results $mash_screen \\
        --assembly_fasta $assembly_fasta \\
        --prefix $prefix \\
        --output ${prefix}.fasta_w_screen_hits.fna

    gzip ${prefix}.fasta_w_screen_hits.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(echo \$(biopython_version.py 2>&1))
        pandas: \$(echo \$(pandas_version.py 2>&1))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta_w_screen_hits.fna
    gzip ${prefix}.fasta_w_screen_hits.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(echo \$(biopython_version.py 2>&1))
        pandas: \$(echo \$(pandas_version.py 2>&1))
    END_VERSIONS
    """
}
