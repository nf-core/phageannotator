process MASH_PASTE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--hd3113c8_4':
        'biocontainers/mash:2.3--hd3113c8_4' }"

    input:
    tuple val(meta) , path(sketch1)
    tuple val(meta2), path(sketch2)

    output:
    tuple val(meta), path("*.msh"), emit: msh
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${meta2.id}"
    """
    mash \\
        paste \\
        ${prefix}\\
        $sketch1 \\
        $sketch2 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$( mash --version )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${meta2.id}"
    """
    touch ${prefix}.msh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$( mash --version )
    END_VERSIONS
    """
}
