process PRODIGAL_PRODIGALGV {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prodigal-gv:2.11.0--he4a0461_2':
        'biocontainers/prodigal-gv:2.11.0--he4a0461_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.prodigalgv.faa.gz")    , emit: faa
    tuple val(meta), path("*.prodigalgv.fna.gz")    , emit: fna
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prodigal-gv \\
        -i $fasta \\
        -d ${prefix}.prodigalgv.fna \\
        -a ${prefix}.prodigalgv.faa \\
        $args

    gzip ${prefix}.prodigalgv.fna
    gzip ${prefix}.prodigalgv.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigalgv: \$(echo \$(prodigal-gv --version 2>&1) | awk 'NF{ print \$NF }')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.prodigalgv.fna
    touch ${prefix}.prodigalgv.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigalgv: \$(echo \$(prodigal-gv --version 2>&1) | awk 'NF{ print \$NF }')
    END_VERSIONS
    """
}
