process PROPAGATE_PROPAGATE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1' :
        'biocontainers/vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(provirus_coords)

    output:
    tuple val(meta), path("${prefix}_propagate/${prefix}_propagate.tsv")    , emit: results
    tuple val(meta), path("${prefix}_propagate/${prefix}_propagate.log")    , emit: log     , optional: true
    path  "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    propagate.py \\
        -f ${fasta} \\
        -v ${provirus_coords} \\
        -r ${reads} \\
        -o ${prefix}_propagate \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        propagate: \$(propagate.py --version 2>&1 | sed 's/^.*PropagAtE v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_propagate
    touch ${prefix}_propagate/${prefix}_propagate.tsv
    touch ${prefix}_propagate/${prefix}_propagate.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        propagate: \$(propagate.py --version 2>&1 | sed 's/^.*PropagAtE v//; s/ .*\$//')
    END_VERSIONS
    """
}
