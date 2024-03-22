process PROPAGATE_DEREPCOORDINATES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' :
        'biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' }"

    input:
    tuple val(meta) , path(coordinates)
    tuple val(meta2), path(clusters)

    output:
    tuple val(meta), path("${prefix}_derep_provirus_coords.tsv")  , emit: provirus_coords
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    derep_coordinates.py \\
        --coordinates ${fasta} \\
        --clusters ${clusters} \\
        --output ${prefix}_derep_provirus_coords.tsv \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(echo \$(biopython_version.py 2>&1))
        numpy: \$(echo \$(numpy_version.py 2>&1))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_derep_provirus_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(echo \$(biopython_version.py 2>&1))
        numpy: \$(echo \$(numpy_version.py 2>&1))
    END_VERSIONS
    """
}
