process VIROMEQC_INSTALL {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b28a1a551d380ce8d57f9d83894ccb9559b44404:08a4cf815fcef0080ede5b3633202cbda8edf59b-0':
        'biocontainers/mulled-v2-b28a1a551d380ce8d57f9d83894ccb9559b44404:08a4cf815fcef0080ede5b3633202cbda8edf59b-0' }"

    output:
    path("${prefix}/")  , emit: viromeqc_index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "viromeqc_index"
    """
    viromeQC.py \\
        --install \\
        --index_dir ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viromeqc: \$(echo \$(viromeQC.py --version 2>&1 | sed 's/^.*viromeQC.py../ //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "viromeqc_index"
    """
    mkdir -p ${prefix}
    touch ${prefix}/SILVA_132_LSURef_tax_silva.clean.1.bt2
    touch ${prefix}/SILVA_132_LSURef_tax_silva.clean.2.bt2
    touch ${prefix}/SILVA_132_LSURef_tax_silva.clean.3.bt2
    touch ${prefix}/SILVA_132_LSURef_tax_silva.clean.4.bt2
    touch ${prefix}/SILVA_132_LSURef_tax_silva.clean.rev.1.bt2
    touch ${prefix}/SILVA_132_LSURef_tax_silva.clean.rev.2.bt2
    touch ${prefix}/SILVA_132_SSURef_Nr99_tax_silva.clean.1.bt2
    touch ${prefix}/SILVA_132_SSURef_Nr99_tax_silva.clean.2.bt2
    touch ${prefix}/SILVA_132_SSURef_Nr99_tax_silva.clean.3.bt2
    touch ${prefix}/SILVA_132_SSURef_Nr99_tax_silva.clean.4.bt2
    touch ${prefix}/SILVA_132_SSURef_Nr99_tax_silva.clean.rev.1.bt2
    touch ${prefix}/SILVA_132_SSURef_Nr99_tax_silva.clean.rev.2.bt2
    touch ${prefix}/amphora_bacteria.dmnd
    touch ${prefix}/amphora_bacteria_294.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viromeqc: \$(echo \$(viromeQC.py --version 2>&1 | sed 's/^.*viromeQC.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
