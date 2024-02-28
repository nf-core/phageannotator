process COBRA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cobra-meta:1.2.3--pyhdfd78af_0':
        'biocontainers/cobra-meta:1.2.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(coverage), path(viral_contigs), path(bam)
    val assembler
    val mink
    val maxk

    output:
    tuple val(meta), path("${prefix}/COBRA_extended.fasta.gz")  , emit: extended_assemblies
    tuple val(meta), path("${prefix}/COBRA_joining_summary.txt"), emit: joining_summary
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cobra-meta \\
        --fasta ${fasta} \\
        --coverage ${coverage} \\
        --query ${viral_contigs} \\
        --mapping ${bam} \\
        --assembler ${assembler} \\
        --mink ${mink} \\
        --maxk ${maxk} \\
        --threads ${task.cpus} \\
        --output ${prefix} \\
        $args

    # if the output is empty, use input file
    if [ -e ${prefix}/COBRA_extended.fasta ]
    then
        gzip ${prefix}/*.fasta
        cat ${prefix}/*fasta.gz > ${prefix}/COBRA_extended.fasta.gz
    else
        gzip -c ${fasta} > ${prefix}/COBRA_extended.fasta
        touch ${prefix}/COBRA_joining_summary.txt
    fi

    gzip ${prefix}/*.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra v//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/COBRA_extended.fasta
    gzip ${prefix}/COBRA_extended.fasta
    touch ${prefix}/COBRA_joining_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra v//' ))
    END_VERSIONS
    """
}
