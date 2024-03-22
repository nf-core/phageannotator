process PROPAGATE_PROPAGATE {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::python=3.8 bioconda::bowtie2 bioconda::samtools bioconda::pysam conda-forge::numpy conda-forge::numba"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1' :
        'biocontainers//vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(provirus_coords)

    output:
    tuple val(meta), path("${prefix}_propagate/${prefix}_propagate.tsv")    , emit: propagate_results
    tuple val(meta), path("${prefix}_propagate/${prefix}_propagate.log")    , emit: log
    path  "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def group = task.ext.prefix ?: "${meta.group}"
    prefix = task.ext.prefix ?: "${meta.group}_${meta.id}"
    """
    if [[ -f ${provirus_coords} ]]
    then
        propagate.py \\
        -f ${fasta} \\
        -v ${provirus_coords} \\
        -r ${reads} \\
        -o ${prefix}_propagate \\
        -t $task.cpus \\
        $args
    else
        mkdir -p ${prefix}_propagate
        echo "No integrated proviruses identified" > ${prefix}_propagate/${prefix}_propagate.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(Propagate --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}
