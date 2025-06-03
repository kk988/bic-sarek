process MERGE_MAFS {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(mafs)

    output:
    tuple val(meta), path('*merge.maf.tsv.gz'), emit: merged_maf
    tuple val(meta), path('*.rda')            , emit: rdata
    path "versions.yml"                       , emit: versions

    script:
    template 'merge_mafs.R'

}
