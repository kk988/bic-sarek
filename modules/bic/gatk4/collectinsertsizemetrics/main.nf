process GATK4_COLLECTINSERTSIZEMETRICS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/juno/bic/depot/singularity/gatk_4.5.0.0_R/gatk_4.5.0.0_R3.6.2.simg':
        '/juno/bic/depot/singularity/gatk_4.5.0.0_R/gatk_4.5.0.0_R3.6.2.simg' }"

    input:
    tuple val(meta), path(cram), path(cram_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.txt"), emit: metrics
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK CollectInsertSizeMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CollectInsertSizeMetrics \\
        -I $cram \\
        -O ${prefix}_is_metrics.txt \\
        -H ${prefix}_is_histogram.pdf \\
        -R $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
