process TAG_VCF {
    tag "${meta.id}_${meta.variantcaller}"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*_tagged.vcf') , emit: vcf
    path "versions.yml"            , emit: versions

    script:
    def output_file = "${vcf.baseName}_tagged.vcf"
    """

    cat ${vcf} | egrep "^##" >> ${output_file}
    echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Name of mutation caller">' >> ${output_file}
    cat ${vcf} | egrep "^#CHROM" >> ${output_file}
    cat ${vcf} | egrep -v "^#" | awk -v tag=${meta.variantcaller} 'BEGIN{OFS="\t"}{$8=$8";CALLER="tag;print $0}' >> ${output_file}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | awk '{print \$4}')
    END_VERSIONS
    """
}
