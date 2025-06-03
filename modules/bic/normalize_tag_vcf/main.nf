process NORMALIZE_TAG_VCF {
    tag "${meta.id}_${meta.variantcaller}"
    label 'process_medium'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
    //    'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(target_bed)

    output:
    tuple val(meta), path('*_tagged.vcf') , emit: vcf
    path "versions.yml"                   , emit: versions

    script:
    def int_file1 = "temp.vcf"
    def int_file2 = "temp2.vcf"
    def int_file3 = "temp3.vcf"
    def output_file = "${vcf.baseName}_tagged.vcf"
    def sort = task.ext.sort_cmd ?: ''
    def post_cmd = task.ext.post_cmd ?: "cp ${int_file2} ${int_file3}"

    // if tbi does not end with .tbi, then use tabix
    def tbi_file = tbi.toString()
    def tbi_ext = tbi_file.substring(tbi_file.lastIndexOf('.'))
    def tbi_is_tbi = tbi_ext == '.tbi'
    def create_input_index = tbi_is_tbi ? "" : "tabix ${vcf}"

    """
    ${create_input_index}

    bcftools view -R ${target_bed} ${vcf} \
        ${sort} \
        | bcftools norm -m- \
        > ${int_file1}

    cat ${int_file1} | egrep "^##" >> ${int_file2}
    echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Name of mutation caller">' >> ${int_file2}
    cat ${int_file1} | egrep "^#CHROM" >> ${int_file2}
    cat ${int_file1} | egrep -v "^#" | awk -v tag=${meta.variantcaller} 'BEGIN{OFS="\t"}{\$8=\$8";CALLER="tag;print \$0}' >> ${int_file2}

    ${post_cmd}

    rm ${int_file1} ${int_file2}
    mv ${int_file3} ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | awk '{print \$4}')
    END_VERSIONS
    """
}
