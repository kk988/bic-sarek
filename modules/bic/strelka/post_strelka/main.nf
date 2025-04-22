process POST_STRELKA {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input: 
    tuple val(meta), path(snv_vcf), path(indel_vcf)
    path(target_bed)

    output:
    tuple val(meta), path("*_strelka.vcf") , emit: vcf
    path "versions.yml"                    , emit: versions

    script:
    """ 
    
    bcftoools concat ${snv_vcf} ${indel_vcf} -a | bgzip -c - > intermediate.vcf.gz
    tabix -p vcf intermediate.vcf.gz

    bcftools view -R $TARGET_BED ${TMP}.vcf.gz \
        | bcftools sort - \
        | bcftools norm -m- \
        | perl -pe 's/NORMAL/'${NORMAL}'/ if /^#C/; s/TUMOR/'${TUMOR}'/ if /^#C/' > ${meta.tumor_id}___${meta.normal_id}_strelka.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS

    """

}