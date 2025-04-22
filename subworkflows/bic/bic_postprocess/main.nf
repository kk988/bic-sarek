
include { POST_STRELKA } from '../../../modules/bic/strelka/post_strelka'
include { TAG_VCF } from '../../../modules/bic/tag_vcf'
include { TAG_VCF as TAG_VCF_FREEBAYES } from '../../../modules/bic/tag_vcf'
include { TAG_VCF as TAG_VCF_VARDICT } from '../../../modules/bic/tag_vcf'

workflow BIC_POSTPROCESSING {
    take:
    strelka_vcf //channel [meta, vcf]
    mutect2_vcf
    freebayes_vcf
    vardict_vcf

    // for strelka
        // 1. Strelka_post
        // 2. tag vcf
        // 3. vcf2maf
    POST_STRELKA(strelka_vcf.groupTuple())
    TAG_VCF(POST_STRELKA.out.vcf)
    VCF_TO_MAF


    // if vcf is mutect2, freebayes, or vardict
        // 1.  this vcf goes through normalize and tag sh
        // note freebayes needs fgrep -v ".:.:." to happen
        // and we remove del events from vardict:  fgrep -v "<DEL>" 
        
        // THen it goes through vcf2maf

    //Every script with a specific tumor/normal sample has vcf2maf.sh script

    // then every MAF for a tumor/normal pairing is merged into one large maf.



}