
include { POST_STRELKA } from '../../../modules/bic/strelka/post_strelka'
include { NORMALIZE_TAG_VCF as TAG_VCF_STRELKA } from '../../../modules/bic/normalize_tag_vcf'
include { NORMALIZE_TAG_VCF as TAG_VCF_FREEBAYES } from '../../../modules/bic/normalize_tag_vcf'
include { NORMALIZE_TAG_VCF as TAG_VCF_VARDICT } from '../../../modules/bic/normalize_tag_vcf'
include { NORMALIZE_TAG_VCF as TAG_VCF_MUTECT2 } from '../../../modules/bic/normalize_tag_vcf'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_STRELKA } from '../../../modules/nf-core/bcftools/concat'
include { VCF2MAF } from '../../../modules/nf-core/vcf2maf/main'

workflow BIC_POSTPROCESSING {
    take:
    strelka_vcf //channel [meta, vcf]
    mutect2_vcf // channel [meta, vcf, tbi]
    freebayes_vcf // channel [meta, vcf, tbi]
    vardict_vcf // channel [meta, vcf, tbi]
    noramlize_vcf_bed // channel [meta, bed]
    fasta // channel [meta, fasta]
    vep_cache // path vep_cache

    main:
    // for strelka
        // 1. bcftools concat
        // 2. normalizes+tag
        // 3. fix sample names
        // 4. vcf2maf
    BCFTOOLS_CONCAT_STRELKA(strelka_vcf.groupTuple())
    strelka_vcf_tbi = Channel.empty().mix(BCFTOOLS_CONCAT_STRELKA.out.vcf.zip(BCFTOOLS_CONCAT_STRELKA.out.tbi)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] })
    TAG_VCF_STRELKA(strelka_vcf_tbi, noramlize_vcf_bed)

    //mutect
    TAG_VCF_MUTECT2(mutect2_vcf, noramlize_vcf_bed)

    //freebayes
    TAG_VCF_FREEBAYES(freebayes_vcf, noramlize_vcf_bed)

    //vardict
    TAG_VCF_VARDICT(vardict_vcf, noramlize_vcf_bed)

    vcf_to_maf_input = Channel.empty().mix(TAG_VCF_STRELKA.out.vcf, TAG_VCF_MUTECT2.out.vcf, TAG_VCF_FREEBAYES.out.vcf, TAG_VCF_VARDICT.out.vcf)
    VCF2MAF(vcf_to_maf_input, fasta.first().map{ _meta, file -> file}, Channel.value(vep_cache))

    // then every MAF for a tumor/normal pairing is merged into one large maf.



}
