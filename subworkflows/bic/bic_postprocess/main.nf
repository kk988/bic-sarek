
include { NORMALIZE_TAG_VCF as TAG_VCF_STRELKA } from '../../../modules/bic/normalize_tag_vcf'
include { NORMALIZE_TAG_VCF as TAG_VCF_FREEBAYES } from '../../../modules/bic/normalize_tag_vcf'
include { NORMALIZE_TAG_VCF as TAG_VCF_VARDICT } from '../../../modules/bic/normalize_tag_vcf'
include { NORMALIZE_TAG_VCF as TAG_VCF_MUTECT2 } from '../../../modules/bic/normalize_tag_vcf'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_STRELKA } from '../../../modules/nf-core/bcftools/concat'
include { VCF2MAF } from '../../../modules/nf-core/vcf2maf/main'
include { MERGE_MAFS } from '../../../modules/bic/merge_mafs'

workflow BIC_POSTPROCESSING {
    take:
    strelka_vcf_grouped  //channel [meta, vcf_files, tbi_files]
    mutect2_vcf_tbi // channel [meta, vcf, tbi]
    freebayes_vcf_tbi // channel [meta, vcf, tbi]
    vardict_vcf // channel [meta, vcf, tbi]
    normalize_tag_bed // channel [meta, bed]
    fasta // channel [meta, fasta]
    vep_cache // path vep_cache
    vep_fasta // path vep_fasta

    main:
    versions = Channel.empty()
    // for strelka
        // 1. bcftools concat
        // 2. normalizes+tag
        // 3. fix sample names
        // 4. vcf2maf
    BCFTOOLS_CONCAT_STRELKA(strelka_vcf_grouped)
    strelka_vcf_tbi = Channel.empty().mix(BCFTOOLS_CONCAT_STRELKA.out.vcf.join(BCFTOOLS_CONCAT_STRELKA.out.tbi)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] })
    TAG_VCF_STRELKA(strelka_vcf_tbi, normalize_tag_bed)

    //mutect
    TAG_VCF_MUTECT2(mutect2_vcf_tbi, normalize_tag_bed)

    //freebayes
    TAG_VCF_FREEBAYES(freebayes_vcf_tbi, normalize_tag_bed)

    //vardict
    TAG_VCF_VARDICT(vardict_vcf, normalize_tag_bed)

    vcf_to_maf_input = Channel.empty().mix(TAG_VCF_STRELKA.out.vcf, TAG_VCF_MUTECT2.out.vcf, TAG_VCF_FREEBAYES.out.vcf, TAG_VCF_VARDICT.out.vcf)

    VCF2MAF(vcf_to_maf_input, vep_fasta, vep_cache)

    // then every MAF for a tumor/normal pairing is merged into one large maf.
    // subtract variant caller from the meta data
    MERGE_MAFS(VCF2MAF.out.maf
        .map { meta, maf -> [ meta - meta.subMap('variantcaller') , maf ] }
        .groupTuple())

    merged_maf = MERGE_MAFS.out.merged_maf
    rdas = MERGE_MAFS.out.rdata



    versions = versions.mix(
        BCFTOOLS_CONCAT_STRELKA.out.versions,
        TAG_VCF_STRELKA.out.versions,
        TAG_VCF_MUTECT2.out.versions,
        TAG_VCF_FREEBAYES.out.versions,
        TAG_VCF_VARDICT.out.versions,
        VCF2MAF.out.versions,
        MERGE_MAFS.out.versions,
    )

    emit:
    merged_maf // channel [meta, maf]
    rdas // channel [meta, rda]
    versions // channel [versions]

}
