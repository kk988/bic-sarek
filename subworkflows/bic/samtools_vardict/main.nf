// This should run samtools (if necessary) to generate bam files from Applybsqr cram output (i think)
// Then we will run vardict on the bam files to generate vcf files
// ALso this is paired, so the tumor and normal should be passed in together


include { VARDICTJAVA } from '../../../modules/nf-core/vardictjava/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_NORM } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_TUM } from '../../../modules/nf-core/samtools/convert/main'

workflow SAMTOOSL_VARDICT {
    take:
    cram_crai    // channel: [meta, normal_cram, normal_crai tumor_cram, tumor_crai]
    fasta        // channel: [meta, fasta]
    fasta_fai    // channel: [meta, fasta_fai]
    intervals    // channel: [meta, bed]


    main:
    versions = Channel.empty()

    // tuple val(meta), path(bams) (Should be a List (tumor/normal)), path(bais) (Should be a List (tumor/normal)), path(bed) ( vardict_bed )
    //tuple val(meta2), path(fasta) (from bic_musvar.config)
    //tuple val(meta3), path(fasta_fai) (from bic_musvar.config)

    // CRAM_TO_BAM(cram_crai, fasta, fasta_fai) // This will convert the cram to bam
    // First set up channel for cram to bam
    CRAM_TO_BAM_NORM(
        cram_crai.map{ meta, normal_cram, normal_crai -> [meta, normal_cram, normal_crai] },
        fasta,
        fasta_fai)

    CRAM_TO_BAM_TUM(
        cram_crai.map{ meta, tumor_cram, tumor_crai -> [meta, tumor_cram, tumor_crai] },
        fasta,
        fasta_fai)

    // This will run vardictjava on the bam files
    normal_input = CRAM_TO_BAM_NORM.out.bam.join(CRAM_TO_BAM_NORM.out.bai)
    tumor_input = CRAM_TO_BAM_TUM.out.bam.join(CRAM_TO_BAM_TUM.out.bai)
    combined_input = (normal_input.join(tumor_input)).combine(intervals)
    vardict_input = combined_input.map{ meta, nb, ni, tb, ti, meta_i, bed ->
        [meta, [nb, tb], [ni, ti], bed]
    }
    VARDICTJAVA(
        vardict_input,
        fasta,
        fasta_fai)

    // get workflow output ready
    vardict_vcf = VARDICTJAVA.out.vcf
    versions = versions.mix(CRAM_TO_BAM_NORM.out.versions)
    versions = versions.mix(CRAM_TO_BAM_TUM.out.versions)
    versions = versions.mix(VARDICTJAVA.out.versions)

    emit:
    vardict_vcf // channel: [meta, vcf]
    versions // channel: [versions]

}
