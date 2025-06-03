// This should run samtools (if necessary) to generate bam files from Applybsqr cram output (i think)
// Then we will run vardict on the bam files to generate vcf files
// ALso this is paired, so the tumor and normal should be passed in together


include { VARDICTJAVA } from '../../../modules/nf-core/vardictjava/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_NORM } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_TUM } from '../../../modules/nf-core/samtools/convert/main'
include { TABIX_TABIX as TABIX_VARDICT } from '../../../modules/nf-core/tabix/tabix/main'

workflow SAMTOOLS_VARDICT {
    take:
    norm_cram    // channels [meta, cram, crai]
    tumor_cram   // channel: [meta, cram, crai]
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
    CRAM_TO_BAM_NORM (
        norm_cram.map{ _key, meta, cram, crai -> [meta, cram, crai]},
        fasta,
        fasta_fai)

    CRAM_TO_BAM_TUM(
        tumor_cram.map{ _key, meta, cram, crai -> [meta, cram, crai]},
        fasta,
        fasta_fai)

    // This will run vardictjava on the bam files
    normal_input = CRAM_TO_BAM_NORM.out.bam.join(CRAM_TO_BAM_NORM.out.bai).map{ meta, bam, bai -> [meta.patient + meta.sample, meta, bam, bai]}
    tumor_input = CRAM_TO_BAM_TUM.out.bam.join(CRAM_TO_BAM_TUM.out.bai).map{ meta, bam, bai -> [meta.patient, meta, bam, bai]}

    // Using the meta from cram_crai to merge the bams together with the "grouped" meta
    norm_join_prep = cram_crai.map{ meta, _nc, _nci, _tc, _tci -> [ meta.patient + meta.normal_id, meta ] }

    norm_input = norm_join_prep.join(normal_input)

    //reset the join key so tumors can be joined
    tumor_join_prep = norm_input.map{ _oldkey, meta1, _meta2, nb, nbi -> [meta1.patient, meta1, nb, nbi]}
    combined_input = (tumor_join_prep.cross(tumor_input)
        .map{ normal, tumor ->
            def meta = [:]

            meta.id         = "${tumor[1].sample}_vs_${normal[1].normal_id}".toString()
            meta.normal_id  = normal[1].normal_id
            meta.patient    = normal[0]
            meta.sex        = normal[1].sex
            meta.tumor_id   = tumor[1].sample

        [meta, [tumor[2], normal[2]], [tumor[3], normal[3]]]}).combine(intervals)

    VARDICTJAVA(
        combined_input,
        fasta,
        fasta_fai)

    // tabix to index
    TABIX_VARDICT(
        VARDICTJAVA.out.vcf
    )

    // get workflow output ready
    vardict_vcf = (VARDICTJAVA.out.vcf).join(TABIX_VARDICT.out.tbi).map{ meta, vcf, tbi -> [meta + [ variantcaller: 'vardict' ], vcf, tbi] }
    versions = versions.mix(CRAM_TO_BAM_NORM.out.versions)
    versions = versions.mix(CRAM_TO_BAM_TUM.out.versions)
    versions = versions.mix(VARDICTJAVA.out.versions)

    emit:
    vardict_vcf // channel: [meta, vcf, tbi]
    versions // channel: [versions]

}
