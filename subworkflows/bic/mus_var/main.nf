// Additional BIC processing
include { SAMTOOLS_VARDICT as BIC_SAMTOOLS_VARDICT } from '../samtools_vardict/main'
include { BIC_POSTPROCESSING                       } from '../bic_postprocess/main'
include { FILTER_MAFS                              } from '../../../modules/bic/filter_mafs/main'
include { GATK4_COLLECTALIGNMENTSUMMARYMETRICS     } from '../../../modules/bic/gatk4/collectalignmentsummarymetrics/main'
include { GATK4_COLLECTINSERTSIZEMETRICS           } from '../../../modules/bic/gatk4/collectinsertsizemetrics/main'
include { GATK4_COLLECTHSMETRICS                   } from '../../../modules/bic/gatk4/collecthsmetrics/main'
include { BIC_SAMPLE_QC                            } from '../../../modules/bic/bic_sample_qc/main'

workflow MUS_VAR{
    take:
    cram_variant_calling_normal_to_cross
    cram_variant_calling_pair_to_cross
    cram_variant_calling_pair
    fasta
    fasta_fai
    intervals_bed_combined
    vcf_strelka
    vcf_mutect2
    vcf_freebayes
    bic_qc_reports

    main:
    versions = Channel.empty()
    params.normalize_vcf_bed = getGenomeAttribute('normalize_vcf_bed')
    params.vep_cache = getGenomeAttribute('vep_cache')
    params.vep_config = getGenomeAttribute('vep_config')
    params.vep_fasta = getGenomeAttribute('vep_fasta')
    params.fasta_bed = getGenomeAttribute('fasta_bed')
    def targets_ilist = getGenomeAttribute('targets_ilist')
    def bait_ilist = getGenomeAttribute('bait_ilist')

    if (!targets_ilist || !file(targets_ilist).exists() || !bait_ilist || !file(bait_ilist).exists()) {
        log.warn "No targets_ilist or bait_ilist found for genome ${params.genome} in the genomes config file, skipping hsmetrics"
    }
    else {
        GATK4_COLLECTHSMETRICS(
            (cram_variant_calling_normal_to_cross
                .mix(cram_variant_calling_pair_to_cross))
                .map{ _meta_pt, meta, cram, crai -> [ meta, cram, crai ] },
            bait_ilist,
            targets_ilist,
            fasta,
            fasta_fai
        )
        versions = versions.mix(GATK4_COLLECTHSMETRICS.out.versions)
        bic_qc_reports = bic_qc_reports.mix(GATK4_COLLECTHSMETRICS.out.metrics.collect{ _meta, report  -> [ report ] })
    }

    // Alignment QC
    //
    GATK4_COLLECTALIGNMENTSUMMARYMETRICS(
        (cram_variant_calling_normal_to_cross
            .mix(cram_variant_calling_pair_to_cross))
            .map{ _meta_pt, meta, cram, crai -> [ meta, cram, crai ] },
        fasta,
        fasta_fai
    )
    versions = versions.mix(GATK4_COLLECTALIGNMENTSUMMARYMETRICS.out.versions)
    bic_qc_reports = bic_qc_reports.mix(GATK4_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics.collect{ _meta, report  -> [ report ] })

    GATK4_COLLECTINSERTSIZEMETRICS(
        (cram_variant_calling_normal_to_cross
            .mix(cram_variant_calling_pair_to_cross))
            .map{ _meta_pt, meta, cram, crai -> [ meta, cram, crai ] },
        fasta,
        fasta_fai
    )
    versions = versions.mix(GATK4_COLLECTINSERTSIZEMETRICS.out.versions)
    bic_qc_reports = bic_qc_reports.mix(GATK4_COLLECTINSERTSIZEMETRICS.out.metrics.collect{ _meta, report  -> [ report ] })

    BIC_SAMPLE_QC(
        bic_qc_reports.collect(),
        params.input,
        params.qc_control_csv
    )

    // BIC variant calling
    //
    BIC_SAMTOOLS_VARDICT(
        cram_variant_calling_normal_to_cross,
        cram_variant_calling_pair_to_cross,
        cram_variant_calling_pair,
        fasta,
        fasta_fai,
        intervals_bed_combined
    )

    versions = versions.mix(BIC_SAMTOOLS_VARDICT.out.versions)

    // end BIC variant calling
    // BIC POST PROCESSING

    normalize_tag_bed = params.normalize_vcf_bed ? Channel.fromPath(params.normalize_vcf_bed).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

    strelka_vcf_grouped = vcf_strelka
        .groupTuple() // Group VCF files by meta
        .map { meta, vcf_files ->
            // Construct the corresponding .tbi paths for each VCF file
            def tbi_files = vcf_files.collect { vcf_file ->
                def tbi_file = file(vcf_file.toString() + ".tbi") // Construct the .tbi path
                tbi_file.exists() ? tbi_file : null
                tbi_file
            }
            [meta, vcf_files, tbi_files] // Return the desired structure
        }

    mutect2_vcf_tbi = vcf_mutect2.map { meta, vcf_file ->
            def tbi_file = file(vcf_file.toString() + ".tbi") // Construct the path to the index file
            if (!tbi_file.exists()) {
                log.info "Index file not found for VCF: ${vcf_file} TBI: ${tbi_file}"
                tbi_file = params.fasta
            }
            [meta, vcf_file, tbi_file] // Add the index file to the channel
        }

    freebayes_vcf_tbi = vcf_freebayes.map { meta, vcf_file ->
            def tbi_file = file(vcf_file.toString() + ".tbi") // Construct the path to the index file
            if (!tbi_file.exists()) {
                log.info "Index file not found for VCF: ${vcf_file} TBI: ${tbi_file}"
                tbi_file = params.fasta
            }
            [meta, vcf_file, tbi_file] // Add the index file to the channel
        }

    vardict_vcf = BIC_SAMTOOLS_VARDICT.out.vardict_vcf

    BIC_POSTPROCESSING(strelka_vcf_grouped,
        mutect2_vcf_tbi,
        freebayes_vcf_tbi,
        vardict_vcf,
        normalize_tag_bed,
        fasta,
        params.vep_cache,
        params.vep_fasta
    )

    rdas = BIC_POSTPROCESSING.out.rdas.map{ _meta, rda -> [ rda ] }.collect()
    FILTER_MAFS(rdas)

    versions = versions.mix(BIC_POSTPROCESSING.out.versions, FILTER_MAFS.out.versions)

    // end of BIC Post Processing

    emit:
    versions // channel [meta, versions]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
// Copied from `main.nf`
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}
