// Additional BIC processing
include { SAMTOOLS_VARDICT as BIC_SAMTOOLS_VARDICT         } from '../samtools_vardict/main'
include { BIC_POSTPROCESSING                               } from '../bic_postprocess/main'
include { FILTER_MAFS                                      } from '../../../modules/bic/filter_mafs/main'

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

    main:
    versions = Channel.empty()
    params.normalize_vcf_bed = getGenomeAttribute('normalize_vcf_bed')
    params.vep_cache = getGenomeAttribute('vep_cache')
    params.vep_config = getGenomeAttribute('vep_config')

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
        params.vep_cache
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
