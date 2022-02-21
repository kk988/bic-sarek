//
// GERMLINE VARIANT CALLING
//

include { BGZIP as BGZIP_DEEPVARIANT_GVCF             } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_DEEPVARIANT_VCF              } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_FREEBAYES                    } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_HAPLOTYPECALLER              } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_DIPLOID                } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SMALL_INDELS           } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SV                     } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_STRELKA                      } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_STRELKA_GENOME               } from '../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_GVCF_DEEPVARIANT       } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_DEEPVARIANT        } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_FREEBAYES          } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_HAPLOTYPECALLER    } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_DIPLOID      } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SMALL_INDELS } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SV           } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA            } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA_GENOME     } from '../../modules/local/concat_vcf/main'
include { DEEPVARIANT                                 } from '../../modules/nf-core/modules/deepvariant/main'
include { FREEBAYES                                   } from '../../modules/nf-core/modules/freebayes/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER    } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING         } from '../../subworkflows/nf-core/joint_germline_variant_calling/main'
include { MANTA_GERMLINE                              } from '../../modules/nf-core/modules/manta/germline/main'
include { STRELKA_GERMLINE                            } from '../../modules/nf-core/modules/strelka/germline/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIP_TIDDIT_SV   } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_DEEPVARIANT_GVCF       } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DEEPVARIANT_VCF        } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_FREEBAYES              } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_HAPLOTYPECALLER        } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_MANTA                  } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_STRELKA                } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TIDDIT_SV                                   } from '../../modules/nf-core/modules/tiddit/sv/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        cram_recalibrated            // channel: [mandatory] cram
        dbsnp                        // channel: [mandatory] dbsnp
        dbsnp_tbi                    // channel: [mandatory] dbsnp_tbi
        dict                         // channel: [mandatory] dict
        fasta                        // channel: [mandatory] fasta
        fasta_fai                    // channel: [mandatory] fasta_fai
        intervals                    // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi         // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz_tbi // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        intervals_bed_combine_gz     // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        num_intervals                // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        // joint_germline               // val: true/false on whether to run joint_germline calling, only works in combination with haplotypecaller at the moment

    main:

    ch_versions = Channel.empty()

    // Remap channel with intervals
    cram_recalibrated_intervals = cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            sample = meta.sample
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? sample + "_" + new_intervals.baseName : sample
            [[ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ], cram, crai, new_intervals]
        }

    // Remap channel with gziped intervals + indexes
    cram_recalibrated_intervals_gz_tbi = cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            sample = meta.sample
            new_bed = bed.simpleName != "no_intervals" ? bed : []
            new_tbi = tbi.simpleName != "no_intervals" ? tbi : []
            id = new_bed ? sample + "_" + new_bed.simpleName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, new_bed, new_tbi]
        }

    // DEEPVARIANT

    //TODO: benchmark if it is better to provide multiple bed files & run on multiple machines + mergeing afterwards || one containing all intervals and run on one larger machine
    // Deepvariant: https://github.com/google/deepvariant/issues/510

    DEEPVARIANT(
        cram_recalibrated_intervals,
        fasta,
        fasta_fai)

    // Only when no intervals
    TABIX_DEEPVARIANT_VCF(DEEPVARIANT.out.vcf)
    TABIX_DEEPVARIANT_GVCF(DEEPVARIANT.out.gvcf)

    // Only when using intervals
    BGZIP_DEEPVARIANT_VCF(DEEPVARIANT.out.vcf)
    BGZIP_DEEPVARIANT_GVCF(DEEPVARIANT.out.gvcf)

    CONCAT_VCF_DEEPVARIANT(
        BGZIP_DEEPVARIANT_VCF.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    CONCAT_GVCF_DEEPVARIANT(
        BGZIP_DEEPVARIANT_GVCF.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    deepvariant_vcf = channel.empty().mix(
        CONCAT_GVCF_DEEPVARIANT.out.vcf,
        CONCAT_VCF_DEEPVARIANT.out.vcf,
        DEEPVARIANT.out.gvcf,
        DEEPVARIANT.out.vcf)

    // FREEBAYES

    // Remap channel for Freebayes
    cram_recalibrated_intervals_freebayes = cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            sample = meta.sample
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? sample + "_" + new_intervals.baseName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, [], [], new_intervals]
        }

    FREEBAYES(
        cram_recalibrated_intervals_freebayes,
        fasta,
        fasta_fai,
        [], [], [])

    // Only when no intervals
    TABIX_FREEBAYES(FREEBAYES.out.vcf)

    // Only when using intervals
    BGZIP_FREEBAYES(FREEBAYES.out.vcf)

    CONCAT_VCF_FREEBAYES(
        BGZIP_FREEBAYES.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    freebayes_vcf = channel.empty().mix(
        CONCAT_VCF_FREEBAYES.out.vcf,
        FREEBAYES.out.vcf)

    // HAPLOTYPECALLER

    HAPLOTYPECALLER(
        cram_recalibrated_intervals,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi)

    // Only when no intervals
    TABIX_HAPLOTYPECALLER(HAPLOTYPECALLER.out.vcf)

    // Only when using intervals
    BGZIP_HAPLOTYPECALLER(HAPLOTYPECALLER.out.vcf)

    CONCAT_VCF_HAPLOTYPECALLER(
        BGZIP_HAPLOTYPECALLER.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    haplotypecaller_gvcf = Channel.empty().mix(
        CONCAT_VCF_HAPLOTYPECALLER.out.vcf,
        HAPLOTYPECALLER.out.vcf)

    // if (joint_germline) {
    //     run_haplotypecaller = false
    //     run_vqsr            = true //parameter?
    //     some feedback from gavin
    //     GATK_JOINT_GERMLINE_VARIANT_CALLING(
    //         haplotypecaller_vcf_gz_tbi,
    //         run_haplotypecaller,
    //         run_vqsr,
    //         fasta,
    //         fasta_fai,
    //         dict,
    //         dbsnp,
    //         dbsnp_tbi,
    //         "joined",
    //          allelespecific?
    //          resources?
    //          annotation?
    //         "BOTH",
    //         true,
    //         truthsensitivity -> parameter or module?
    //     )
    //     ch_versions = ch_versions.mix(GATK_JOINT_GERMLINE_VARIANT_CALLING.out.versions)
    // }

    // MANTA
    // TODO: Research if splitting by intervals is ok, we pretend for now it is fine.
    // Seems to be the consensus on upstream modules implementation too

    MANTA_GERMLINE(
        cram_recalibrated_intervals_gz_tbi,
        fasta,
        fasta_fai)

    // Figure out if using intervals or no_intervals
    MANTA_GERMLINE.out.candidate_small_indels_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_small_indels_vcf}

    MANTA_GERMLINE.out.candidate_sv_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_sv_vcf}

    MANTA_GERMLINE.out.diploid_sv_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_diploid_sv_vcf}

    // Only when using intervals
    BGZIP_MANTA_DIPLOID(MANTA_GERMLINE.out.diploid_sv_vcf)

    CONCAT_VCF_MANTA_DIPLOID(
        BGZIP_MANTA_DIPLOID.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_MANTA_SMALL_INDELS(MANTA_GERMLINE.out.candidate_small_indels_vcf)

    CONCAT_VCF_MANTA_SMALL_INDELS(
        BGZIP_MANTA_SMALL_INDELS.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_MANTA_SV(MANTA_GERMLINE.out.candidate_sv_vcf)

    CONCAT_VCF_MANTA_SV(
        BGZIP_MANTA_SV.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    manta_vcf = Channel.empty().mix(
        CONCAT_VCF_MANTA_DIPLOID.out.vcf,
        CONCAT_VCF_MANTA_SMALL_INDELS.out.vcf,
        CONCAT_VCF_MANTA_SV.out.vcf,
        manta_diploid_sv_vcf.no_intervals,
        manta_small_indels_vcf.no_intervals,
        manta_sv_vcf.no_intervals)

    // STRELKA
    // TODO: Research if splitting by intervals is ok, we pretend for now it is fine.
    // Seems to be the consensus on upstream modules implementation too

    STRELKA_GERMLINE(
        cram_recalibrated_intervals_gz_tbi,
        fasta,
        fasta_fai)

    // Figure out if using intervals or no_intervals
    STRELKA_GERMLINE.out.vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{strelka_vcf}

    STRELKA_GERMLINE.out.genome_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{strelka_genome_vcf}

    // Only when using intervals
    BGZIP_STRELKA(STRELKA_GERMLINE.out.vcf)

    CONCAT_VCF_STRELKA(
        BGZIP_STRELKA.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_STRELKA_GENOME(STRELKA_GERMLINE.out.genome_vcf)

    CONCAT_VCF_STRELKA_GENOME(
        BGZIP_STRELKA_GENOME.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    strelka_vcf = Channel.empty().mix(
        CONCAT_VCF_STRELKA.out.vcf,
        CONCAT_VCF_STRELKA_GENOME.out.vcf,
        strelka_genome_vcf.no_intervals,
        strelka_vcf.no_intervals)

    // if (tools.contains('tiddit')) {
    //     TODO: Update tiddit on bioconda, the current version does not support cram usage, needs newest version:
    //     https://github.com/SciLifeLab/TIDDIT/issues/82#issuecomment-1022103264
    //     Issue opened, either this week or end of february

    //     TIDDIT_SV(
    //         cram_recalibrated,
    //         fasta,
    //         fasta_fai
    //     )

    //     TABIX_BGZIP_TIDDIT_SV(TIDDIT_SV.out.vcf)
    //     tiddit_vcf_gz_tbi = TABIX_BGZIP_TIDDIT_SV.out.gz_tbi
    //     tiddit_ploidy     = TIDDIT_SV.out.ploidy
    //     tiddit_signals    = TIDDIT_SV.out.signals
    //     tiddit_wig        = TIDDIT_SV.out.wig
    //     tiddit_gc_wig     = TIDDIT_SV.out.gc_wig

    //     ch_versions = ch_versions.mix(TABIX_BGZIP_TIDDIT_SV.out.versions)
    //     ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)
    // }

    ch_versions = ch_versions.mix(BGZIP_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(BGZIP_DEEPVARIANT_VCF.out.versions)
    ch_versions = ch_versions.mix(BGZIP_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(BGZIP_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(BGZIP_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(BGZIP_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(BGZIP_STRELKA.out.versions)
    ch_versions = ch_versions.mix(CONCAT_GVCF_DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(CONCAT_VCF_DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(CONCAT_VCF_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(CONCAT_VCF_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(CONCAT_VCF_STRELKA.out.versions)
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)
    ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)
    ch_versions = ch_versions.mix(TABIX_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(TABIX_DEEPVARIANT_VCF.out.versions)
    ch_versions = ch_versions.mix(TABIX_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(TABIX_HAPLOTYPECALLER.out.versions)

    emit:
    deepvariant_vcf
    freebayes_vcf
    haplotypecaller_gvcf
    manta_vcf
    strelka_vcf

    versions = ch_versions
}
