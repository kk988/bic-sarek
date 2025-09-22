# bic-sarek: Changelog

All notable changes to this project made through MSK will be added here.

## v3.5.1_bic_1.0.0

### Added

- Added MusVar run fork which includes
  - Vardict caller
  - CollectHSMetrics
  - CollectAlignmentSummaryMetrics
  - CollectInsertSizeMetrics
  - BicSampleQC (report creation with all metrics)
  - BicPostPRocessing
    - NormalizeTagVCF
    - VCF2MAF
    - MergeMAFs
  - FilterMAF

- Added 2 local genomes (GRCm38_local and GRCm39_local)
- Added targets for M-IMPACT
- Set up targets for Twist

### Changed

### Fixed

### Removed

### Dependencies

| Dependency    | Old version | New version           |
| ------------- | ----------- | --------------------- |
| VardictJava   |             | 1.8.3.                |
| gatk4         |             | 4.5.0.0               |
| vcf2maf       |             | 1.6.22                |
| VEP (vcf2maf) |             | 102 (mm38) 113 (mm39) |

