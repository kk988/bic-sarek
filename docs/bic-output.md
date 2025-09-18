# nf-core/sarek: Output <!-- omit in toc -->

## Introduction <!-- omit in toc -->

This output information is IN ADDITION to the original output.md found in this same directory. This should only include changes made to the bic-sarek pipeline made by BIC.

### MusVar Pipeline

nf-sarek goes through the normal alignment, and variant calling steps. At the end of the variant calling steps, an additional caller is added:

[Vardict](#vardict)

As well as additional QC steps:
- [CollectHSMetrics](#collecthsmetrics)
- [CollectAlignmentSummaryMetrics](#collectalignmentsummarymetrics)
- [CollectInsertSizeMetrics](#collectinsertsizemetrics)

These QC steps as well as MarkDuplicates output (from nf-sarek) are sent to
[BicSampleQc](#bicsampleqc)

Variant calls from `strelka`, `vardict`, `mutect2`, and `freebayes` are sent through post processing steps which involves
- [NormalizeTagVCF](#normalizetagvcf)
- [VCF2MAF](#vcf2maf)
- [MergeMAFs](#mergemaf)
- [FilterMAF](#filtermaf)


### Vardict

[VardictJava](https://github.com/AstraZeneca-NGS/VarDictJava) is the Java port of the VarDict variant discovery program.

<details markdown="1">
<summary>Output files for all tumor-normal pairs</summary>

**Output directory: `{outdir}/variant_calling/vardict/<tumor_normal_pair>/`**
- `<tumor_normal_pair>.vcf.gz` - variants for tumor-normal pair
</details>

### CollectHSMetrics

[CollectHSMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard) Collects hybrid-selection (HS) metrics for a SAM or BAM file.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/metrics/`**
- `{sample}_hs_metrics.txt`
</details>

### CollectAlignmentSummaryMetrics

[CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard) Collects alignment summary metrics for a given alignment.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/metrics/`**
- `{sample}.as.txt`
</details>

### CollectInsertSizeMetrics

[CollectInsertSizeMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard) provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/metrics/`**
- `{sample}_is_metrics.txt`
- `{sample}_is_histogram.pdf`
</details>

### BicSampleQc

This is an R script that combines important metrics gathered into a plot and excel for easy visualization of sample and experiment quality

<details markdown="1">
<summary>Output files include all samples</summary>

**Output directory: `{outdir}/bic/`**
- `qcTable.xlsx`
- `qcPlot.pdf`
</details>

### NormalizeTagVCF

Set of commands normalizes the vcf files from all the various callers, and then uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to normalize indels. This script allows for steps before or after the vcf normalization to handle all of the variant callers we use in MusVar analysis.

<details markdown="1">
<summary>Output files for all tumor normal pairs</summary>

**Output directory: `{outdir}/variant_calling/{caller}/{tumor_normal_pair}/`**
- `*_tagged.vcf`
</details>

### VCF2MAF

[VCF2MAF](https://github.com/mskcc/vcf2maf) Converts vcf to maf format, and optionally annotates variants using VEP. vcf2maf and VEP versions differ depending on reference version, so feel free to look into the config file for details.

<details markdown="1">
<summary>Output files for each tumor normal pair for each caller </summary>

**Output directory: `{outdir}/vcf2maf/`**
- `{tumor_normal_pair}_{caller}.maf`
</details>

### MergeMAF

Quick script that will combine all MAF results for callers associated with tumor/normal pair

<details markdown="1">
<summary>Output files include all caller results for each sample</summary>

**Output directory: `{outdir}/variant_calling/`**
- `<tumor_normal_pair>merge.maf.tsv.gz` - Merged MAF gzipped
- `<tumor_normal_pair>merge.maf.rda` - R data for this session
</details>

### FilterMAF

Rscript to filter maf data using rdata object from [MergeMAF](#mergemaf) output. Combines output from all merged MAFs into two outputs, mutationReport excel, as well as a merged unfiltered MAF file.

<details markdown="1">
<summary>Output files include all samples</summary>

**Output directory: `{outdir}/post/reports/`**
- `*_mutationReport_MusVarV1.xlsx` - mutation report excel
- `*_MergedUnFiltered_MAF.txt` - All samples and callers maf results
</details>
