#!/usr/bin/env Rscript

library(tidyverse)

merged_mafs=map(unlist(strsplit("${mafs}", " ")),read_tsv,col_types=cols(.default="c"),comment="#") %>%
  bind_rows %>%
  type_convert(locale=locale(grouping_mark="")) %>%
  mutate(Chromosome=factor(Chromosome,levels=c(1:19,"X","Y","MT"))) %>%
  arrange(Chromosome,Start_Position,CALLER) %>%
  mutate(ETAG=paste0(
            Chromosome,":",Start_Position,":",End_Position,":",
            Reference_Allele,":",Tumor_Seq_Allele2
            )
        )

saveRDS(merged_mafs,paste0('${meta.id}',"merge.maf.rda"),compress=T)
write_tsv(merged_mafs,paste0('${meta.id}',"merge.maf.tsv.gz"))

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
