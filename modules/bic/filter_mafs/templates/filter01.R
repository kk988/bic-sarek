
#!/usr/bin/env Rscript

#x11 = function (...) grDevices::x11(...,type='cairo')
library(tidyverse)

callerOrder <- c("mutect2", "freebayes", "strelka", "vardict")

format_maf<-function(maf) {
    maf %>%
        select(1,5,6,9,11,13,16,HGVSp_Short,t_vaf,t_alt_count,n_vaf,n_depth,CALLERS,FILTERS) %>%
        mutate(Tumor_Sample_Barcode=fct_reorder(Tumor_Sample_Barcode,as.numeric(gsub(".*RR","",Tumor_Sample_Barcode)))) %>%
        arrange(Tumor_Sample_Barcode,desc(t_vaf),Start_Position)
}

mm <- map(unlist(strsplit('${rdas}', " ")), readRDS) %>%
    map(mutate,PUBMED=as.character(PUBMED)) %>%
    bind_rows %>%
    type_convert(locale=locale(grouping_mark=""))

# remove calls where caller is "freebayes" and t_ref_count is "." and t_alt_count is "."
# - These are null calls that were not properly filtered out in previous steps
# then
# change vcf qual to numeric, assign to QUAL, change t_alt_count and n_alt_count are numeric
mm <- mm %>%
    filter(!(CALLER == "freebayes" & t_ref_count == "." & t_alt_count == ".")) %>%
    mutate(
        QUAL = as.numeric(na_if(vcf_qual, ".")),
        t_alt_count = as.numeric(t_alt_count),
        n_alt_count = as.numeric(n_alt_count)
    )


# freebayes filter - add pass in certain cases
# sort mutation data for some reason
mm <- mm %>%
    mutate(t_vaf=t_alt_count/t_depth,n_vaf=n_alt_count/n_depth) %>%
    mutate(FILTER=case_when(
                    CALLER=="freebayes" & QUAL>15000 & t_vaf>5*n_vaf ~ "PASS",
                    T ~ FILTER
                )
    ) %>%
    mutate(CALLER=factor(CALLER,levels=callerOrder)) %>%
    arrange(ETAG,CALLER,Tumor_Sample_Barcode)

# grab unchanging columns in the 21+ rows (will remove in next step)
val1Cols <- mm %>%
  select(-(1:20)) %>%
  select(where(~ n_distinct(.) == 1)) %>%
  colnames()

# remove columns that are the same for all rows
# clean up status column
# Group by ETAG and Tumor_Sample_Barcode
# Summarize the data for each group
# ungroup the data again
m1 <- mm %>%
    select(-any_of(val1Cols)) %>%
    mutate(STATUS=ifelse(is.na(STATUS),"",STATUS)) %>%
    group_by(ETAG,Tumor_Sample_Barcode) %>%
    mutate(
        N.CALL=n(),
        N.PASS=sum(FILTER=="PASS"),
        CALLERS=paste(CALLER,collapse=";"),
        FILTERS=paste(FILTER,collapse="|"),
        STATUS=paste(STATUS,collapse="|")
        ) %>%
    ungroup


#
# Grab all events for each ETAG and Tumor and order them
maf0 <- m1 %>% distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    arrange(Chromosome,Start_Position,Tumor_Sample_Barcode)

# Remove silent events and population events (dbSNP)
maf0f <- maf0 %>%
    filter(!is.na(HGVSp_Short) & !grepl("=\$",HGVSp_Short))

# high confidence evens with at least 1 PASS
maf1 <- maf0f %>%
    filter(t_vaf >= 0.05 & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    filter(N.PASS>=1)

# er, also collect the high confidence events with at least 2 PASS
maf3 <- maf1 %>% filter(N.PASS>=2)

# selecting and renaming columns we want to keep
tbl3 <- maf3 %>%
    select(Sample=Tumor_Sample_Barcode,Gene=Hugo_Symbol,
        Type=Variant_Classification,dbSNP=dbSNP_RS,Alteration=HGVSp_Short,
        MAF=t_vaf,t_depth,t_alt_count,n_depth,n_alt_count,CALLERS,FILTERS)

class(tbl3[["MAF"]])="percentage"

library(openxlsx)
wb=createWorkbook()
addWorksheet(wb,sheetName="HCEvents")
addWorksheet(wb,sheetName="HighSensMAF")
writeDataTable(wb,sheet=1,tbl3,tableStyle="none",withFilter=F)
setColWidths(wb,sheet=1,cols=1:ncol(tbl3),widths="auto")
writeDataTable(wb,sheet=2,maf1,tableStyle="none",withFilter=F)

projNo=grep("^Proj",strsplit(getwd(),"/")[[1]],value=T)
if(length(projNo)==0) {
    projNo=""
}

rFile=paste(projNo,"mutationReport","MusVarV1.xlsx", sep="_")

saveWorkbook(wb,rFile,overwrite=T)

mFile=paste(projNo,"MergedUnFiltered","MAF.txt", sep="_")
write_tsv(maf0,mFile)

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
