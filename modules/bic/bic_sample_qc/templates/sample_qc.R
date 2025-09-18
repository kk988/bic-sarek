require(tidyverse)

# constants
as.metrics=c("SAMPLE", "PCT_PF_READS_ALIGNED", "MEAN_ALIGNED_READ_LENGTH",
    "PCT_READS_ALIGNED_IN_PAIRS", "PCT_PF_READS_IMPROPER_PAIRS",
    "STRAND_BALANCE", "PCT_CHIMERAS", "PCT_ADAPTER", "PCT_SOFTCLIP")

md.metrics=c("SAMPLE","PERCENT_DUPLICATION")

hs.metrics=c("SAMPLE",
"FOLD_ENRICHMENT",
"MEAN_TARGET_COVERAGE",
"PCT_TARGET_BASES_100X",
"PCT_USABLE_BASES_ON_TARGET",
"ZERO_CVG_TARGETS_PCT")

metrics0=c(
    "FOLD_ENRICHMENT", "MEAN_TARGET_COVERAGE",
    "PCT_CHIMERAS", "PCT_PF_READS_ALIGNED", "PCT_PF_READS_IMPROPER_PAIRS",
    "PCT_READS_ALIGNED_IN_PAIRS", "PCT_SOFTCLIP", "PCT_TARGET_BASES_100X",
    "PERCENT_DUPLICATION", "ZERO_CVG_TARGETS_PCT"
)

boxWidth=.4

machines=list(
    NovaSeq6000=c("DIANA","MICHELLE","RUTH"),
    NovaSeqX=c("BONO","FAUCI2")
)

# grab command line arguments
args=commandArgs(trailingOnly=T)
input_file=args[1]
qc_control_csv = args[2]

# functions
read_as_metrics<-function(ff) {
    read_tsv(ff,comment="#",show_col_types=F,progress=F) %>%
        filter(CATEGORY=="PAIR") %>%
        select(SAMPLE,everything()) %>%
        select(-LIBRARY,-READ_GROUP,-PCT_HARDCLIP)
}

read_md_metrics<-function(ff) {
    read_tsv(ff,comment="#",show_col_types=F,progress=F,n_max=1, col_types=cols(.default="c")) %>%
        rename(SAMPLE=LIBRARY)
}

read_hs_metrics<-function(ff) {
    sid=basename(ff) %>% str_replace("_hs_metrics.txt","")
    read_tsv(ff,comment="#",show_col_types=F,progress=F,n_max=1,col_types=cols(.default="c")) %>%
        mutate(SAMPLE=sid)
}

quibble2 <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}

get_qc_table<-function(tbl) {
    qvals=c(.05/2,.25,.5,.75,(1-.05/2))
    tbl %>%
        select(all_of(names(which(map_vec(tbl,class)=="numeric")))) %>%
        gather(Metric,V) %>%
        group_by(Metric) %>%
        reframe(quibble2(V,qvals)) %>%
        mutate(V_q=paste("q",V_q, sep="_")) %>%
        mutate(V_q=factor(V_q,levels=paste("q",qvals, sep="_"))) %>%
        spread(V_q,V)#        filter(q_0<q_1 & q_0<q_0.25 & q_1>q_0.75)
}

# if .INCLUDED is true, halt to avoid re-sourcing
if(exists(".INCLUDED") && .INCLUDED ) halt(".#INCLUDE")
.INCLUDED=TRUE

################################################
# start!

# manipulate input file to determine machine, machineType projId, and tumor status, etc
manifest=read_csv(input_file) %>%
    mutate(igoId=str_extract(fastq_1,"_IGO_([^/]*)/",group=1)) %>%
    mutate(machine=str_extract(fastq_1,"FASTQ.([^_]+)_",group=1)) %>%
    select(-fastq_1,-fastq_2,-lane) %>%
    mutate(machType=case_when(machine %in% machines$NovaSeq6000 ~ "NovaSeq6000", machine %in% machines$NovaSeqX ~ "NovaSeqX", T ~ "Unknown")) %>%
    mutate(projId=gsub("_\\d+","",igoId)) %>%
    mutate(sid=paste(patient,sample,sep="_")) %>%
    distinct %>%
    mutate(status=ifelse(status==0,"Normal","Tumor"))

# grab aligntment summary metrics
asm=fs::dir_ls(".", recur=T,regex=".as.txt$") %>%
    map(read_as_metrics,.progress=T) %>%
    bind_rows %>%
    select(all_of(as.metrics)) %>%
    gather(Metric,Value,-SAMPLE) %>%
    rename(sid=SAMPLE) %>%
    left_join(manifest)

# grab markduplicates metrics
mdm=fs::dir_ls(".",recur=T,regex=".metrics$") %>%
    map(read_md_metrics,.progress=T) %>%
    bind_rows %>%
    type_convert %>%
    select(all_of(md.metrics)) %>%
    gather(Metric,Value,-SAMPLE) %>%
    rename(sample=SAMPLE) %>%
    left_join(manifest)

# grab hs metrics
hsm=fs::dir_ls(".",recur=T,regex="_hs_metrics.txt$") %>%
    map(read_hs_metrics,.progress=T) %>%
    bind_rows %>%
    type_convert %>%
    select(all_of(hs.metrics)) %>%
    gather(Metric,Value,-SAMPLE) %>%
    rename(sample=SAMPLE) %>%
    left_join(manifest)



qcControl=read_csv(qc_control_csv) %>% filter(Metric %in% metrics0) %>% gather(Q,V,-Metric)
qcControl=qcControl %>% mutate(QQ=as.numeric(gsub("q_","",Q))) %>% mutate(OO=abs(log(QQ/(1-QQ)))/3+.5)

dq=bind_rows(list(asm,mdm,hsm)) %>% filter(Metric %in% metrics0) %>% distinct(sid,Metric,.keep_all=T)
pg=dq %>% ggplot(aes(status,Value,color=status)) + theme_light(16) + geom_hline(aes(yintercept=V,color=Q),data=qcControl,linewidth=qcControl$OO,alpha=.3) + geom_boxplot(outlier.shape=NA,width=boxWidth) + facet_wrap(~Metric,scale="free") + geom_jitter(width=boxWidth/2,size=3,alpha=.5) + scale_color_manual(values=c("darkblue","darkred","grey30","darkred","grey30","darkred","darkgreen")) + theme(panel.grid.minor = element_blank(),panel.grid.major.x=element_blank())
pdf(file="qcPlot.pdf",width=1.5*14,height=1.5*8.5); print(pg); dev.off()

qcTbl=dq %>% select(patient,sample,status,Metric,Value) %>% spread(Metric,Value) %>% arrange(patient,status)
openxlsx::write.xlsx(qcTbl,"qcTable.xlsx")

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

