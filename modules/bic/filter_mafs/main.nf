process FILTER_MAFS {
    tag "filter_mafs"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/juno/bic/depot/singularity/r_xlsx_tidyverse/r_xlsx_tidyverse.simg' :
        '/juno/bic/depot/singularity/r_xlsx_tidyverse/r_xlsx_tidyverse.simg' }"

    input:
    path(rdas)

    output:
    path('*MusVarV1.xlsx'), emit: mutation_xlsx
    path('*MAF.txt')      , emit: unfiltered_maf
    path "versions.yml"   , emit: versions

    script:
    template 'filter01.R'

}