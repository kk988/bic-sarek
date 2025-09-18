process BIC_SAMPLE_QC {
    tag "QC"
    label 'process_high' // change to low when testing is complete - using high to hang the pipeline if this fails.

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/juno/bic/depot/singularity/r_xlsx_tidyverse/r_xlsx_tidyverse.simg' :
        '/juno/bic/depot/singularity/r_xlsx_tidyverse/r_xlsx_tidyverse.simg' }"

    input:
    path reports
    path input
    path qc_control_csv

    output:
    path("qcPlot.pdf")   , emit: plot
    path("qcTable.xlsx") , emit: table
    path("versions.yml") , emit: versions

    script:


    """
    Rscript --vanilla ${moduleDir}/templates/sample_qc.R  ${input} ${qc_control_csv}
    """
}
