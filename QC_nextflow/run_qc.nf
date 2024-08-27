process rna_qc_noiseq {
    cache 'lenient'
    publishDir params.out + "/results", mode: 'copy'
     
    input:
    tuple val(ID), path(counts), path(metadata)

    output:                                            
    tuple path("*.rds"), path("*.png"), path("*.txt"), emit: counts_and_metadata 
    path("pca_preqc.png") , emit: pca_plot_ori
    path("*Ori.png"), emit: pre_qc_plots
    path("*Final.png"), emit: post_qc_plots
    path("pca_postqc.png"), emit: pca_plot_post
    path("*.log"), emit: R_sesion_info
   
    script:
    """
    Rscript /home/paulinapg/redesROSMAP/RNAseq_QC_NOISeq/QC_nextflow/scripts/2.pre-promRNA_parallel_copy.R ${counts} ${metadata}
    """
}
