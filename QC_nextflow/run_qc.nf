process rna_qc_noiseq {
    cache 'lenient'
    publishDir params.out + "/results", mode: 'copy'
   //errorStrategy 'ignore'
     
    input:
    tuple val(ID), path(counts), path(metadata)

    output:                                            
    path("${ID}/pca_preqc.png"), emit: pca_plot_ori
    path("${ID}/*Ori.png"), emit: pre_qc_plots
    path("${ID}/*Final.png"), emit: post_qc_plots
    path("${ID}/pca_postqc.png"), emit: pca_plot_post
    path("${ID}/*")

    script:
    """
    mkdir ${ID}
    cd ${ID}
    Rscript /home/paulinapg/redesROSMAP/RNAseq_QC_NOISeq/QC_nextflow/scripts/2.pre-promRNA_parallel.R ${counts} ${metadata}
    """
}
