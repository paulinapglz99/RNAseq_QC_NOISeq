process rna_qc_noiseq {
    cache 'lenient'
    publishDir params.out, mode: 'copy'
     
    input:
    tuple val(ID), path(counts), path(metadata)

    output:                                            
    tuple path("*.rds"), path("*.png"), path("*.txt"), emit: cuack 
   
    script:
    """
    Rscript /home/paulinapg/redesROSMAP/QC_nextflow/scripts/2.pre-promRNA_parallel_copy.R ${counts} ${metadata}
    """
}
