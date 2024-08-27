#!/usr/bin/env nextflow
// Workflow    : RNA-seq QC
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Versión     : 0.1

nextflow.enable.dsl=2

include { rna_qc_noiseq } from "/home/paulinapg/redesROSMAP/RNAseq_QC_NOISeq/QC_nextflow/run_qc.nf"

// Imprimir algunos directorios importantes
println " "
println "RNA-seq quality control"
println "Nombre del proyecto: $params.project_name"
println "Información de los datos: $params.data_info"
println "Directorio de salida: $params.out"
println " "

workflow {

// Data preprocessing
   Channel.fromPath("${params.data_info}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def ID = "${row.ID}"
                         def counts  = file("${row.counts}")
                         def metadata  = file("${row.metadata}")
                 return [ ID, counts, metadata ]
               }
          .set { set_data}
    
   rna_qc_noiseq(set_data)
}
