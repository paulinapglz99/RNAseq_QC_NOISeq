# 1.1prepro-mRNA.R
#Script for annotation, bias detection, correction (QC) and normalization of RNA-seq data
#By paulinapglz.99@gmail.com, adapted from https://github.com/CSB-IG/SGCCA/blob/main/prepro-mRNA.R
#For further details of how NOISeq works, go to https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
#For further details of how ComBat works, go to https://rdrr.io/bioc/sva/man/ComBat.html

#Libraries --- ---

#Here we use the NOISeq, edgeR and EDAseq packages for bias detection and correction

pacman::p_load('tidyverse', 
               'biomaRt',
               'NOISeq',
               'edgeR', 
               'EDASeq', 
               "sva",
               "ggplot2")

#Seed --- ---

set.seed(10)

#Functions  --- ---

#Homogenize metadata for sex
homogenize_exclude <- function(metadata) {
  unique_values <- unique(metadata$exclude)
  
  if ("TRUE" %in% unique_values && any(is.na(unique_values))) {
    # Caso 1: TRUE para excluir, NA para mantener
    metadata <- metadata %>% mutate(exclude = ifelse(is.na(exclude), FALSE, TRUE))
  } else {
    # Caso 2: FALSE para mantener, TRUE para excluir
    metadata <- metadata %>% mutate(exclude = ifelse(exclude == FALSE, FALSE, TRUE))
  }
  
  return(metadata)
}

#Delete duplicates in count data
del_dupl <- function(counts) {
  # Verificar genes repetidos
  repeated_values <- counts %>%
    group_by(feature) %>%
    filter(n() > 1) %>%
    distinct(feature) %>%
    pull(feature)
  
  # Ver las filas duplicadas
  repeated_rows <- counts[counts$feature %in% repeated_values, ]
  
  # Ordenar y calcular la mediana de los valores duplicados
  repeated_rows <- repeated_rows[order(repeated_rows$feature),]
  repeated_rows <- repeated_rows %>%
    group_by(feature) %>%
    summarize(across(everything(), median, na.rm = TRUE))
  
  # Eliminar las filas duplicadas y agregar las filas con la mediana calculada
  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  
  # Convertir las columnas seleccionadas a enteros
  counts <- counts %>% mutate(across(-feature, as.integer))
  
  return(counts)
}

#PCA 
pca <- function(counts, factors, title = "PCA Scatterplot coloured by sequencing batch", 
                                   save_plot = TRUE, filename = "PCA_plot.png") {
  pca_res <- prcomp(t(counts), retx = TRUE, center = TRUE, scale. = FALSE)  # prcomp for pca
  pca_df <- pca_res$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID')   # Convert to data frame

  # variance table
  variance_table <- data.frame(
    PC = 1:length(pca_res$sdev),
    Variance_Percentage = pca_res$sdev^2 / sum(pca_res$sdev^2) * 100,
    cumulative_percentage = cumsum(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100))
  
  # plot
  plot <- pca_df %>% 
    ggplot() +
    aes(x = PC1, y = PC2, 
        colour = as.factor(factors$sequencingBatch)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID)) +
    labs(title = title,
         subtitle = "PC1 vs PC2", 
         x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
    theme_minimal()
  
  # Guardar el gráfico si save_plot es TRUE
  if (save_plot) {
    ggsave(filename = filename, plot = plot, width = 10, height = 8, dpi = 300)
  }
  
  # Retornar el dataframe del PCA y la tabla de varianza por si se necesitan
  return(list(pca_df = pca_df, variance_table = variance_table, plot = plot))
}

#Get the data --- ---

# #Read counts data, this was already filtered by 1.QC_pre_analysis.R
# 
# counts <- readRDS(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.rds') %>%
#   as.data.frame()
# dim(counts)
# 
# #Obtain factors from metadata --- --- 
# 
# metadata <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
# dim(metadata)

args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
metadata_file <- args[2]

#Fix metadata

metadata <- homogenize_exclude(metadata)

#Replace msex with sex
names(metadata) <- ifelse(names(metadata) == "msex", "sex", names(metadata))

#Exclude specimens by metadata

metadata <- metadata %>% filter(exclude == FALSE)
metadata <- metadata %>% filter(if_all(c(ceradsc, sex), ~ !is.na(.)))
dim(metadata)

#Acondicionar conteos

counts <- counts[-c(1:4), ]
dim(counts)

#Delete duplicates

counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"))

counts <- del_dupl(counts)
dim(counts)

#Counts

counts <- counts %>% dplyr::select(1, all_of(metadata$specimenID))
dim(counts)

#Set names to rows 
rownames(counts) <- NULL
counts <- counts %>% column_to_rownames(var = "feature")

#I do this to make sure the rowlength of factors match with the counts columns, it is like metadata filtering

factors <- data.frame("specimenID" = colnames(counts))   

#Left join with metadata
 
factors <- factors %>% left_join(metadata, by = "specimenID") %>%
  dplyr::select(c(specimenID, ceradsc, sex, sequencingBatch))
dim(factors)

#Filter again counts

counts <- counts %>% dplyr::select(1,intersect(colnames(counts), factors$specimenID))
dim(counts)

#Annotation with ensembl --- ---

#Generate mart object

mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")

#We create myannot, with GC content, biotype, info for length & names per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  rownames(counts),  #annotate the genes in the count matrix 
                 mart = mart)

#to have control
myannot <- myannot %>%  rename(feature = ensembl_gene_id)

#Add length column

myannot$length <- abs(myannot$end_position-myannot$start_position)
dim(myannot)

#Create NOISeq object bias detect and bias correction --- ---

#For NOISeq, order of factors$specimenIDs and  colnames(coutns) must match
identical(colnames(counts), factors$specimenID)
#[1] TRUE

#Names of features characteristics

mylength <- setNames(myannot$length, myannot$feature)

mygc <- setNames(myannot$percentage_gene_gc_content, myannot$feature)

mybiotype <-setNames(myannot$gene_biotype, myannot$feature)

#Create NOISeq object

noiseqData <- NOISeq::readData(data = counts,
                               factors = factors,           #variables indicating the experimental group for each sample
                               gc = mygc,                   #%GC in myannot
                               biotype = mybiotype,         #biotype
                               length =  mylength)          #gene length

# 0) Diagnostic of data

#each sample "s" is compared to a reference "r" (which can be arbitrarily chosen).
#by computing M values=log2(counts = countsr). 

#Confidence intervals for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 
#"cd" means "Cumulative Distribution."

mycd <- dat(noiseqData, type = "cd", norm = F) #slow

#[1] "Warning: 6321 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 01_120405"

#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])

#1)check expression bias for counts
#This object contains the count distribution for each biological group and also the 
#percentage of features with counts per million higher than 0, 1, 2, 5 or 10, per each 
#sample independently and in at least one of the samples (total). 

#Use a same factor variable to all bias detection

mycountsbio <- dat(noiseqData, 
                   type =  "countsbio",  
                   norm = F,      #T when already normalized counts as input
                   factor = NULL) #When NULL, all factors are considered
#Plots

png("CountsOri.png")
explo.plot(mycountsbio,
           plottype = "boxplot", #type of plot
           samples = 1:20)  #only showing the first 15 samples
dev.off()

#2)check for low count genes

png("lowcountsOri.png")
explo.plot(mycountsbio,
           plottype = "barplot",
           samples = 1:20)
dev.off()

#3)check for transcript composition bias

#Plot for Mvalues

png("MvaluesOri.png")
explo.plot(mycd,samples=sample(1:100,10))
dev.off()

#4)check for length & GC bias

#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature

myGCcontent <- dat(noiseqData,
                   k = 0,            #A feature is considered to be detected if the corresponding number of read counts is > k. 
                   type = "GCbias", 
                   factor = NULL)

png("GCbiasOri.png",width=1000)
explo.plot(myGCcontent,
           samples = 1:12,   
           toplot = "global")
dev.off()

#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation

mylengthbias <- dat(noiseqData, 
                    k = 0,
                    type = "lengthbias",
                    factor = NULL)

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."

#Plot length bias

png("lengthbiasOri.png", width=1000)
explo.plot(mylengthbias, 
           samples = 1:12, 
           toplot = "global")
dev.off()

#PCA pre-QC --- ---

pca.pre <- pca(counts, factors, title = "PCA Scatterplot coloured by sequencing batch - preQC", 
               save_plot = TRUE, filename = "PCA_plot_preQC.png")

# Acceder al gráfico y mostrarlo si es necesario
pca.pre$plot

#################SOLVE BIASES###################################

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 

#Filter Genes By Expression Level
#Determine which genes have sufficiently large counts to be retained in a statistical analysis.

countMatrixFiltered <- filtered.data(counts,
                                     factor = "ceradsc",       #using all factors
                                     norm = F,            #counts are not normalized 
                                     method = 1,          #Method 1 (CPM) removes those features that have an average expression per condition less than cpm value and a coefficient of variation per condition higher than cv.cutoff (in percentage) in all the conditions
                                     cpm = 10,             #Cutoff for the counts per million value
                                     p.adj = "fdr")       #Method for the multiple testing correction
dim(countMatrixFiltered)

#Create new annotation with filtered matrix

myannot <- myannot %>% filter(feature %in% rownames(countMatrixFiltered))
dim(myannot)

myannot <- myannot[!duplicated(myannot$feature), ]
dim(myannot)

#Not all genes are annotated in the matrix

countMatrixFiltered <- countMatrixFiltered[rownames(countMatrixFiltered) %in% myannot$feature,]
dim(countMatrixFiltered)

#Feature data

featureData <-  data.frame("feature" = rownames(countMatrixFiltered)) %>% 
  left_join(myannot, by = "feature")  %>% column_to_rownames(var = "feature")
dim(featureData)

#Pheno data

phenoData <- data.frame("specimenID" = colnames(countMatrixFiltered)) %>% 
  left_join(factors, by = "specimenID") %>% column_to_rownames(var = "specimenID")
dim(phenoData)

##Create EDA object --- --

mydataEDA <- newSeqExpressionSet(
  counts = as.matrix(countMatrixFiltered),
  featureData = featureData,
  phenoData = phenoData)

#If you're re-normalizing, there's a warning, 
#Warning message:
#  In validityMethod(object) : 'counts' contains non-integer numbers

#order for less bias

#for gc content
gcFull <- withinLaneNormalization(mydataEDA, 
                                  "percentage_gene_gc_content",
                                  which = "full")#corrects GC bias 
# 
# #for length
# lFull <- withinLaneNormalization(gcFull, 
#                                  "length", 
#                                  which = "full")#corrects length bias 

# Normalization --- ---
#TMM normalization adjusts library sizes based on the assumption that most genes are not differentially expressed.

norm_count <- NOISeq::tmm(normCounts(gcFull),
                          long = 1000,  # If long == 1000, no length correction is applied (no matter the value of parameter lc). 
                          lc = 0, # If lc = 0, no length correction is applied
                          k = 0) # By default, k = 0. 

noiseqData_norm_count <- NOISeq::readData(data = norm_count, 
                                          factors = factors)

#cd Diagnostic test for length and gc correction

mycd_lessbias <- NOISeq::dat(noiseqData_norm_count,
                             type = "cd",
                             norm = TRUE)

#Table diagnostic

table(mycd_lessbias@dat$DiagnosticTest[,  "Diagnostic Test"])

#Solve batch effect --- --- 

lessbatch <- ComBat_seq(counts = exprs(noiseqData_norm_count), batch = factors$sequencingBatch)

#############################FINAL QUALITY CHECK#######################################################

pca.post <- pca(lessbatch, factors = factors, 
                title = "PCA Scatterplot coloured by library batch, post batch correction",
                save_plot = TRUE, filename = "PCA_plot_post.png")

#Create new noiseq object with normalized counts 

noiseqData_final <- NOISeq::readData(lessbatch,
                                     factors = factors,
                                     gc = mygc,
                                     biotype = mybiotype,
                                     length = mylength)

#Check for bias with renormalized

mycountsbio_final <- dat(noiseqData_final, 
                         type = "countsbio", 
                         norm=T, 
                         factor = NULL)

#Plot final plots

png("CountsFinal.png")
explo.plot(mycountsbio_final,
           plottype = "boxplot",
           samples = 1:15)
dev.off()

#Low counts 

png("lowcountsFinal.png")
explo.plot(mycountsbio_final,
           plottype = "barplot",
           samples = 1:20)
dev.off()

#Plot for Mvalues

png("MvaluesFinal.png")
explo.plot(mycd_lessbias,
           samples=sample(1:ncol(counts),10))
dev.off()

#calculate final GC bias

myGCcontent_final <- dat(noiseqData_final,
                         k = 0, 
                         type = "GCbias", 
                         factor = NULL,
                         norm = T)

#Plot final GC bias

png("GCbiasFinal.png",width=1000)
explo.plot(myGCcontent_final, 
           plottype = "boxplot",
           samples = 1:12)
dev.off()

#calculate final length bias

mylenBias <- dat(noiseqData_final, 
                 k = 0, 
                 type = "lengthbias", 
                 factor = NULL,
                 norm=T)

#Plot final length bias

png("lengthbiasFinal.png",width=1000)
explo.plot(mylenBias, samples = 1:12)
dev.off()

#Final count matrix --- ---

final_counts <- lessbatch
dim(final_counts)

#Final metadata --- ---

final_metadata <- data.frame(
  "specimenID" = colnames(final_counts))   

final_metadata <- final_metadata %>% 
  left_join(metadata, by = "specimenID")
dim(final_metadata)

#Finally, save counts table --- ---

saveRDS(final_counts, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_filteredQC_counts_DLPFC.rds")

#Save filtered metadata --- ---

vroom::vroom_write(final_metadata, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_filteredQC_DLPFC.txt")

#Save annotation --- --- 

vroom::vroom_write(myannot, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_filteredQC_annotation_DLPFC.txt")

#END