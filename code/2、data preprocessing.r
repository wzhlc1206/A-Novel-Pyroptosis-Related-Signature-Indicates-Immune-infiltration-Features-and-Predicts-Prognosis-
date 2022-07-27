suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggplotify))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))


######## 数据载入
rawdata <- readr::read_tsv("TCGA-BRCA.htseq_fpkm-uq.tsv.gz")

if(T){
probeMap <- read.table("gencode.v22.annotation.gene.probeMap",sep="\t" , header = T)
gset <- rawdata %>% inner_join(probeMap, by=c("Ensembl_ID"="id")) %>%
  dplyr::select(gene, starts_with("TCGA"))

# 对重复的基因取最大值合并
if(any(duplicated(gset$gene))){
  suppressMessages(library(data.table))
  gset <- as.data.table(gset, keep.rownames=F)
  gset <- gset[,keyby=gene,lapply(.SD, max)]
  gset <- as.data.frame(gset)
  rownames(gset) <- gset$gene
  gset <- gset[,-1]
}

######## 表型及临床数据整合
phenodata <- readr::read_tsv("TCGA-BRCA.GDC_phenotype.tsv.gz")
surdata <- readr::read_tsv("TCGA-BRCA.survival.tsv")

phen_surv <- phenodata %>%
  inner_join(surdata, by=c("submitter_id.samples"="sample")) %>%
  dplyr::select(submitter_id.samples,age_at_index.demographic,gender.demographic,tumor_grade.diagnoses,tumor_stage.diagnoses,pathologic_M,pathologic_N,pathologic_T,breast_cancer_surgery_margin_status, margin_status, menopause_status,metastatic_breast_carcinoma_estrogen_receptor_status,metastatic_breast_carcinoma_lab_proc_her2_neu_immunohistochemistry_receptor_status, metastatic_breast_carcinoma_progesterone_receptor_status, person_neoplasm_cancer_status, vital_status.demographic,OS,OS.time)

phen_surv$group <- ifelse(as.numeric(substring(phen_surv$submitter_id.samples,14,15)) < 10, "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))

phen_surv$Survival <- ifelse(phen_surv$group == "Tumor" & phen_surv$OS == 1 & phen_surv$OS.time < 365*5, "Short",
                             ifelse(phen_surv$group == "Tumor", "Long", NA)) %>% factor(.,levels = c("Short","Long"))

phen_surv$TNM <- ifelse(grepl('Tis', phen_surv$pathologic_T) & grepl('N0', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M), '0', 
                        ifelse((grepl('T1', phen_surv$pathologic_T) & grepl('N0', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)) | (grepl('T0|T1', phen_surv$pathologic_T) & grepl('N1mi', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)), 'I',
                               ifelse((grepl('T0|T1', phen_surv$pathologic_T) & grepl('N1', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)) |
                                        (grepl('T2', phen_surv$pathologic_T) & grepl('N0|N1', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)) |
                                        (grepl('T3', phen_surv$pathologic_T) & grepl('N0', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)), 'II',
                                      ifelse((grepl('T0|T1|T2|T3', phen_surv$pathologic_T) & grepl('N2', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)) |
                                               (grepl('T3', phen_surv$pathologic_T) & grepl('N1', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)) |
                                               (grepl('T4', phen_surv$pathologic_T) & grepl('N0|N1|N2', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)) | (grepl('N3', phen_surv$pathologic_N) & grepl('M0', phen_surv$pathologic_M)), 'III',
                                             ifelse(grepl('M1', phen_surv$pathologic_M), 'IV', NA)))))

common_sample <- intersect(colnames(gset), phen_surv$submitter_id.samples)
phen_surv <- phen_surv[phen_surv$submitter_id.samples %in% common_sample, ]
gset <- gset[ ,common_sample]

phen_surv <- phen_surv[match(common_sample, phen_surv$submitter_id.samples), ]
phen_surv <- dplyr::select(phen_surv,submitter_id.samples,group,everything())
}


suppressMessages(library(limma))
gset <- gset[rowSums(gset>1) > 0.1*ncol(gset),]
expr_nor <- normalizeBetweenArrays(gset)

# 焦亡基因
pyro <- readLines('pyroptosis_genes.txt')
gse_pyro <- intersect(pyro, rownames(expr))

library(future.apply)
plan(multicore)
cor_df <- future.apply::future_lapply(gse_pyro,function(x){
  df <- future.apply::future_sapply(
    as.data.frame(t(expr_nor)), 
          function(i){
              test <- cor.test(i,as.numeric(expr_nor[x,]), method="pearson")
              c(signif(test$estimate, 5), signif(test$p.value, 5))
          }) %>% t()
  df <- data.frame(Gene1=x, Gene2=rownames(df),
                   cor=df[,1], p=df[,2])
})
cor_df <- do.call(rbind, cor_df)
cor_df_sub <- subset(cor_df, abs(cor) > 0.6 & p < 1e-10)
gse_pyro <- c(gse_pyro, unique(cor_df_sub$Gene2)) %>% unique()

