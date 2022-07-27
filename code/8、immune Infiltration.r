##########################
####### 数据载入
TCGA_infiltration <- readr::read_csv("infiltration_estimation_for_tcga.csv.gz") %>% as.data.frame()
rownames(TCGA_infiltration) <- TCGA_infiltration$cell_type
BRCA_infiltration <- TCGA_infiltration[unique(substring(colnames(expr_nor_tum),1,15)), grepl('cell_type|CIBERSORT$', colnames(TCGA_infiltration))]
BRCA_infiltration$Survival <- phen_surv_tum[match(rownames(BRCA_infiltration), substring(phen_surv_tum$submitter_id.samples,1,15)),'Survival',drop=T]
BRCA_infiltration %<>% reshape2::melt(id=c('cell_type','Survival'))
BRCA_infiltration$variable <- gsub('_CIBERSORT','',BRCA_infiltration$variable)
colnames(BRCA_infiltration) <- c("sample", "Survival", "cell_type", "value")

suppressMessages(library(rstatix))
F7_dat <- BRCA_infiltration

stat.test <- F7_dat %>%
  group_by(cell_type) %>%
  t_test(value ~ Survival, p.adjust.method = 'fdr')
stat.test$Type <- ifelse(stat.test$statistic > 0, 'UP', 'DOWN')
F7_dat <- rbind(data.frame(Dataset='TCGA',
                stat.test[,c('cell_type','Type', 'p')]),
                data.frame(Dataset='GSE14520',
                F7_gse_dat[,c('cell_type','Type', 'p')]))

F7A <- ggplot(F7_dat, aes(Dataset, y=cell_type)) +
  geom_point(aes(col=Type, size=-log10(p))) +
  ylab('Cell Type') +
  theme_my_pub(border = T, major = T)

BRCA_infiltration <- t(gsva_es) %>% as.data.frame() %>%
  mutate(id=substring(rownames(.), 1, 15)) %>%
  right_join(BRCA_infiltration, by=c('id'='sample')) %>%
  dplyr::select(id, Survival, cell_type, value, mySet)
colnames(BRCA_infiltration) <- c("sample", "Survival", "cell_type", "value", "Complex")

cor_cell_type <- c('B cell naive', "Macrophage M0", "T cell CD4+ memory activated", 'Mast cell resting','Myeloid dendritic cell activated')

F7_ls <- list()
for(i in cor_cell_type){
  dat <- BRCA_infiltration[BRCA_infiltration$cell_type==i, ]
  gg <- ggpubr::ggscatterhist(data = dat, x='value', y='Complex',
    size=2, margin.params = list(fill=c("#00bfc4")),
    ggtheme = theme_my_pub(border = T, major = T),
    xlab=i, add = "reg.line",
    add.params = list(color = "#f8766d", fill = "lightgray"),
    conf.int = TRUE, cor.coef = TRUE, 
    cor.coeff.args = list(method = "pearson",
                          label.x.npc = .5,
                          label.y.npc = .95,
                          label.sep = ", "),
    print = F)
  gg <- as.ggplot(~print(gg))
  F7_ls <- c(F7_ls, list(gg))
}
