#############  PCA
pca_gene <- apply(expr_nor, 1, mad) %>% sort(decreasing = T) %>%
  names()
pca <- prcomp(t(expr_nor[pca_gene[1:5000], ]), scale.= T)
F1_A <- factoextra::fviz_pca_ind(pca, label="none", 
                                 habillage = phen_surv$group,
                                 addEllipses = T, ggtheme=theme_my_pub())

tum <- phen_surv[phen_surv$group == 'Tumor',]$submitter_id.samples
pca_gene <- apply(expr_nor[,tum], 1, mad) %>% sort(decreasing = T) %>%
  names()
pca <- prcomp(t(expr_nor[pca_gene[1:5000], tum]), scale.= T)
F1_B <- factoextra::fviz_pca_ind(pca, label="none", 
                                 habillage = phen_surv[phen_surv$submitter_id.samples %in% tum,]$Survival,
                                 addEllipses = T, ggtheme=theme_my_pub())

F1_A + F1_B + plot_annotation(tag_levels = 'A')


########## 差异分析
pval_cutoff = 0.05
LFC_cutoff = log2(1.5)
deg_gp <- DE_limma(expr_nor, phen_surv, 'group',
                   ref.group = 'Normal', normalization = F)
deg_gp <- DE_table(deg_gp, pval_cutoff = pval_cutoff,
                   LFC_cutoff = LFC_cutoff, padj = F)
deg_sur <- DE_limma(expr_nor[,tum], phen_surv[phen_surv$submitter_id.samples %in% tum,], 'Survival', ref.group = 'Short', normalization = F)
deg_sur <- DE_table(deg_sur, pval_cutoff = pval_cutoff,
                   LFC_cutoff = LFC_cutoff, padj = F)
deGene_gp <- deg_gp[deg_gp$reg != 'NOT',] %>% rownames()
deSur_gp <- deg_sur[deg_sur$reg != 'NOT',] %>% rownames()

F2A <- ggplot(deg_gp, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=reg)) + theme_my_pub()
my_write.table(deg_gp, file='Table/TableS1_DEG_tumor_normal.tsv')

F2B <- ggplot(deg_sur, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=reg)) + theme_my_pub()
my_write.table(deg_sur, file='Table/TableS1_DEG_survival.tsv')

F2A + F2B + plot_annotation(tag_levels = 'A')

deg_pyro <- intersect(gse_pyro, intersect(deSur_gp, deGene_gp))
deg_pyro
