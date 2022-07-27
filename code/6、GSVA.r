########## GSVA
# 准备list注释或gmt注释
geneset <- list(mySet=ROC_data %>% dplyr::select(Gene, AUC) %>% 
                  distinct() %>% arrange(desc(AUC)) %>%
                  .[,1] %>% gsub(',.*', '', .) %>% .[1:5])

gsva_es <- GSVA::gsva(as.matrix(expr_nor_tum), geneset, #parallel.sz=20,
                min.sz=2, max.sz=500, verbose=F,
                BPPARAM=BiocParallel::SnowParam(20L))

tcga_sur_data <- read.table('TCGA_BRCA_survival.txt', sep="\t", header = T)
tcga_sur_data <- tcga_sur_data[,-(2:4)]

phen_surv_tum <- phen_surv_tum %>%
  mutate(nid=substring(submitter_id.samples,1,15)) %>%
  inner_join(tcga_sur_data, by=c("nid"="sample")) %>%
  dplyr::select(-nid)

identical(colnames(gsva_es), phen_surv_tum$submitter_id.samples)
complex_sur <- data.frame(phen_surv_tum[,c('submitter_id.samples',
                          'OS', 'OS.time','DSS','DSS.time',
                          'DFI','DFI.time','PFI','PFI.time')],
                          Complex=gsva_es[1,])
cox_model <- coxph(Surv(OS.time, OS)~Complex, data = complex_sur)
summary(cox_model)
fp <- complex_sur$Complex
names(fp) <- rownames(complex_sur)

mycolors <- colorRampPalette(c('#00bfc4', "white", '#f8766d'),
                             bias=1.2)(100)
ann_colors <- list(Type = c('Low'='#00bfc4', 'High'='#f8766d'))

# TCGA_Subtype
TCGA_Subtype <- read.table("../TCGASubtype.20170308.tsv.gz", sep="\t", header=T, quote='', comment.char='', check.names=F)
phen_surv_tum$Subtype <- TCGA_Subtype[match(substr(phen_surv_tum$submitter_id.samples, 1, 15), TCGA_Subtype$sampleID), 'Subtype_mRNA']

plot.fp <- list()
for(i in c('all', unique(phen_surv_tum$Subtype))){
  if(i=='all')
    sub_samp <- phen_surv_tum$submitter_id.samples
  else
    sub_samp <- phen_surv_tum$submitter_id.samples[phen_surv_tum$Subtype==i] %>% na.omit()
  
  fp0 <- fp[sub_samp]
  plot.cur <- single_SUR(complex_sur_data[sub_samp, ],
                         x=3,cut_fun='median', risk.table = T) 
  plot.cur$plot <- plot.cur$plot + 
    theme_my_pub(legend.pos = 'top-right')  +
    scale_color_manual(values = c('Low'='#00bfc4', 'High'='#f8766d'))
  
  fp0_dat <- data.frame(Patients=1:length(fp0),
                        Risk=as.numeric(sort(fp0)),
                        Type=ifelse(1:length(fp0)>=median(order(fp0)),
                                    'High', 'Low'))
  
  plot.point <- ggplot(fp0_dat,aes(x=Patients,y=Risk)) +
    geom_point(aes(col=Type)) +
    geom_vline(xintercept = median(order(fp0)), lty='dashed') +
    geom_hline(yintercept = median(fp0), lty='dashed') +
    theme_my_pub(border = T, legend.pos = 'top-left',
                 legend.title=element_blank())
  
  sur_dat <- data.frame(Patients=1:length(fp0),
                        Time=complex_sur[names(sort(fp0 )),'OS.time'],
                        e=complex_sur[names(sort(fp0)),'OS'])
  sur_dat$e <- ifelse(sur_dat$e==0,'Alive','Dead')
  
  plot.sur <- ggplot(sur_dat,aes(x=Patients,y=Time)) +
    geom_point(aes(col=e)) +
    geom_vline(xintercept = median(order(fp0)), lty='dashed') +
    theme_my_pub(border = T, legend.pos = 'top-left',
                 legend.title=element_blank())
  
  exp_dat <- expr_nor_tum[rf_selected_gene, names(sort(fp0))]
  tmp <- t(scale(t(exp_dat)))
  tmp[tmp > 2] <- 2
  tmp[tmp < -2] <- -2
  anno_df <- fp0_dat[,'Type', drop=F]
  rownames(anno_df) <- colnames(tmp)
  plot.h <- pheatmap(tmp,col= mycolors,show_colnames = F,
                     cluster_cols = F, cluster_rows = T,
                     annotation_col = anno_df,
                     annotation_colors = ann_colors, run_draw = F)
  
  roc_dat <- single_ROC(phen_surv_tum[match(sub_samp, phen_surv_tum$submitter_id.samples),]$Survival,
                        'Complex',
                        t(complex_sur[sub_samp,]))
  
  plot.roc <- ggplot(data=roc_dat, mapping=aes(x=FPR, y=TPR)) +
    geom_step(aes(color=Gene), size=1, direction = "mid") +
    geom_segment(aes(x=0, xend=1, y=0, yend=1))  +  
    xlab("FPR") + ylab("TPR") + #coord_fixed(1) +
    xlim(0,1) + ylim(0,1) +
    annotate('text', x=0.5, y=0.25,
             label=paste('AUC=', round(unique(roc_dat$AUC),3))) + 
    theme_my_pub(legend.pos = 'none',
                 legend.title = element_blank())
  
  plot.ls <- plot_grid(plot.cur, plot.point, plot.sur, plot.h, plot.roc,
                       align = 'vh', axis='lbr', ncol = 1)
  plot.fp[[i]] <- plot.ls
}
