############## 随机森林
# 数据准备
suppressMessages(library(randomForest))
suppressMessages(library(caret))


identical(colnames(expr_nor_tum), phen_surv_tum$submitter_id.samples)
cox_gene <- gsub('_', '-', cox_gene)
expr_rf <- scale(t(expr_nor_tum[cox_gene, ]))

# 训练
trControl <- trainControl(method="repeatedcv", number=10, repeats=3)
set.seed(1234)
rf_tune <- train(x=expr_rf, y=phen_surv_tum$Survival, method="rf",
                  trControl=trControl, importance=T, 
      tuneGrid=expand.grid(mtry=seq(2,ncol(expr_rf),by = 1)))

# 最终模型设定
trControl <- trainControl(method="none", classProbs = T)
set.seed(1234)
rf_final <- train(x=expr_rf, y=phen_surv_tum$Survival, method="rf",
       trControl=trControl, tuneGrid = rf_tune$bestTune,
       importance=T)

##### RF筛选合并相关基因
rf_imp <- importance(rf_final$finalModel) %>% as.data.frame()
rf_imp <- rf_imp[order(rf_imp[,'MeanDecreaseGini'],decreasing=T), ]
rf_imp$Gene <- rownames(rf_imp)

rf_gene_n <- 10
rf_imp$Top <- c(rep('yes',rf_gene_n),
                rep('no', nrow(rf_imp)-rf_gene_n))
rf_selected_gene <- rownames(rf_imp)[1:rf_gene_n]

F4A <- ggpubr::ggdotchart(rf_imp, x = "Gene",
                          y = "MeanDecreaseGini",
           color = "Top",
           #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
           sorting = "descending", #"ascending",
           add = "segments",
           dot.size = 4,
           ggtheme = theme_my_pub(),
           xlab=""
) + scale_y_continuous(expand = c(0,0,.05,.05)) +
  theme(legend.position = 'none')
F4A


# heatmap
suppressMessages(library(ComplexHeatmap))
heatmap_scale <- t(scale(t(expr_nor_tum[rf_selected_gene, ])))
left_df <- data.frame(Type=deg_sur[rf_selected_gene, ]$reg)
annotation_left <- rowAnnotation(df=left_df, 
                        col=list(Type=c('UP'='#f8766d',
                                        'DOWN'='#00ba38')), 
                       gap = unit(2, "mm"),
                       show_annotation_name = F,
                       simple_anno_size = unit(8,'pt'),
                       show_legend = T)

ht <- Heatmap(heatmap_scale, name = "mat", color_space = "RGB",
              cluster_columns = T, cluster_rows = T,
              #top_annotation = annotation_top,
              left_annotation = annotation_left,
              width = ncol(heatmap_scale)*unit(.5, "pt"),
              height = nrow(heatmap_scale)*unit(20, "pt"),
              show_column_names = F, show_row_names = T,
              row_names_gp = grid::gpar(fontsize = 8),
              show_heatmap_legend = T
)

# survival
suppressMessages(library(survival))
suppressMessages(library(survminer))

sur_gene <- gsub('-', '_', rf_selected_gene)
sur_data <- cox_data[,c('OS', 'OS.time', sur_gene)]
pl <- lapply(3:ncol(sur_data), function(x){
  p <- single_SUR(sur_data,x=x,cut_fun='median', risk.table = T)
  p$plot <- p$plot + theme_my_pub(legend.pos = 'top-right')
  p
})
wrap_plots(pl, nrow=3, byrow=T) + 
  plot_annotation(tag_levels=list(LETTERS[3:25]))

# ROC曲线绘制
suppressMessages(library(pROC))
ROC_data <- single_ROC(phen_surv_tum$Survival,
                       rf_selected_gene,
                       expr_nor_tum)

ROC_data$Gene <- paste0(ROC_data$Gene, ', AUC=',
                        round(ROC_data$AUC,3))
ggplot(data=ROC_data, mapping=aes(x=FPR, y=TPR)) +
  geom_step(aes(color=Gene), size=1, direction = "mid") +
  geom_segment(aes(x=0, xend=1, y=0, yend=1))  +  
  xlab("False positive rate") + 
  ylab("True positive rate") + coord_fixed(1) + xlim(0,1) + ylim(0,1) +
  theme_my_pub(legend.pos = c(.75,.3),
               legend.title = element_blank())
