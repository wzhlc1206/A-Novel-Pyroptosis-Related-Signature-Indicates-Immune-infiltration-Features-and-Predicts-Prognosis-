#########
suppressMessages(library(survival))
suppressMessages(library(survminer))

identical(colnames(expr_nor), phen_surv$submitter_id.samples)
phen_surv_tum <- phen_surv[phen_surv$group == 'Tumor',]
expr_nor_tum <- expr_nor[,phen_surv$group == 'Tumor']
identical(colnames(expr_nor_tum), phen_surv_tum$submitter_id.samples)

cox_nor <- t(expr_nor_tum[deg_pyro,]) %>% as.data.frame
cox_data <- cbind(phen_surv_tum[, c('OS', 'OS.time')], cox_nor)

colnames(cox_data) <- gsub('-', '_', colnames(cox_data))
covariates <- gsub('-', '_', deg_pyro)
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS.time, OS)~', x)))

univ_models <- lapply(univ_formulas, function(x) coxph(x, data=cox_data))

## Extract data
univ_results <- lapply(univ_models, function(x) {
  char <- names(x$assign)
  x <- summary(x)
  p.value <- round(x$wald["pvalue"], 2)
  wald.test <- round(x$wald["test"], 2)
  beta <- round(x$coef[1], 2);#coeficient beta
  HR0 <- round(x$coef[2], 2);#exp(beta)
  CI5 <- round(x$conf.int[,"lower .95"], 2)
  CI95 <- round(x$conf.int[,"upper .95"],2)
  HR <- paste0(HR0, " (", CI5, "-", CI95, ")")
  res <- data.frame('Characteristics'=char, 'HR'=HR0,
                    'CI5'=CI5, 'CI95'=CI95,
                    "HR (95% CI for HR)"=HR, 'beta'=beta,
                    "wald.test"=wald.test, "p.value"=p.value,
                    check.names=F, row.names=F)
  return(res)
})

res <- do.call(rbind, univ_results)
res <- res[!is.na(res$p.value), ]
result <- res[order(res$p.value),]
result
table(result$p.value < 0.05)
write.table(result, file = "univariate_COX.txt", quote = FALSE,
            sep = "\t", row.names = TRUE, col.names = TRUE)

suppressMessages(library(forestplot))
for_res <- rbind(c("Characteristics", NA, NA, NA, "HR (95%CI)",NA,NA,"P Value"), result[result$p.value<0.05, ])
for(i in 2:4) {for_res[, i] = as.numeric(for_res[, i])}
forestplot(for_res[,c(1,5,8)],  # 列显示为变量 HR(CI) p数值形式
                                            mean=for_res[,2],   # 第2列为HR，森林图的小方块
                                            lower=for_res[,3], 
                                            upper=for_res[,4],  
                                            zero=1,            #零线或参考线为HR=1即x轴的垂直线
                                            boxsize=0.4,       #设置小黑块的大小
                                            graph.pos="right",
                                            hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                                            "2" = gpar(lty=2),
                                                            "21" = gpar(lwd=2,lty=1)),
                                            graphwidth = unit(.3,"npc"),
                                            #xlab="HR (95% CI)",
                                            #xticks=c(0,1,2,4,6,8)
                                            txt_gp=fpTxtGp(label=gpar(cex=.8),
                                                           ticks=gpar(cex=.8), 
                                                           xlab=gpar(cex=1), 
                                                           title=gpar(cex=1)),
                                            lwd.zero=1,
                                            lwd.ci=1.6,
                                            lwd.xaxis=1, 
                                            lty.ci=1,
                                            ci.vertices =T,
                                            ci.vertices.height=0.2, 
                                            #clip=c(0.1,8),
                                            #lineheight=unit(8, 'mm'), 
                                            line.margin=unit(8, 'mm'),
                                            colgap=unit(6, 'mm'),
                                            col=fpColors(zero = "#e22e2a",
                                                         box = '#1a37d9', 
                                                         lines = '#1a37d9'),
                                            fn.ci_norm="fpDrawDiamondCI", # 菱形
                                            size = .6
) #森林图放在最右侧


#########
sur_data <- cox_data[,c('OS', 'OS.time', 
                        rownames(result[result$p.value<0.05, ]))]
pl <- lapply(3:ncol(sur_data), function(x) 
  single_SUR(sur_data,x=x,cut_fun='median', plot = F))

pl_df <- do.call(rbind, pl)
table(pl_df$pval < 0.05)
cox_gene <- pl_df$Var1[pl_df$pval < 0.05]

pl <- lapply(pl_df$Var1[pl_df$pval < 0.05], function(x){
  p <- single_SUR(sur_data,x=x,cut_fun='median', risk.table = T)
  p$plot <- p$plot + theme_my_pub(legend.pos = 'top-right')
  p
})

F3_1 <- plot_grid(plotlist = pl[1:6], labels = LETTERS[2:7],
                  nrow = 2, byrow = T)
F3_1 <- plot_grid(fig, F3_1, labels = c('A'), rel_widths = c(1,3))
F3_2 <- plot_grid(plotlist = pl[7:14], labels = LETTERS[8:15],
                  nrow = 2, byrow = T)
plot_grid(F3_1, F3_2, nrow = 2)
