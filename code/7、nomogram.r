###################################################
phen_surv_tum$age <- ifelse(phen_surv_tum$age_at_index.demographic > 60, '>60','<=60')
phen_surv_tum$stage <- gsub('.* (i+).*', '\\1', phen_surv_tum$tumor_stage.diagnoses) %>%
  gsub('not reported', NA, .)
colnames(phen_surv_tum)[4] <- 'gender'

F6_dat <- data.frame(complex_sur[,c(2:3,10)], 
                      gender=ifelse(phen_surv_tum$gender=='female',0,1),
                      stage=ifelse(phen_surv_tum$stage=='i',0,1),
              TNM=ifelse(phen_surv_tum$TNM%in% c('III','IV'), 1, 0),
                      age=ifelse(phen_surv_tum$age=='>60',1,0),
              her2=ifelse(phen_surv_tum$Subtype=='Her2', 1, 0))
F6_dat$Complex <- ifelse(F6_dat$Complex > median(F6_dat$Complex, na.rm = T), 1, 0)

univ_formulas <- sapply(colnames(F6_dat)[-(1:2)], function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
univ_models <- lapply(univ_formulas, function(x) coxph(x, data=F6_dat))
## Extract data
univ_results <- lapply(univ_models, function(x) {
  char <- names(x$assign)
  x <- summary(x)
  p.value <- round(x$wald["pvalue"], 3)
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
result <- res[!is.na(res$p.value), ]
result

suppressMessages(library(forestplot))
for_res <- rbind(c("Characteristics", NA, NA, NA, "HR (95%CI)",NA,NA,"P Value"), result)
for_res_char <- c("Characteristics", 'Complex (high vs. low)', 'gender (male vs. female)', 'stage (II&III vs. I)', 'TNM (III/IV vs I/II)', 'age (>60 vs. <=60)','Her2 (yes vs no)')

for_res$Characteristics <- for_res_char
for_res[-1,] <- for_res[-1,][order(for_res[-1,]$p.value),]
for(i in 2:4) {for_res[, i] = as.numeric(for_res[, i])}

forestplot(for_res[,c(1,5,8)],  
             mean=for_res[,2],   # 第2列为HR，森林图的小方块
             lower=for_res[,3], 
             upper=for_res[,4],  
             zero=1,            
             boxsize=0.1,       #设置小黑块的大小
             graph.pos="right",
             #mar=unit(c(2,6,6,6), "mm"),
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "8" = gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
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
             ci.vertices=T,
             ci.vertices.height=0.1, 
             new_page=T,
             #clip=c(0.1,8),
             ineheight=unit(8, 'mm'), 
             line.margin=unit(8, 'mm'),
             colgap=unit(6, 'mm'),
             col=fpColors(zero = "#e22e2a",
                          box = '#1a37d9', 
                          lines = '#1a37d9'),
             fn.ci_norm="fpDrawDiamondCI", # 菱形
             size = .2,
) #森林图放在最右侧

cox_multi <- as.formula(paste0('Surv(OS.time, OS)~', paste0(colnames(F6_dat)[-(1:2)],collapse='+')))
coxph_multi <- coxph(cox_multi, data=F6_dat)
result <- cox_multi_HR(coxph_multi)
result

for_res <- rbind(c("Characteristics", NA, NA, NA, "HR (95%CI)",NA,NA,"P Value"), result)
for_res$Characteristics <- for_res_char
for_res[-1,] <- for_res[-1,][order(for_res[-1,]$p.value),]
for(i in 2:4) {for_res[, i] = as.numeric(for_res[, i])}

forestplot(for_res[,c(1,5,8)],  
             mean=for_res[,2],   # 第2列为HR，森林图的小方块
             lower=for_res[,3], 
             upper=for_res[,4],  
             zero=1,            
             boxsize=0.1,       #设置小黑块的大小
             graph.pos="right",
             #mar=unit(c(2,6,6,6), "mm"),
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "8" = gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
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
             ci.vertices=T,
             ci.vertices.height=0.1, 
             new_page=T,
             #clip=c(0.1,8),
             ineheight=unit(8, 'mm'), 
             line.margin=unit(8, 'mm'),
             colgap=unit(6, 'mm'),
             col=fpColors(zero = "#e22e2a",
                          box = '#1a37d9', 
                          lines = '#1a37d9'),
             fn.ci_norm="fpDrawDiamondCI", # 菱形
             size = .2,
) #森林图放在最右侧


###############  列线图  ##########
suppressMessages(library(rms))
dd <- datadist(F6_dat)
options(datadist="dd")

F6_sur <- Survival(cph(cox_multi, data=F6_dat, x=T, y=T, surv=T))
F6_nom <- nomogram(F6_cph, fun=list(function(x) F6_sur(365, x),
                                function(x) F6_sur(365*3, x)),
                    fun.at=c(.05, seq(.1,.9,by=.1)),
                    funlabel=c('1 year survival', 
                               '3 year survival'))

F6_cph_1 <- cph(cox_multi, data=F6_dat, x=T, y=T, 
                 surv=T, time.inc = 365)
F6_cph_3 <- cph(cox_multi, data=F6_dat, x=T, y=T, 
                 surv=T, time.inc = 365*3)
cal_1 <- calibrate(F6_cph_1,u=365,cmethod='hare',m=30,B=1000)
cal_3 <- calibrate(F6_cph_3,u=365*3,cmethod='hare',m=30,B=1000)
##绘制1,3年生存期校准曲线
if(T){
pdf('Figure/Figure7C.pdf', width = 10, height = 8)
op <- par()
par(mar=c(2,.5,1,1))
par(fig=c(0,6,0,10)/10)
plot(F6_nom, cex.axis=.8, cex.var=.98)
par(mar=c(4,4,1,1))
par(fig=c(6.3,10,5.3,10)/10)
par(new=T)
plot(cal_1,lwd=2,lty=1, subtitles = F,
     errbar.col=c(rgb(0,118,192,
                      maxColorValue = 255)), ##设置一个颜色
     xlab='Predicted Probability of 1-year Survival',#便签
     ylab='Actual 1-year Survival',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0.9,1),ylim = c(.9,1))
par(mar=c(4,4,1,1))
par(fig=c(6.3,10,0.3,5)/10)
par(new=T)
plot(cal_3,lwd=2,lty=1, subtitles = F,
     errbar.col=c(rgb(0,118,192,
                      maxColorValue = 255)), ##设置一个颜色
     xlab='Predicted Probability of 3-year Survival',#便签
     ylab='Actual 3-year Survival',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0.6,1),ylim = c(0.6,1))
dev.off()
par(op)
}
