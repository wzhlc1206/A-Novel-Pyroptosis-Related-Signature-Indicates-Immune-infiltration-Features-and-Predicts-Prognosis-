###############
suppressMessages(library(maftools))
maf <- read.maf("TCGA.BRCA.mutect.somatic.maf.gz")

maf <- subsetMaf(maf, genes = rf_selected_gene)
oncoplot(maf, top=20, drawRowBar = T, showTitle = F)

laml.titv <- titv(maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)

lollipopPlot(maf, gene = 'FREM1', 
                                refSeqID='NM_001177704',
             showMutationRate = TRUE,
             labelPos = 'all')
lollipopPlot(maf, gene = 'PAK7', 
                                showMutationRate = TRUE,
                                labelPos = 'all')

rainfallPlot(maf, detectChangePoints = T, pointSize = .6)

