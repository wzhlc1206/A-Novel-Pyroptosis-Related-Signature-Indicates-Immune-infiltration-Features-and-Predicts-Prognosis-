suppressMessages(library(clusterProfiler))
#基因名称转换，返回的是数据框
gene <- bitr(rf_selected_gene, fromType="SYMBOL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")

# 根据相关gene的功能富集分析推测该gene的功能
go_bp <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="bp")
kegg_en <- enrichKEGG(gene = gene$ENTREZID, organism="hsa",
                      minGSSize=1, keyType="kegg",
                      use_internal_data=T)
enrich_df <- rbind(go_bp@result[c(1:7,9),],
                   kegg_en@result[c(1,2),])
F8A <- ggplot(enrich_df, aes(GeneRatio, Description)) +
  geom_point(aes(col=p.adjust, size=Count)) +
  scale_size(range = c(3,6)) +
  scale_color_gradient(low="red", high="blue",
                       guide = guide_colorbar(reverse=T)) +
  ylab(NULL) + 
  theme_my_pub(border = T, major = T)+
  theme(axis.text.y = element_text(size = 12))

# GSEA
gsea_gene_list <- deg_gp$logFC
names(gsea_gene_list) <- rownames(deg_gp)
gsea_gene_list <- sort(gsea_gene_list, decreasing = T)

geneset <- read.gmt("~/ref/msigdb/msigdb.v7.4.symbols.gmt")
geneset <- geneset[grepl('^GOBP_', geneset$term), ]
geneset$term <- gsub('^GOBP_', '', geneset$term)
geneset$term <- gsub('_', ' ', geneset$term)
geneset$term <- stringr::str_wrap(geneset$term, width = 32)

gsea_res <- GSEA(gsea_gene_list, TERM2GENE=geneset)
F8B <- enrichplot::gseaplot2(gsea_res, 1:5, base_size = 10)

go_bp_term <- go_bp@result$Description
gsea_term <- gsea_res@result[gsea_res@result$p.adjust < 0.1, ]$Description %>%
  gsub('\n','',.) %>% tolower()

merge_term <- intersect(go_bp_term, gsea_term)
merge_term_id <- go_bp@result[go_bp@result$Description %in% merge_term, ]$ID
go_id_term <- go_bp@result[go_bp@result$Description %in% merge_term, 2, drop=F]

write.table(go_bp@result, 'Table/TableS2_GO_BP.tsv', sep = '\t',
            col.names = T, row.names = F, quote = F)
write.table(kegg_en@result, 'Table/TableS2_KEGG.tsv', sep = '\t',
            col.names = T, row.names = F, quote = F)
write.table(gsea_res@result, 'Table/TableS2_GSEA.tsv', sep = '\t',
            col.names = T, row.names = F, quote = F)


suppressMessages(library(GOSemSim))
d <- godata('org.Hs.eg.db',ont="BP",computeIC=FALSE)
go_sim_df <- expand.grid(merge_term_id, merge_term_id)
go_sim_df$sim <- apply(go_sim_df, 1, function(x)
  GOSemSim::goSim(x[1], x[2], semData=d))
go_sim_df <- reshape2::dcast(go_sim_df, Var1 ~ Var2,
                   value.var="sim")
rownames(go_sim_df) <- go_sim_df$Var1
go_sim_df <- go_sim_df[,-1]
go_dist <- dist(go_sim_df, method="euclidean")
go_hclust <- hclust(go_dist, method="complete")

suppressMessages(library(ggtree))
den <- as.dendrogram(go_hclust)
clus <- cutree(go_hclust, 3)
go_g <- split(names(clus), clus)
go_p <- ggtree(go_hclust)
clades <- sapply(go_g, function(n) MRCA(go_p, n))
go_p <- groupClade(go_p, clades, group_name='subtree') +
  aes(color=subtree)
go_d <- data.frame(label = names(clus), 
                Term = go_id_term[names(clus),])

F8C <- go_p %<+% go_d + 
  geom_tiplab(aes(label=Term, col=subtree), size=4, hjust=0) +
  scale_color_brewer(palette='Set1', breaks=1:4) +
  coord_cartesian(clip = 'off') +
  theme(legend.position='none',
        plot.margin=margin(6,150,6,6))

