rm(list = ls())

suppressMessages(library(GEOquery))

##### 数据载入
options('download.file.method.GEOquery' = 'libcurl')
gse_1 <- getGEO('GSEXXXX',GSEMatrix=T, destdir = "./")  # GSE

