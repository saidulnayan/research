# Second Part : SEMANTIC SIMILIARY 

# 1. Set working directory
setwd("C:/Users/Nayan/Desktop/2nd text")

# 2. Load required libraries
library(GOSemSim)
library(readtext)
library(stringr)
library(factoextra)
library(dendextend)
library(corrplot)
library(RColorBrewer)

# 3. List of selected datasets correspondence files 
path <- c("BlaC1_GSE27448_correspondence.txt",
          "BlaC2_GSE37817_correspondence.txt",
          "BlaC3_GSE37815_correspondence.txt",
          "BlaC4_GSE42089_correspondence.txt",
          "BloC1_GSE56495_correspondence.txt",
          "BloC2_GSE61629_correspondence.txt",
          "BloC3_GSE47018_correspondence.txt",
          "BoC01_GSE36001_correspondence.txt",
          "BoC02_GSE56001_correspondence.txt",
          "BraC1_GSE50161_correspondence.txt",
          "BraC2_GSE31262_correspondence.txt",
          "BraC3_GSE19404_correspondence.txt",
          "BraC4_GSE39182_correspondence.txt",
          "BreC1_GSE157737_correspondence.txt",
          "BreC2_GSE86374_correspondence.txt",
          "BreC3_GSE70905_correspondence.txt",
          "BreC4_GSE71053_correspondence.txt",
          "BreC5_GSE72644_correspondence.txt",
          "CC001_GSE164191_correspondence.txt",
          "CC002_GSE141174_correspondence.txt",
          "CC003_GSE46824_correspondence.txt",
          "CC004_GSE37182_correspondence.txt",
          "EC001_GSE38129_correspondence.txt",
          "EC002_GSE100942_correspondence.txt",
          "EC003_GSE92396_correspondence.txt",
          "EC004_GSE29001_correspondence.txt",
          "EC005_GSE20347_correspondence.txt",
          "EC006_GSE17351_correspondence.txt",
          "EC007_GSE161533_correspondence.txt",
          "GC001_GSE13911_correspondence.txt",
          "GC002_GSE19826_correspondence.txt",
          "GC003_GSE103236_correspondence.txt",
          "GC004_GSE26899_correspondence.txt",
          "GC005_GSE81948_correspondence.txt",
          "HNC01_GSE51985_correspondence.txt",
          "HNC02_GSE74530_correspondence.txt",
          "HNC03_GSE78060_correspondence.txt",
          "LiC01_GSE101685_correspondence.txt",
          "LiC02_GSE120123_correspondence.txt",
          "LiC03_GSE112790_correspondence.txt",
          "LiC04_GSE121248_correspondence.txt",
          "LiC05_GSE69715_correspondence.txt",
          "LiC06_GSE84402_correspondence.txt",
          "LiC07_GSE89377_correspondence.txt",
          "LuC01_GSE118370_correspondence.txt",
          "LuC02_GSE10072_correspondence.txt",
          "LuC03_GSE134381_correspondence.txt",
          "LuC04_GSE103888_correspondence.txt",
          "LuC05_GSE103527_correspondence.txt",
          "OC001_GSE146553_correspondence.txt",
          "OC002_GSE29220_correspondence.txt",
          "OC003_GSE31682_correspondence.txt",
          "OC004_GSE38666_correspondence.txt",
          "OC005_GSE36668_correspondence.txt",
          "PaC01_GSE45765_correspondence.txt",
          "PaC02_GSE45757_correspondence.txt",
          "PaC03_GSE40098_correspondence.txt",
          "PrC01_GSE69223_correspondence.txt",
          "PrC02_GSE11682_correspondence.txt",
          "SC001_GSE7553_correspondence.txt",
          "SC002_GSE6520_correspondence.txt",
          "SC003_GSE53462_correspondence.txt",
          "T2D01_GSE55650_correspondence.txt",
          "T2D02_GSE16415_correspondence.txt",
          "T2D03_GSE71416_correspondence.txt",
          "T2D04_GSE64998_correspondence.txt",
          "T2D05_GSE59363_correspondence.txt",
          "T2D06_GSE29231_correspondence.txt",
          "T2D07_GSE29226_correspondence.txt",
          "T2D08_GSE29221_correspondence.txt")
files <- readtext(path)
files_split <- strsplit(files[,1],"_")
id <- unlist(lapply(files_split, function(x) paste0(x[1],"_",x[2])))

# 4. Create list of genes corresponding to each dataset
go_term <- lapply(files[,2], function(x) str_extract_all(x,"GO:.{7}"))
go_term_sub <- lapply(go_term, function(x) x[[1]][1:5])
names(go_term_sub) <- id
gene_grasp <- function(text_gene){
  aux_1 <- str_replace_all(text_gene,"GO:.{7} |\n|genes:", "")
  aux_2 <- strsplit(aux_1, " ")
  aux_3 <- strsplit(aux_2[[1]][2:6],",")
  aux_4 <- unique(unlist(aux_3))
  return(aux_4)
}
gene_term_sub <- lapply(files[,2],gene_grasp)
names(gene_term_sub) <- id

# saveRDS(gene_term_sub, file = "gene_term_sub_5.Rda")

# 5. Select the Gene ontology and prepare the annotation
hsGO <- godata('org.Hs.eg.db', ont="BP")
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)

# 6. Create semantic similarity matrices (GO terms, genes)
len <- length(id)
go_sem_sim_mat = gene_sem_sim_mat <- matrix(data = 0, nrow = len, ncol = len)
rownames(go_sem_sim_mat) = rownames(gene_sem_sim_mat) <- id
colnames(go_sem_sim_mat) = colnames(gene_sem_sim_mat) <- id

for(k in 1:len){
  for(kk in 1:len){
    go_sem_sim_mat[k,kk] <- mgoSim(go_term_sub[[k]], go_term_sub[[kk]], 
                                   semData=hsGO, measure="Wang", combine="BMA")
    gene_sem_sim_mat[k,kk] <- clusterSim(gene_term_sub[[k]], gene_term_sub[[kk]], 
                                         semData=hsGO2, measure="Wang", combine="BMA")
    cat("k =",k,"kk =",kk,"per =",
        round(((k-1)*len+kk)/len^2,digits = 3)*100,"%\n")
  }
}

#saveRDS(go_sem_sim_mat,file="go_sem_sim_mat_5.Rda")
#saveRDS(gene_sem_sim_mat,file="gene_sem_sim_mat_5.Rda")


# 7. Plot the semantic similarity matrices
library(DOSE)

pdf(file = "go_mat_5.pdf", width = 10)
simplot(go_sem_sim_mat,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=2,
        font.size=12, xlab="", ylab="")
dev.off()

pdf(file = "gene_mat_5.pdf", width = 10)
simplot(gene_sem_sim_mat,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=2,
        font.size=12, xlab="", ylab="")
dev.off() 

do_id <-  c("DOID:9352","DOID:11054","DOID:9119","DOID:3347","DOID:3068","DOID:1612",
            "DOID:9256","DOID:4914", "DOID:3717","DOID:11934","DOID:684",
            "DOID:5409","DOID:2394","DOID:1793","DOID:10283","DOID:1909")
do_ac <- c("T2D","UBC","AML","OS","GBM","BRCA","CRCA","EAC",
           "AGS","HNC","HCC","SCLC","OC","PaCa","PCa","MM")

do_sem_sim_mat <- doSim(do_id,do_id,measure="Wang")
rownames(do_sem_sim_mat) = colnames(do_sem_sim_mat) <- do_ac

#saveRDS(do_sem_sim_mat,file="do_sem_sim_mat.Rda")

pdf(file = "do_mat.pdf", width = 10)
simplot(do_sem_sim_mat,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=2,
        font.size=20, xlab="", ylab="")
dev.off()

# 8. Create KEGG Enrichment graph
library(clusterProfiler)

corrisp <- lapply(gene_term_sub,function(x) bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
corrisp <- lapply(corrisp,"[","ENTREZID")
corrisp <- lapply(corrisp, unlist)
names(corrisp) = substr(names(corrisp),1,5)


ccluster <- compareCluster(geneCluster = corrisp, fun = "enrichKEGG")
pdf(file = "enrich_KEGG_5.pdf", width = 15)
dotplot(ccluster, font.size = 9)
dev.off()

