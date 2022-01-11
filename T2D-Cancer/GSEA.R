setwd("C:/Users/Nayan/Desktop/Research paper/Analyzed dataset/Bladder cancer/xlsx/bladderc2_GSE37817")

# 2. Load libraries for script one
library(RCurl)
library(GEOquery)
library(limma)
library(topGO)
library(genefilter)
library(Biobase)

gene <- readxl::read_excel("bladderc2_GSE37817.xlsx")



subtable_result1 <- subset(gene, select=c("Gene.symbol","logFC","P.Value","adj.P.Val"))

#DEGs (Pval= 0.05)
pvaltab1 <- (subtable_result1[subtable_result1$P.Value<0.05,])




geneList <- subtable_result1$logFC

names(geneList) <- subtable_result1$Gene.symbol


topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}


Disease_X <- new("topGOdata",
                 description = "Disease study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "Symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 10)
# Description of the experiment
description(Disease_X)

#DEGs(Pval=0.05 and logFC=abs(1))
n_sg <- sum(topDiffGenes(geneList))
print(n_sg)

sg <- sigGenes(Disease_X)



#To see the number of annotated genes

resultFisher <- runTest(Disease_X, algorithm = "classic", statistic = "fisher")


resultKS <- runTest(Disease_X, algorithm = "classic", statistic = "ks")



allRes <- GenTable(Disease_X, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)



printGraph(Disease_X,resultFisher,firstSigNodes = 5, fn.prefix = "GSE37817", useInfo = "all", pdfSW = TRUE)
terms <- allRes$GO.ID

# print(terms)
genes <- genesInTerm(Disease_X,terms)


#To list the genes annotated to a set of specified GO terms 
str(genes)
for (i in 1:length(terms))
{
  term <- terms[i]
  genes_term <- genes[term][[1]]
  # find the genes that are in the list of genes of interest
  fact <- genes_term %in% sg
  genes_term_2 <- genes_term[fact == TRUE]
  genes_term_2 <- paste(genes_term_2, collapse=',')
  cat(term,"genes:",genes_term_2,"/n", append = TRUE, file = "bladderc2_GSE37817_correspondence.txt" )
}


write.csv(allRes,"GSE37817_top30GO.csv")

