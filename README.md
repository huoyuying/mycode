# mycode

##  select features from 151673-151676
DLPFC_151673 <- readRDS("F:\\DLPFC\\benchmark\\d151673.rds")
lib_size <-  ceiling(median(colSums(DLPFC_151673))/1000)*1000
DLPFC_151673 <- CreateSeuratObject(counts = DLPFC_151673@assays$Spatial@counts)
DLPFC_151673 <- NormalizeData(DLPFC_151673, normalization.method = 'RC', scale.factor = lib_size)
DLPFC_151673 <- FindVariableFeatures(DLPFC_151673, selection.method = 'vst', nfeatures = 3000)
hvg151673 <- VariableFeatures(DLPFC_151673)

DLPFC_151674 <- readRDS("F:\\DLPFC\\benchmark\\d151674.rds")
lib_size <-  ceiling(median(colSums(DLPFC_151674))/1000)*1000
DLPFC_151674 <- CreateSeuratObject(counts = DLPFC_151674@assays$Spatial@counts)
DLPFC_151674 <- NormalizeData(DLPFC_151674, normalization.method = 'RC', scale.factor = lib_size)
DLPFC_151674 <- FindVariableFeatures(DLPFC_151674, selection.method = 'vst', nfeatures = 3000)
hvg151674 <- VariableFeatures(DLPFC_151674)

DLPFC_151675 <- readRDS("F:\\DLPFC\\benchmark\\d151675.rds")
lib_size <-  ceiling(median(colSums(DLPFC_151675))/1000)*1000
DLPFC_151675 <- CreateSeuratObject(counts = DLPFC_151675@assays$Spatial@counts)
DLPFC_151675 <- NormalizeData(DLPFC_151675, normalization.method = 'RC', scale.factor = lib_size)
DLPFC_151675 <- FindVariableFeatures(DLPFC_151675, selection.method = 'vst', nfeatures = 3000)
hvg151675 <- VariableFeatures(DLPFC_151675)

DLPFC_151676 <- readRDS("F:\\DLPFC\\benchmark\\d151676.rds")
lib_size <-  ceiling(median(colSums(DLPFC_151676))/1000)*1000
DLPFC_151676 <- CreateSeuratObject(counts = DLPFC_151676@assays$Spatial@counts)
DLPFC_151676 <- NormalizeData(DLPFC_151676, normalization.method = 'RC', scale.factor = lib_size)
DLPFC_151676 <- FindVariableFeatures(DLPFC_151676, selection.method = 'vst', nfeatures = 3000)
hvg151676 <- VariableFeatures(DLPFC_151676)

list <- hvg151673
list <- append(list,hvg151674)
list <- append(list,hvg151675)
list <- append(list,hvg151676)
list <- unique(list)

list
write.csv(list,"allmergelist_151673_3.csv")

DLPFC_151673_exp <- data.frame(DLPFC_151673@assays$Spatial@counts)
merge_151673_exp <- DLPFC_151673_exp[1,]

v <- 1:7319
for(i in v){
  n <- which(grepl(list[i],rownames(DLPFC_151673_exp)))
  merge_151673_exp[i,] <- DLPFC_151673_exp[n,]
  print("row successful")
}

setwd("F:\\DLPFC\\151673\\input\\gene selection\\result")

