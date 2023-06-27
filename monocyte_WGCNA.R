HP_PEC_NToma_csv <-read.csv("/in/monocyte_scRNAseq/new_data/HP_PEC_NToma_matrix.csv", header = TRUE, row.names=1)
HP_PEC_PToma_csv <-read.csv("/in/monocyte_scRNAseq/new_data/HP_PEC_PToma_matrix.csv", header = TRUE, row.names=1)
IL4_PEC_NToma_csv <-read.csv("/in/monocyte_scRNAseq/new_data/IL4_PEC_NToma_matrix.csv", header = TRUE, row.names=1)
IL4_PEC_PToma_csv <-read.csv("/in/monocyte_scRNAseq/new_data/IL4_PEC_PToma_matrix.csv", header = TRUE, row.names=1)
SS_PEC_NToma_csv <-read.csv("/in/monocyte_scRNAseq/new_data/SS_PEC_NToma_matrix.csv", header = TRUE, row.names=1)
SS_PEC_PToma_csv <-read.csv("/in/monocyte_scRNAseq/new_data/SS_PEC_PToma_matrix.csv", header = TRUE, row.names=1)
ThioIL4PCD206_PEC_PToma_csv <-read.csv("/in/monocyte_scRNAseq/new_data/ThioIL4PCD206_PEC_PToma_matrix.csv", header = TRUE, row.names=1)
ThioPCD206_PEC_PToma_csv <-read.csv("/in/monocyte_scRNAseq/new_data/ThioPCD206_PEC_PToma_matrix.csv", header = TRUE, row.names=1)


inter_gene <- Reduce(intersect, list(row.names(HP_PEC_NToma_csv), row.names(HP_PEC_PToma_csv), 
                    row.names(IL4_PEC_NToma_csv), row.names(IL4_PEC_PToma_csv),
                    row.names(SS_PEC_NToma_csv), row.names(SS_PEC_PToma_csv),
                    row.names(ThioIL4PCD206_PEC_PToma_csv), row.names(ThioPCD206_PEC_PToma_csv),
                    row.names(HP_PEC_NToma_csv), row.names(HP_PEC_NToma_csv)))

library(Seurat)
library(SeuratObject)
HP_PEC_NToma <- CreateSeuratObject(counts = HP_PEC_NToma_csv[inter_gene,], project = "HP_PEC_NToma")
HP_PEC_PToma <- CreateSeuratObject(counts = HP_PEC_PToma_csv[inter_gene,], project = "HP_PEC_PToma")
IL4_PEC_NToma <- CreateSeuratObject(counts = IL4_PEC_NToma_csv[inter_gene,], project = "IL4_PEC_NToma")
IL4_PEC_PToma <- CreateSeuratObject(counts = IL4_PEC_PToma_csv[inter_gene,], project = "IL4_PEC_PToma")
SS_PEC_NToma <- CreateSeuratObject(counts = SS_PEC_NToma_csv[inter_gene,], project = "SS_PEC_NToma")
SS_PEC_PToma <- CreateSeuratObject(counts = SS_PEC_PToma_csv[inter_gene,], project = "SS_PEC_PToma")
ThioIL4PCD206_PEC_PToma <- CreateSeuratObject(counts = ThioIL4PCD206_PEC_PToma_csv[inter_gene,], project = "ThioIL4PCD206_PEC_PToma")
ThioPCD206_PEC_PToma <- CreateSeuratObject(counts = ThioPCD206_PEC_PToma_csv[inter_gene,], project = "ThioPCD206_PEC_PToma")
HP <- merge(x = HP_PEC_NToma,
           y = HP_PEC_PToma)
IL4 <- merge(x = IL4_PEC_NToma,
            y = IL4_PEC_PToma)
SS <- merge(x = SS_PEC_NToma,
            y = SS_PEC_PToma)
Thio <- merge(x = ThioIL4PCD206_PEC_PToma, 
              y = ThioPCD206_PEC_PToma)
all_1 <- merge(x = HP, y = IL4)
all_2 <- merge(x = SS, y = Thio)
all <- merge(x = all_1, y = all_2)

library(UCell)
all_scaled <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all_scaled <- FindVariableFeatures(object = all_scaled, selection.method = "vst")
all_scaled <- ScaleData(all_scaled, features = rownames(all_scaled))
all_scaled <- RunPCA(all_scaled, features = VariableFeatures(object = all_scaled))

genes.use <- rownames(all)
datExpr <- as.data.frame(GetAssayData(all, assay='RNA', slot='data')[genes.use,])
datExpr <- as.data.frame(t(datExpr))
all_scaled <- AddMetaData(all_scaled, metadata=datExpr)
all_scaled <- SmoothKNN(all_scaled, 
                        signature.names=c(inter_gene, c("X1700034P13Rik", "X1700006J14Rik", "X2010309G21Rik",
                                                        "H2.Oa", "H2.DMb2", "H2.Ab1", "H2.Aa", "H2.Eb1", "X3830403N18Rik")),
                        reduction="pca", k=20, suffix = "_smooth")


gene_list_smooth = c()
gene_list = c()
for (gene in c(inter_gene, c("X1700034P13Rik", "X1700006J14Rik", "X2010309G21Rik",
                             "H2.Oa", "H2.DMb2", "H2.Ab1", "H2.Aa", "H2.Eb1", "X3830403N18Rik"))){
  if (gene %in% c("1700034P13Rik", "1700006J14Rik", "2010309G21Rik",
                  "H2-Oa", "H2-DMb2", "H2-Ab1", "H2-Aa", "H2-Eb1", "3830403N18Rik")){
    print(gene)
  }else{
    gene_list_smooth = c(gene_list_smooth, paste0(gene,"_smooth"))
    gene_list = c(gene_list, gene)
  }}


datExpr_smooth <- as.data.frame(all_scaled@meta.data[gene_list_smooth])
colnames(datExpr_smooth) <- gene_list

library(WGCNA)
nclusters <- length(unique(all@active.ident))
group <- as.factor(all@active.ident)



powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr_smooth, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower= sft$powerEstimate

nSets = 1
setLabels = 'Mono'
shortLabels = setLabels
multiExpr <- list()
multiExpr[['Mono']] <- list(data=datExpr_smooth)
checkSets(multiExpr)


net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              randomSeed = 42,
                              corType = "pearson",
                              power = softPower,
                              consensusQuantile = 0.3,
                              networkType = "signed",
                              TOMType = "signed",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0.8,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 4,
                              pamStage=FALSE,
                              mergeCutHeight = 0.1,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
table(mergedColors)
moduleColors = as.character(moduleLabels)



library(caret)
dummy <- dummyVars(" ~ .", data=data.frame(all@active.ident, row.names = rownames(datExpr_smooth)))
datTraits <-  data.frame(predict(dummy, newdata=data.frame(all@active.ident, row.names = rownames(datExpr_smooth))))
colnames(datTraits) <- unique(all@active.ident)

nGenes = ncol(datExpr_smooth);
nSamples = nrow(datExpr_smooth);
MEs0 = moduleEigengenes(datExpr_smooth, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-condition relationships"))


moduleModuleCor = cor(MEs, MEs, use = "p");
moduleModulePvalue = corPvalueStudent(moduleModuleCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleModuleCor, 2), "\n(",
                   signif(moduleModulePvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleModuleCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleModuleCor,
               xLabels = names(MEs),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Module relationships"))


TOM = TOMsimilarityFromExpr(datExpr_smooth, power = softPower);
annot = read.csv(file = "/in/monocyte_scRNAseq/GeneAnnotation.csv");
head(annot)

modules = "turquoise"
probes = names(datExpr_smooth)
inModule = is.finite(match(mergedColors, modules))
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modProbes, modProbes)


saa3_top20 <- names(sort(modTOM["Saa3",])[c(1:20)])
saa3_modTOM <- modTOM[c("Saa3", saa3_top20),c("Saa3", saa3_top20)]
cyt = exportNetworkToCytoscape(saa3_modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = c("Saa3", saa3_top20))

barplot(sort(modTOM["Saa3",])[c(1:20)],
        horiz = TRUE, las = 1, main = 'Top 20 co-expressed genes of Saa3')


modules = "brown"
probes = names(datExpr_smooth)
inModule = is.finite(match(mergedColors, modules))
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modProbes, modProbes)

retnla_top20 <- names(sort(modTOM["Retnla",])[c(1:20)])
retnla_modTOM <- modTOM[c("Retnla", retnla_top20),c("Retnla", retnla_top20)]
cyt = exportNetworkToCytoscape(retnla_modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = c("Retnla", retnla_top20))
barplot(sort(modTOM["Retnla",])[c(1:20)],
        horiz = TRUE, las = 1, main = 'Top 20 co-expressed genes of Retnla')


