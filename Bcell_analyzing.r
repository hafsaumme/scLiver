

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
library(stringr)
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("ggplot2")
library("ggpubr")


DefaultAssay(healthy_cirrhotic_Plasma) <- "RNA"
dim(healthy_cirrhotic_Plasma)
healthy_cirrhotic_Plasma$sample=healthy_cirrhotic_Plasma$Liver
healthy_cirrhotic_Plasma <- SCTransform(healthy_cirrhotic_Plasma, method = "glmGamPoi", 
                                        vars.to.regress = "percent.mt",variable.features.n = 3000,
                                        verbose = FALSE)

DefaultAssay(obj.bcell.hcc) <- "RNA"
obj.bcell.hcc$sample=obj.bcell.hcc$liver
obj.bcell.hcc <- SCTransform(obj.bcell.hcc, method = "glmGamPoi", 
                             vars.to.regress = "percent.mt",variable.features.n = 3000,
                             verbose = FALSE)
dim(obj.bcell.hcc)


obj.bcell.icca$sample=obj.bcell.icca$liver
DefaultAssay(obj.bcell.icca) <- "RNA"
obj.bcell.icca <- SCTransform(obj.bcell.icca, method = "glmGamPoi", 
                              vars.to.regress = "percent.mt",variable.features.n = 3000,
                              verbose = FALSE)
dim(obj.bcell.icca) 



healthy_cirrhotic_Bcell$sample=healthy_cirrhotic_Bcell$Liver
DefaultAssay(healthy_cirrhotic_Bcell) <- "RNA"
dim(healthy_cirrhotic_Bcell) 

obj.bcell.cirrhotic<- subset(healthy_cirrhotic_Bcell, subset = orig.ident == 'cirrhotic')
dim(obj.bcell.cirrhotic) 
obj.bcell.cirrhotic <- SCTransform(obj.bcell.cirrhotic, method = "glmGamPoi", 
                                   vars.to.regress = "percent.mt",variable.features.n = 3000,
                                   verbose = FALSE)

obj.bcell.healthy<- subset(healthy_cirrhotic_Bcell, subset = orig.ident == 'healthy')
dim(obj.bcell.healthy)
obj.bcell.healthy <- SCTransform(obj.bcell.healthy, method = "glmGamPoi", 
                                 vars.to.regress = "percent.mt",variable.features.n = 3000,
                                 verbose = FALSE)

# Read RDS file of Bcells object or do integration all bcells from four source

#RDS file of Bcells

obj.all.four=(readRDS("obj.all.four.bcell.rds"))

#....................INTEGRATION on hcc , icca, healthy and cirrhotic datasets...............

list= c(healthy_cirrhotic_Plasma, obj.bcell.healthy,obj.bcell.cirrhotic,
        obj.bcell.icca,obj.bcell.hcc)


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 4000)

list <- PrepSCTIntegration(object.list = list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
obj.all.four <- IntegrateData(anchorset = anchors,normalization.method = "SCT")

DefaultAssay(obj.all.four ) <- "integrated"



obj.all.four@meta.data$Liver=NULL
obj.all.four@meta.data$liver=NULL
obj.all.four@meta.data$integrated_snn_res.0.5=NULL
obj.all.four@meta.data$seurat_clusters=NULL
obj.all.four@meta.data$newCellName=NULL

obj.all.four@meta.data$orig.ident[which(str_detect(obj.all.four@meta.data$orig.ident,"HCC"))] <- 'HCC'
obj.all.four@meta.data$orig.ident[which(str_detect(obj.all.four@meta.data$orig.ident,"iCCA"))] <- 'iCCA'
obj.all.four@meta.data$customclassif[which(str_detect(obj.all.four@meta.data$Type,"B cell"))] <- 'Bcell'
obj.all.four@meta.data$Type=NULL


#Dimensionality Reduction

obj.all.four=RunPCA(obj.all.four, features = VariableFeatures(object = obj.all.four))

ElbowPlot(obj.all.four, ndims = 50)

pct <- obj.all.four[["pca"]]@stdev / sum(obj.all.four[["pca"]]@stdev) * 100

cumu <- cumsum(pct)

co1 <- which(cumu > 90 & pct < 5)[1]

co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
pcs <- min(co1, co2)

pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

library(ggplot2)
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()



obj.all.four=RunUMAP(obj.all.four, reduction = "pca", dims = 1:pcs)
DimPlot(obj.all.four, reduction = "umap", group.by = "orig.ident")

obj.all.four <- FindNeighbors(obj.all.four, reduction = "pca", dims = 1:pcs)

obj.all.four= FindClusters(obj.all.four,resolution = 0.2)

DimPlot(obj.all.four, reduction = "umap")

DefaultAssay(obj.all.four ) <- "RNA"

obj.all.four@assays$RNA@data=NormalizeData(obj.all.four@assays$RNA@counts, normalization.method = "LogNormalize", scale.factor = 10000)

obj.all.four.genes <- rownames(obj.all.four)

obj.all.four@assays$RNA@scale.data= ScaleData(obj.all.four@assays$RNA@data, features = obj.all.four.genes)


#==============sctype====================

library(HGNChelper)
library(openxlsx)

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "materials/hafsa/reproduce/ScTypeDB_full.xlsx"; tissue = "Liver_Bcell"

gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = obj.all.four@assays$RNA@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 


# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(obj.all.four@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(obj.all.four@meta.data[obj.all.four@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj.all.four@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


obj.all.four@meta.data$customclassif = NA
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  obj.all.four@meta.data$customclassif[obj.all.four@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


DimPlot(obj.all.four, reduction = "umap", label = TRUE, repel = TRUE) 

#do integration of B cells from all sources or read prepared RDS file from directory
obj.all.four=(readRDS("./rstudio-export/obj.all.four.bcell.rds"))
# Rename all identities
obj.all.four <- RenameIdents(object = obj.all.four, 
                             "0" = "Naïve Bcells",
                             "1" = "Naïve Bcells",
                             "2" = "Memory Bcells",
                             "3" = "Plasma cells",
                             "4" = "Memory Bcells",
                             "5" = "Plasma cells",
                             "6" = "Memory Bcells"
)

#DimPlot

DimPlot(obj.all.four, reduction = "umap", label = F, repel = T)+
  theme(axis.title = element_text(size=11))+ggtitle("")+labs(color="")

DimPlot(obj.all.four, reduction = "umap", group.by="orig.ident", label = F, repel = T)+
  theme(axis.title = element_text(size=11))+ggtitle("")+labs(color="")


#Differential Gene

DefaultAssay(obj.all.four) <- "RNA"

obj.all.four_marker=FindAllMarkers(obj.all.four, only.pos = TRUE, test.use = "wilcox",logfc.threshold = 0.25)

#DotPlot of known marker genes

DotPlot(object=obj.all.four,features = c("CD79A","IGHM","CD79B","CD37","CD52","MS4A1","CD27","CD38","CD40","JCHAIN",
                                         "IGHA1","IGHG4","IGHA2","IGHG3","IGHG2"))+
  theme(axis.text.x = element_text(angle=60,hjust=1))




####============================Trajectory Analysis==================================

library(monocle)


data <- as(as.matrix(obj.all.four@assays$RNA@scale.data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data =obj.all.four@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
HSMM2 <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        #lowerDetectionLimit = 0.5,
                        expressionFamily =uninormal())


ordering_genes<-(obj.all.four_marker$gene)

HSMM2 <- setOrderingFilter(HSMM2, ordering_genes)

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM2 <- reduceDimension(HSMM2, norm_method="none",
                         reduction_method="DDRTree",
                         max_components=2,
                         scaling=TRUE,
                         verbose=TRUE,
                         pseudo_expr=0)

pData(HSMM2)$celltype=obj.all.four@meta.data$customclassif
pData(HSMM2)$source=obj.all.four@meta.data$orig.ident

HSMM2 <- orderCells(HSMM2)



# .....trajectory plot

plot_cell_trajectory(HSMM2, 
                     color_by = "customclassif",
                     theta = -15,
                     show_branch_points = T,
                     show_tree = TRUE,
                     cell_size =1)

plot_cell_trajectory(HSMM2, color_by = "Pseudotime",cell_size = 1) 

plot_cell_trajectory(HSMM2, color_by = "celltype",cell_size = 1)

plot_cell_trajectory(HSMM2, color_by = "source",cell_size = 1)


# DifferentialGeneTest for pseudotime

to_be_tested <- row.names(subset(fData(HSMM2),
                                 gene_short_name %in% obj.all.four_marker$gene))
cds_subset <- HSMM2[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval <0.01)

sig_gene_names <- row.names(subset(diff_test_res, qval <0.01))


#branchHeatmap

plot_genes_branched_heatmap(HSMM2[sig_gene_names,],
                            
                            num_clusters = 5,
                            cores = 1)

t=plot_genes_branched_heatmap(HSMM2[sig_gene_names,],
                              
                              num_clusters = 5,
                              cores = 1,
                              return_heatmap = T)

cluster_1_3_5=filter((t$annotation_row), Cluster%in%c(1,3,5))

#------------Pseudotemporal Expression Pattern

plot_genes_in_pseudotime(HSMM2[c("IGKC","JCHAIN","IGHG1","IGHG2", "IGHG3","IGHG4"),],ncol=2,color_by = "orig.ident")


#===============pathway analysis GO and KEGG==============================================

library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(GO.db)
library(org.Hs.eg.db)
library(data.table)
require(DOSE)
library(DOSE)
library(clusterProfiler)
library(purrr)


#==Preparing genes for pathway analysis

SYMBOL2EG <-
  eval(parse(text = sprintf(
    'org.%s.egSYMBOL2EG', 'Hs'
  )))

clusters_genes=filter(obj.all.four_marker,gene %in%rownames(cluster_1_3_5))
logFC_score <- clusters_genes$avg_log2FC
logFC_score
de_genes<- (clusters_genes)$gene
length(x = de_genes)
names(logFC_score) <- de_genes
logFC_score

genes <- intersect(de_genes, mappedkeys(SYMBOL2EG)) #  Get the gene symbol that are mapped to an entrez gene identifiers
#length(x = genes)


logFC_score<- logFC_score[genes] # access logFC score 
#length(x = logFC_score) 
gene_entrez <-genes %>% SYMBOL2EG[.] %>% as.list %>% map( ~ .[1]) %>% purrr::simplify()
#length(x = gene_entrez)

names(logFC_score) <- gene_entrez
logFC_score # input for rank in fgsea


gse <- gseGO(geneList=sort(logFC_score, decreasing = T), 
             ont ="BP", 
             pvalueCutoff = 0.05,
             verbose = TRUE,
             scoreType = "pos",
             OrgDb = org.Hs.eg.db
)



gse=gse%>%arrange(desc(enrichmentScore))
go_result= gse@result[1:20,]
go_result$GeneRatio=go_result$enrichmentScore


library(ggplot2)

ggplot(go_result,aes(
  x=reorder(Description,GeneRatio), 
  y=GeneRatio,
  color = p.adjust, size = GeneRatio
)) +
  geom_point()+
  scale_color_gradient(low = "red", high = "blue")+coord_flip()+theme_bw()+
  
  scale_x_discrete(labels = function(x) str_wrap(x, width =70))+
  theme(axis.text=element_text(color="black",size=13),
        panel.border = element_rect(color="black",fill=NA,linewidth =0.5),
  )+
  labs(title = "",
       x="",
       y = 'GeneRatio')





# Prognosis value of B cell subtypes===============

library(readxl)

subtypes_marker = obj.all.four_marker

LIHC_mrna <- read_excel("LIHC_mrna.xlsx")

subtypes_marker=filter(subtypes_marker, gene %in% LIHC_mrna$`Updated Name`)

subtypes_marker = subtypes_marker %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

naive_gene=subtypes_marker[1:30,]
memory_gene=subtypes_marker[31:60,]
plasma_gene=subtypes_marker[61:90,]




memoryGene_coeff=filter(LIHC_mrna,`Updated Name`%in% memory_gene$gene)
naiveGene_coeff=filter(LIHC_mrna,`Updated Name`%in% naive_gene$gene)
plasmaGene_coeff=filter(LIHC_mrna,`Updated Name`%in% plasma_gene$gene)

memoryGene_coeff_v=mean(memoryGene_coeff$`Cox coefficient`)
naiveGene_coeff_v=mean(naiveGene_coeff$`Cox coefficient`)
plasmaGene_coeff_v=mean(plasmaGene_coeff$`Cox coefficient`)


memoryCoeff=memoryGene_coeff[,c("Cox coefficient","Updated Name")]
memoryCoeff$subtype= "Memory Bcell"

naiveCoeff=naiveGene_coeff[,c("Cox coefficient","Updated Name")]
naiveCoeff$subtype= "Naive Bcell"

plasmaCoeff=plasmaGene_coeff[,c("Cox coefficient","Updated Name")]
plasmaCoeff$subtype= "Plasma cell"

df= rbind(memoryCoeff,naiveCoeff,plasmaCoeff)
df$coefficient=df$`Cox coefficient`

library(ggpubr)

ggbarplot(
  df, x = "subtype", y = "coefficient",
  add = c("mean_sd"),
  fill = c("#807F7F"),ylab="Cox coefficient",xlab="",width = 0.4
)+geom_jitter(position = position_jitter(height = .1, width = .08),size=1)+ylim(-0.4,0.4)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust=0.4))+geom_hline(yintercept=0)


#======  Downloading TCGA LIHC data from TCGA database

query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query = query_TCGA)

tcga_data_lihc = GDCprepare(query_TCGA)

#====== Preparing TCGA dataframe with features===============

tcga_data_lihc$deceased <- ifelse(tcga_data_lihc$vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
tcga_data_lihc$overall_survival <- ifelse(tcga_data_lihc$vital_status == "Alive",
                                          tcga_data_lihc$days_to_last_follow_up,
                                          tcga_data_lihc$days_to_death)

tcga_data_lihc= tcga_data_lihc [,!is.na(tcga_data_lihc$overall_survival)]

#data normalization

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = assay(tcga_data_lihc),
                              colData = colData(tcga_data_lihc),
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# vst 
vsd <- vst(dds, blind=FALSE)

bulk.mtx=as.data.frame(assay(vsd)) %>%
  rownames_to_column(var = 'gene_id') %>%  
  left_join(., as.data.frame(tcga_data_lihc@rowRanges)[,c("gene_id","gene_name")], by = "gene_id") 


bulk.mtx$gene_id=NULL
rownames(bulk.mtx)=make.unique(bulk.mtx$gene_name)
bulk.mtx$gene_name=NULL
bulk.mtx_t=t(bulk.mtx)
bulk.mtx_t=rownames_to_column(as.data.frame(bulk.mtx_t),"sample_id")

clinical_data=as.data.frame(colData(tcga_data_lihc))
clinical_data=clinical_data[,c("sample_type","ajcc_pathologic_stage","gender","vital_status", "days_to_last_follow_up", "days_to_death","deceased","overall_survival")]

clinical_data$stage=clinical_data$ajcc_pathologic_stage
clinical_data$stage[which(str_detect(clinical_data$stage,"Stage III"))]="Stage III"
clinical_data$stage[which(str_detect(clinical_data$stage,"Stage IV"))]="Stage IV"
clinical_data$stage[which(is.na(clinical_data$stage))]="Normal"
clinical_data=rownames_to_column(clinical_data,"sample_id")
merge_lihc=merge(clinical_data, bulk.mtx_t,by="sample_id")



#survival curve MS4A1 gene

tumor_data=merge_lihc[(merge_lihc$sample_type == "Primary Tumor"),]

gene="MS4A1"

median_value <- median(tumor_data[, gene])

tumor_data[,gene] <- ifelse(tumor_data[,gene] >= median_value, "HIGH", "LOW")

tumor_data$Survival_years= tumor_data$overall_survival/365

fit <-survfit(as.formula(paste('Surv(Survival_years, deceased)~', gene)),data=tumor_data)

ggsurvplot(fit, data=tumor_data,surv.median.line = "none", pval=T, risk.table=F,
           xlab = "Time in years",legend.title="",legend = c(0.8, 0.9),legend.labs = 
             c("High expression", "Low expression"))$plot+ggtitle(gene)+
        theme(plot.title = element_text(hjust=0.5, face = "bold",size=13))


# Boxplot of MS4A1 gene between normal and tumor samples

ggboxplot(merge_lihc, 
          x = "sample_type", 
          y = "MS4A1", 
          fill="sample_type",
          outlier.shape = NA,
          bxp.errorbar = T,
          bxp.errorbar.width = 0.2,
          xlab=FALSE,ylab="log2(TPM + 1)",title = "MS4A1")+ 
  theme(legend.position = "none",plot.title = element_text(hjust=0.5, face = "bold",size=14))+
  geom_jitter(position = position_jitter(height = 4, width = .08),size=0.6)+
  geom_pwc(method = "t_test",label="p.format",bracket.nudge.y = 0.4, size=0.5,tip.length = 0.04)



# Boxplot for comparing stage of tumor

my_comparisons <- list( c("Stage I", "Stage II"),c("Stage I", "Stage III"), c("Stage I", "Stage IV") )


ggboxplot((tumor_data[which(!str_detect(tumor_data$stage,"Normal")),]), x = "stage", y = "MS4A1",
          fill = "stage")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)
