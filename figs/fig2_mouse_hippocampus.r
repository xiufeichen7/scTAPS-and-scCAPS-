library(scales)
library(Seurat)
library(cowplot)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(viridis)
setwd("/Users/chengjingfei/Desktop/mouse_neuron/")
mycolors <- c("#7994C6","#e01a30","#039fcc","#7d2574")
#### 4a. genome methylation ####
qc <- read.csv("stats/all_stats.info",header=TRUE)
exclude_cells <- qc[qc$proper_nmap<500000|qc$proper_nmap>3000000,1] # exclude sample by abnormal proper_nmap reads
qc$lib <- gsub("_N.*","",qc$smp)

qc$lib <- factor(qc$lib,levels=c("C211","C213","C206","C204"))
qc <- qc[order(qc$lib),]
theme_set(theme_light(base_size = 18))
qc %>%
  select(lib, chr_rC) %>%
  group_by(lib) %>%
  summarise(mean = mean(chr_rC))
# 1 C211  24.3 "NeuN+ aged"
# 2 C213  22.0 "NeuN+"
# 3 C206   9.72 "NeuN- aged"
# 4 C204   9.29 "NeuN-"

p <- qc[!qc$smp%in%exclude_cells,] %>%
  ggplot(aes(x = lib, y = chr_rC, color = lib)) +
  coord_flip() +
  labs(x = "", y = "hmC level%") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12,color="black"),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 1, alpha = 0.25, width = 0.15) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_manual(values=mycolors)+
  scale_x_discrete(breaks=c("C204","C206","C213","C211"), 
                   labels=c('NeuN-','NeuN- aged','NeuN+','NeuN+ aged'))
ggsave(p, width =5, height = 3, filename = "plots/chr_sta.pdf")
ggsave(p, width =5, height = 3, filename = "plots/chr_sta.png")
#### 4b. clustering based on gene methylation ####
pos <- read.table("resource/gencode.vM1.annotation.protein_coding.bed",header=FALSE)
pos$id <- gsub("\\..*","",pos$V4)
colnames(pos)[1:3] <- c("CHROM","START","END")

genemeth <- read.table("meth/all_sample.CpG.gene.txt",header=TRUE,sep="\t")
genemeth <- merge(pos[,c(1:3,4,5)],genemeth, by=c("CHROM","START","END"))
rownames(genemeth) <- paste0(genemeth$V5,":",genemeth$V4)
genemeth <- genemeth[,-c(1:5)]
colnames(genemeth) <- gsub("_rev_","_",colnames(genemeth))
genemeth <- genemeth[,!colnames(genemeth)%in%c(paste0(exclude_cells,"_MOD"),
                                               paste0(exclude_cells,"_TOTAL"))]
min_cov1 <- 5
gene_aC <- genemeth %>% dplyr::select(grep("_TOTAL$",colnames(genemeth), value = TRUE))
sel_genes <- rownames(gene_aC)[apply(gene_aC,1,min)>min_cov1]
length(sel_genes)


mC1<- genemeth %>% dplyr::select(grep("_MOD$",colnames(genemeth), value = TRUE))
aC1 <- genemeth %>% dplyr::select(grep("_TOTAL$",colnames(genemeth), value = TRUE))
meth1 <- mC1/aC1



bin100k <- read.table("meth/all_sample.CpG.genome_bin100k.txt",header=TRUE,sep="\t")
rownames(bin100k) <- paste0(bin100k$CHROM,":",bin100k$START)
bin100k <- bin100k[,-c(1:3)]
colnames(bin100k) <- gsub("_rev_","_",colnames(bin100k))
bin100k <- bin100k[,!colnames(bin100k)%in%c(paste0(exclude_cells,"_MOD"),
                                            paste0(exclude_cells,"_TOTAL"))]

min_cov2 <- 3
max_cov2 <- 3000
bin100k_aC <- bin100k %>% dplyr::select(grep("_TOTAL$",colnames(bin100k), value = TRUE))
sel_bin100k <- rownames(bin100k_aC)[apply(bin100k_aC,1,min)>min_cov2 & apply(bin100k_aC,1,max)<max_cov2]
length(sel_bin100k)

mC2<- bin100k %>% dplyr::select(grep("_MOD$",colnames(bin100k), value = TRUE))
aC2 <- bin100k %>% dplyr::select(grep("_TOTAL$",colnames(bin100k), value = TRUE))
meth2 <- mC2/aC2

# meth_sel <- rbind(meth1[rownames(meth1)%in%sel_genes,], meth2[rownames(meth2) %in% sel_bin100k,])
meth_sel <- meth1[rownames(meth1)%in%sel_genes,]
dim(meth_sel)
##### 1. neuron vs. non neuron #####
# mycolors <- c("#e01a30","#7d2574","#039fcc","#7994C6")
# "C204","C206","C213","C211"
mycolors <- c("#7d2574","#e01a30")
cmp <- "NeuN_vs_nonNeuN"　
obj <- CreateSeuratObject(
  meth_sel[,grep("C204|C213",colnames(meth_sel))],
  project = "CreateSeuratObject"
)
obj$orig.ident <- gsub("C204","NeuN-",obj$orig.ident)
obj$orig.ident <- gsub("C206","NeuN- aged",obj$orig.ident)
obj$orig.ident <- gsub("C213","NeuN+",obj$orig.ident)
obj$orig.ident <- gsub("C211","NeuN+ aged",obj$orig.ident)

###### 1.1 cluster ###### 
obj <- NormalizeData(obj)#, normalization.method = "RC")
obj <- FindVariableFeatures(obj)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:20, k.param = 15)
obj <- FindClusters(obj, resolution = 0.1)
obj <- RunUMAP(obj, dims = 1:20)
obj <- RunTSNE(obj, dims = 1:20)


p1 <- DimPlot(obj, reduction = "pca") + DimPlot(obj, reduction = "pca", group.by = "orig.ident")
p2 <- DimPlot(obj, reduction = "umap") + DimPlot(obj, reduction = "umap", group.by = "orig.ident")
p3 <- DimPlot(obj, reduction = "tsne") + DimPlot(obj, reduction = "tsne", group.by = "orig.ident")

p <- plot_grid(p1 & scale_colour_manual(values = mycolors),
          p2 & scale_colour_manual(values = mycolors),
          p3 & scale_colour_manual(values = mycolors),nrow=3 ) 

ggsave(paste0("plots/",cmp,"_cluster.pdf"),p,width = 10,height = 10)
ggsave(paste0("plots/",cmp,"_cluster.png"),p,width = 10,height = 10)
###### 1.2 find markers ####
cluster.markers <- FindMarkers(obj, ident.1 = c(0), ident.2 = c(1), min.pct = 0.25)
cluster.markers$name <- gsub(":.*","",rownames(cluster.markers))
write.csv(cluster.markers,paste0("plots/",cmp,"_markers.csv"),quote = FALSE)

# avg_log2FC < 0, Higher in NeuN+
# avg_log2FC > 0, Higher in NeuN-
mark1 <- gsub(":.*","", cluster.markers %>% filter(avg_log2FC<0 & p_val_adj<0.05) %>% rownames()) #NeuN+
mark2 <- gsub(":.*","", cluster.markers %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% rownames()) #NeuN-



dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023","Tabula_Muris")
library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  enriched_mark1 <- enrichr(mark1, dbs)
  enriched_mark2 <- enrichr(mark2, dbs)
}

###### 1.2.1 enrich tabula_Muris ###### 
sel_n <- 1
enrich_term <- rbind(head(enriched_mark1$Tabula_Muris,sel_n),
      head(enriched_mark2$Tabula_Muris,sel_n))
enrich_term$group <- c(rep("NeuN+",sel_n),rep("NeuN-",sel_n))
enrich_term$termshort <- gsub("\\(.*\\)","",enrich_term$Term)
enrich_term$logp <- -log(enrich_term$Adjusted.P.value,10)
p <- enrich_term %>% 
  ggplot(aes(x = termshort, y = logp,fill=group)) + 
  geom_bar(stat = "identity",width = 0.5) +
  coord_flip() +
  ylab("-log(P_value,10)") + 
  xlab("") + 
  ggtitle("Enriched ontology in Tabula_Muris") +
  scale_fill_manual(values=mycolors)
ggsave(paste0("plots/",cmp,"_top1_onotology.pdf"),p,width = 10,height = 2.5)
ggsave(paste0("plots/",cmp,"_top1_onotology.png"),p,width = 10,height = 2.5)

###### 1.2.2 enrich go ###### 
sel_n <- 10
enrich_term <- head(enriched_mark1$GO_Biological_Process_2023,sel_n)
enrich_term$termshort <- gsub("\\(.*\\)","",enrich_term$Term)
enrich_term$id <- gsub("\\)","",gsub(".*\\(GO","GO",enrich_term$Term))
enrich_term <- enrich_term[order(enrich_term$Adjusted.P.value,decreasing = TRUE),]
enrich_term$logp <- -log(enrich_term$Adjusted.P.value,10)
selterm <- enrich_term$termshort
enrich_term$termshort <- factor(enrich_term$termshort, levels=selterm)
selid<- enrich_term$id
enrich_term$id <- factor(enrich_term$id, levels=selid)
enrich_term <- enrich_term[order(enrich_term$termshort),]
enrich_term <- enrich_term[order(enrich_term$id),]


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, max(enrich_term$logp)))

p1 <- enrich_term %>%
  ggplot(aes(x = termshort, y = logp, 
             color = logp, size = Odds.Ratio)) + 
  # facet_grid(rows = vars(group))+
  geom_point() +
  coord_flip() +
  sc +
  ylab("-log(P_value,10)") + 
  xlab("") + 
  ggtitle("GO enrichment analysis") +
  ylim(0,max(enrich_term$logp))
p2 <- enrich_term %>%
  ggplot(aes(x = id, y = logp, 
             color = logp, size = Odds.Ratio)) + 
  # facet_grid(rows = vars(group))+
  geom_point() +
  coord_flip() +
  sc +
  ylab("-log(P_value,10)") + 
  xlab("") + 
  ggtitle("GO enrichment analysis") +
  ylim(0,max(enrich_term$logp))
p <- plot_grid(p1,p2,nrow=1,rel_widths = c(2,1))

ggsave(paste0("plots/",cmp,"_top100_go_bp.NeuN.pdf"),p,width = 18,height = 4)
ggsave(paste0("plots/",cmp,"_top100_go_bp.NeuN.png"),p,width = 18,height = 4)


sel_n <- 10
enrich_term <- head(enriched_mark2$GO_Biological_Process_2023,sel_n)
enrich_term$termshort <- gsub("\\(.*\\)","",enrich_term$Term)
enrich_term$id <- gsub("\\)","",gsub(".*\\(GO","GO",enrich_term$Term))
enrich_term <- enrich_term[order(enrich_term$Adjusted.P.value,decreasing = TRUE),]
enrich_term$logp <- -log(enrich_term$Adjusted.P.value,10)
selterm <- enrich_term$termshort
enrich_term$termshort <- factor(enrich_term$termshort, levels=selterm)
selid<- enrich_term$id
enrich_term$id <- factor(enrich_term$id, levels=selid)
enrich_term <- enrich_term[order(enrich_term$termshort),]
enrich_term <- enrich_term[order(enrich_term$id),]


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, max(enrich_term$logp)))

p1 <- enrich_term %>%
  ggplot(aes(x = termshort, y = logp, 
             color = logp, size = Odds.Ratio)) + 
  # facet_grid(rows = vars(group))+
  geom_point() +
  coord_flip() +
  sc +
  ylab("-log(P_value,10)") + 
  xlab("") + 
  ggtitle("GO enrichment analysis") +
  ylim(0,max(enrich_term$logp))
p2 <- enrich_term %>%
  ggplot(aes(x = id, y = logp, 
             color = logp, size = Odds.Ratio)) + 
  # facet_grid(rows = vars(group))+
  geom_point() +
  coord_flip() +
  sc +
  ylab("-log(P_value,10)") + 
  xlab("") + 
  ggtitle("GO enrichment analysis") +
  ylim(0,max(enrich_term$logp))
p <- plot_grid(p1,p2,nrow=1,rel_widths = c(2,1))
ggsave(paste0("plots/",cmp,"_top100_go_bp.non_NeuN.pdf"),p,width = 18,height = 4)
ggsave(paste0("plots/",cmp,"_top100_go_bp.non_NeuN.png"),p,width = 18,height = 4)

sapply(names(enriched_mark1), 
       function (x) write.csv(enriched_mark1[[x]], file=paste("plots/",cmp,"_", x, ".NeuN.csv", sep=""),quote = FALSE))
sapply(names(enriched_mark2), 
       function (x) write.csv(enriched_mark2[[x]], file=paste("plots/",cmp,"_", x, ".non_NeuN.csv", sep=""),quote = FALSE))

###### 1.2.3 plot marker genes #####
# Cnksr3|SEMA4D|TNR|DOCK5|MOB3B
# RBFOX3|CNTNAP2|celf4|grm1|syt1
marker_1 <- grep("Cnksr3|SEMA4D|DOCK5|MOB3B",rownames(cluster.markers),ignore.case = TRUE,value=TRUE)
p1 <- FeaturePlot(obj, features = marker_1,
                 reduction = "tsne",ncol=4, pt.size = 0.5) & scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))
marker_2 <- grep("RBFOX3|CNTNAP2|syt1:|grm1",rownames(cluster.markers),ignore.case = TRUE,value=TRUE)
p2 <- FeaturePlot(obj, features = marker_2,
                  reduction = "tsne",ncol=4, pt.size = 0.5) & scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))
p <- plot_grid(p1,p2,nrow=2)
ggsave(paste0("plots/",cmp,"_markers.png"),p,width = 15,height = 7)
ggsave(paste0("plots/",cmp,"_markers.pdf"),p,width = 15,height = 7)

##### 2. Aging neuron #####
# mycolors <- c("#e01a30","#7d2574","#039fcc","#7994C6")
# "C204","C206","C213","C211"
mycolors <- c("#e01a30","#7994C6")
cmp <- "Aging_NeuN"
obj <- CreateSeuratObject(
  meth_sel[,grep("C211|C213",colnames(meth_sel))],
  project = "CreateSeuratObject"
)
obj$orig.ident <- gsub("C204","NeuN-",obj$orig.ident)
obj$orig.ident <- gsub("C206","NeuN- aged",obj$orig.ident)
obj$orig.ident <- gsub("C213","NeuN+",obj$orig.ident)
obj$orig.ident <- gsub("C211","NeuN+ aged",obj$orig.ident)
###### 2.1 cluster ######
obj <- NormalizeData(obj)#, normalization.method = "RC")
obj <- FindVariableFeatures(obj)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:10, k.param = 15)
obj <- FindClusters(obj, resolution = 0.1)
obj <- RunUMAP(obj, dims = 1:10)
obj <- RunTSNE(obj, dims = 1:10)


p1 <- DimPlot(obj, reduction = "pca") + DimPlot(obj, reduction = "pca", group.by = "orig.ident")
p2 <- DimPlot(obj, reduction = "umap") + DimPlot(obj, reduction = "umap", group.by = "orig.ident")
p3 <- DimPlot(obj, reduction = "tsne") + DimPlot(obj, reduction = "tsne", group.by = "orig.ident")

p <- plot_grid(p1 & scale_colour_manual(values = mycolors),
               p2 & scale_colour_manual(values = mycolors),
               p3 & scale_colour_manual(values = mycolors),nrow=3 ) 

ggsave(paste0("plots/",cmp,"_cluster.pdf"),p,width = 10,height = 10)
ggsave(paste0("plots/",cmp,"_cluster.png"),p,width = 10,height = 10)

###### 2.2 find markers #####
cluster.markers <- FindMarkers(obj, ident.1 = c(0), ident.2 = c(1), min.pct = 0.25)
cluster.markers$name <- gsub(":.*","",rownames(cluster.markers))
write.csv(cluster.markers,paste0("plots/",cmp,"_markers.csv"),quote = FALSE)

###### 2.2.1 relate with age markers #####
age_gene <- read.table("plots/age_correlated_gene_2023_cell.txt",header=TRUE,sep="\t")

# avg_log2FC < 0, Higher in Aged
# avg_log2FC > 0, Higher in Young
mark1 <- gsub(":.*","", cluster.markers %>% filter(avg_log2FC<0)  %>% rownames()) #NeuN+
mark2 <- gsub(":.*","", cluster.markers %>% filter(avg_log2FC>0) %>% rownames()) #NeuN-
age_gene[age_gene$Gene.Symbol%in%c(mark1),] %>%
  ggplot(aes(x="up",y=Spearman.s.correlation.coefficient))+
  geom_violin() 
age_gene[age_gene$Gene.Symbol%in%c(mark2),] %>%
  ggplot(aes(x="down",y=Spearman.s.correlation.coefficient))+
  geom_violin()
dat_m <- merge(cluster.markers,age_gene,by.x=("name"),by.y=("Gene.Symbol")) %>%
  distinct(name, .keep_all = TRUE)
p <- ggplot(dat_m, aes(x=Spearman.s.correlation.coefficient,y=-avg_log2FC))+
  geom_point() +
  xlab("Spearman's correlation \n between gene expression and ageing") +
  ylab("hmC log2FC (Aged / Young) ") +
  geom_smooth() +
  annotate(geom="text", x=-0.4, y=0.4, label="cor: 0.079 \n P-value: 0.6331",color="blue")
p

cor.test(dat_m$Spearman.s.correlation.coefficient,-dat_m$avg_log2FC)

# Pearson's product-moment correlation
# 
# data:  dat_m$Spearman.s.correlation.coefficient and -dat_m$avg_log2FC
# t = 0.48131, df = 37, p-value = 0.6331
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2426772  0.3848193
# sample estimates:
#        cor 
# 0.07887955 

ggsave(paste0("plots/",cmp,"_correlation_with_ageing_gene_expr.pdf"),p,width = 6,height = 6)
ggsave(paste0("plots/",cmp,"_correlation_with_ageing_gene_expr.png"),p,width = 6,height = 6)


marker_1 <- grep("Trim2|grip1|oxr1|hcn1|aldh1l1|foxo3|megf10|camk4|tox2|shank2|runx1|app",rownames(cluster.markers),ignore.case = TRUE,value=TRUE)
p1 <- FeaturePlot(obj, features = marker_1,
                  reduction = "tsne",ncol=4, pt.size = 0.5) & scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))
ggsave(paste0("plots/",cmp,"_markers.png"),p1,width = 15,height = 10)
ggsave(paste0("plots/",cmp,"_markers.pdf"),p1,width = 15,height = 10)

###### 2.2.2 markers #####
marker_1 <- grep("Tcf4|App|Syt1|kcnip4|cdon|cntn4|prkca|nrg2|hdac9|usp7",rownames(cluster.markers),ignore.case = TRUE,value=TRUE)
p1 <- FeaturePlot(obj, features = marker_1,
                  reduction = "tsne",ncol=4, pt.size = 0.5) & scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))
ggsave(paste0("plots/",cmp,"_markers.png"),p1,width = 15,height = 7.5)
ggsave(paste0("plots/",cmp,"_markers.pdf"),p1,width = 15,height = 7.5)


##### 3. Aging non neuron #####
# mycolors <- c("#e01a30","#7d2574","#039fcc","#7994C6")
# "C204","C206","C213","C211"
mycolors <- c("#7d2574","#039fcc")
cmp <- "Aging_non_NeuN"
obj <- CreateSeuratObject(
  meth_sel[,grep("C204|C206",colnames(meth_sel))],
  project = "CreateSeuratObject"
)
obj$orig.ident <- gsub("C204","NeuN-",obj$orig.ident)
obj$orig.ident <- gsub("C206","NeuN- aged",obj$orig.ident)
obj$orig.ident <- gsub("C213","NeuN+",obj$orig.ident)
obj$orig.ident <- gsub("C211","NeuN+ aged",obj$orig.ident)

obj <- NormalizeData(obj)#, normalization.method = "RC")
obj <- FindVariableFeatures(obj)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:10, k.param = 15)
obj <- FindClusters(obj, resolution = 0.1)
obj <- RunUMAP(obj, dims = 1:10)
obj <- RunTSNE(obj, dims = 1:10)


p1 <- DimPlot(obj, reduction = "pca") + DimPlot(obj, reduction = "pca", group.by = "orig.ident")
p2 <- DimPlot(obj, reduction = "umap") + DimPlot(obj, reduction = "umap", group.by = "orig.ident")
p3 <- DimPlot(obj, reduction = "tsne") + DimPlot(obj, reduction = "tsne", group.by = "orig.ident")

p <- plot_grid(p1 & scale_colour_manual(values = mycolors),
               p2 & scale_colour_manual(values = mycolors),
               p3 & scale_colour_manual(values = mycolors),nrow=3 ) 

ggsave(paste0("plots/",cmp,"_cluster.pdf"),p,width = 10,height = 10)
ggsave(paste0("plots/",cmp,"_cluster.png"),p,width = 10,height = 10)

###### 3.2 find markers #####
cluster.markers <- FindMarkers(obj, ident.1 = c(0), ident.2 = c(1), min.pct = 0.25)
cluster.markers$name <- gsub(":.*","",rownames(cluster.markers))
write.csv(cluster.markers,paste0("plots/",cmp,"_markers.csv"),quote = FALSE)


###### 3.2.1 relate with age markers #####
age_gene <- read.table("plots/age_correlated_gene_2023_cell.txt",header=TRUE,sep="\t")

# avg_log2FC < 0, Higher in Aged
# avg_log2FC > 0, Higher in Young
mark1 <- gsub(":.*","", cluster.markers %>% filter(avg_log2FC<0)  %>% rownames()) #NeuN+
mark2 <- gsub(":.*","", cluster.markers %>% filter(avg_log2FC>0) %>% rownames()) #NeuN-
age_gene[age_gene$Gene.Symbol%in%c(mark1),] %>%
  ggplot(aes(x="up",y=Spearman.s.correlation.coefficient))+
  geom_violin() 
age_gene[age_gene$Gene.Symbol%in%c(mark2),] %>%
  ggplot(aes(x="down",y=Spearman.s.correlation.coefficient))+
  geom_violin()
dat_m <- merge(cluster.markers,age_gene,by.x=("name"),by.y=("Gene.Symbol")) %>%
  distinct(name, .keep_all = TRUE)
p <- ggplot(dat_m, aes(x=Spearman.s.correlation.coefficient,y=-avg_log2FC))+
  geom_point() +
  xlab("Spearman's correlation \n between gene expression and ageing") +
  ylab("hmC log2FC (Aged / Young) ") +
  geom_smooth() +
  annotate(geom="text", x=-0.4, y=0.4, label="cor: 0.334 \n P-value: 2.7E-05",color="blue")
p

cor.test(dat_m$Spearman.s.correlation.coefficient,-dat_m$avg_log2FC)

# Pearson's product-moment correlation
# 
# data:  dat_m$Spearman.s.correlation.coefficient and -dat_m$avg_log2FC
# t = 4.3311, df = 149, p-value = 2.713e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1845156 0.4690639
# sample estimates:
#       cor 
# 0.3343888

ggsave(paste0("plots/",cmp,"_correlation_with_ageing_gene_expr.pdf"),p,width = 6,height = 6)
ggsave(paste0("plots/",cmp,"_correlation_with_ageing_gene_expr.png"),p,width = 6,height = 6)
p <- ggplot(dat_m, aes(x=Spearman.s.correlation.coefficient,y=-avg_log2FC,label=name))+
  geom_point() +
  xlab("Spearman's correlation \n between gene expression and ageing") +
  ylab("hmC log2FC (Aged / Young) ") +
  geom_smooth() +
  annotate(geom="text", x=-0.4, y=0.4, label="cor: 0.334 \n P-value: 2.7E-05",color="blue")+
geom_text(hjust=0, vjust=0)
ggsave(paste0("plots/",cmp,"_correlation_with_ageing_gene_expr.label.pdf"),p,width = 10,height = 10)
ggsave(paste0("plots/",cmp,"_correlation_with_ageing_gene_expr.label.png"),p,width = 10,height = 10)

###### 3.2.2 markers #####
marker_1 <- grep("Trim2|grip1|oxr1|hcn1|foxo3|megf10|camk4|tox2|shank2|runx1|app|Args|Edil3|Epha3",rownames(cluster.markers),ignore.case = TRUE,value=TRUE)
p1 <- FeaturePlot(obj, features = marker_1,
                  reduction = "tsne",ncol=4, pt.size = 0.5) & scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))
ggsave(paste0("plots/",cmp,"_markers.png"),p1,width = 15,height = 10)
ggsave(paste0("plots/",cmp,"_markers.pdf"),p1,width = 15,height = 10)


# ###### 3.3.3 Plot marker gene number ######
# aging_neuron_marker <- read.csv("plots/Aging_NeuN_markers.csv")
# aging_non_neuron_marker <- read.csv("plots/Aging_non_NeuN_markers.csv")
# aging_neuron_marker[grep("Syt1:|Kcnip4",aging_neuron_marker$X),]
# 
# dmg <- rbind(nrow(aging_neuron_marker %>% filter(p_val_adj<0.05)),
#              nrow(aging_non_neuron_marker %>% filter(p_val_adj<0.05))) %>%
#   as.data.frame()
# rownames(dmg) <- c("neuron","non-neuron")
# ggplot(dmg,aes(x=rownames(dmg),y=V1)) + 
#   geom_bar(stat="identity") +
#   coord_flip() +
#   ylab("Marker genes ") + 
#   xlab("") + 
#   scale_fill_manual(values=mycolors)
# ggsave(paste0("plots/",cmp,"_marker_num.pdf"),p,width = 5,height = 2)
# ggsave(paste0("plots/",cmp,"_marker_num.png"),p,width = 5,height = 2)
