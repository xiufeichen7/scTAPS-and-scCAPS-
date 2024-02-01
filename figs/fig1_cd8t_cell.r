library(scales)
library(Seurat)
library(cowplot)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(data.table)
library(tidyr)
library(corrplot)
setwd("/Users/chengjingfei/Desktop/fig1/")
#### 1. plot QC ####
qc <- read.csv("all_stats.new.info",header=TRUE)
exclude_cells <- qc[qc$proper_nmap<500000|qc$proper_nmap>3000000,1] # exclude sample by abnormal proper_nmap reads
qc$lib <- gsub("_N.*","",qc$smp)

qc$pct_q10_nmap <- qc$q10_nmap/qc$nclean
qc$pct_chr_nC <- qc$chr_nC/24314257
qc$per_base <- qc$chr_aC/qc$chr_nC
theme_set(theme_light(base_size = 18))

p1 <- qc %>%
  ggplot(aes(x = lib, y = pct_q10_nmap*100, fill = lib)) +
  geom_violin() +
  geom_jitter(size=0.5, width = 0.1) +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=2), size=5)+
  theme_bw() +
  ylab("q10 mapping rate%") +
  xlab("lib")+
  theme(legend.position="None") +scale_fill_manual(values=c("#A3DEF4"))+scale_y_continuous(limits=c(0,100),expand = c(0, 0))
p2 <- qc %>%
  ggplot(aes(x = lib, y = pct_chr_nC*100, fill = lib)) +
  geom_violin() +
  geom_jitter(size=0.5, width = 0.1) +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=8), size=5)+
  theme_bw() +
  ylab("genomic CpG coverage%") +
  xlab("lib")+
  theme(legend.position="None") +scale_fill_manual(values=c("#A3DEF4"))+scale_y_continuous(limits=c(0,25),expand = c(0, 0))

p3 <- qc %>%
  ggplot(aes(x=lib,y=chr_rC,fill=lib))+ geom_violin() +
  geom_jitter(size=0.5, width = 0.1) +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=0.1), size=3.5)+
  theme_bw() +
  ylab("genome methylation") +
  xlab("lib")+
  theme(legend.position="None") +
  scale_fill_manual(values=c("#A3DEF4"))+scale_y_continuous(limits=c(0,100),expand = c(0, 0))
p <- plot_grid(p1,p2,p3, labels = c('A', 'B','C'), ncol=3,label_size = 8)
ggsave(p, width =6, height = 4,filename = "sctaps_caps_mapping_taps.pdf")


p4 <- qc %>%
  ggplot(aes(x=nclean/1000000, y=as.numeric(chr_nC)/1000000,fill=lib)) +
  geom_smooth(method = "lm", se = FALSE,color="grey40")+
  geom_point(aes(size=1, color=lib))+theme_bw()+scale_color_manual(values=c("#A3DEF4"))+
  scale_y_continuous(limits=c(0,4.5),expand = c(0, 0))+
  scale_x_continuous(limits=c(0,10.5),expand = c(0, 0))+theme(legend.position="None") 
ggsave("sctaps_caps_perbase_cov_sctaps.pdf",p4, width=4,height=4)

conversion <- data.frame(value=c(sum(qc$lambda_mC)/sum(qc$lambda_aC),
                              sum(qc$hmC_mC)/sum(qc$hmC_aC),
                              sum(qc$unmeth2kb_mC)/sum(qc$unmeth2kb_aC)),
                            type=c("mC","hmC","uC"),
                         lib="scTAPS")
pdf("sctaps_caps_conversion_taps.pdf", width=4,height=4)
conversion %>%
  ggplot(aes(x=type,y=value*100,fill=lib))+ geom_bar(stat="identity") +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=2), size=3.5) +
  theme_bw() +
  ylab("modification level%") +
  xlab("lib")+
  theme(legend.position="None") +
  scale_fill_manual(values=c("#A3DEF4")) +scale_y_continuous(limits=c(-1,100),expand = c(0, 0))
dev.off()


#### 2. plot correlation in 100k ####
get_upper_tri <- function(CorMat){
  CorMat[upper.tri(CorMat)]<- NA
  return(CorMat)
}

qc[order(qc$nclean,decreasing = TRUE),] %>%head()

cov_min <- 50
cov_max <- 30000
dat <- data.frame()
for(i in list.files(pattern="[c,C]*bin*")){
  temp <- fread(i) %>% as.data.frame()
  temp$ratio <- temp[,4]/temp[,5]
  colnames(temp)[1:3] <- c("chr","start","end")
  temp <- temp[temp[,5] <cov_max & temp[,5]>cov_min,]
  temp$smp <- i
  dat <- rbind(dat,temp[,c(1:3,6:7)])
}
dat_w <- dat %>%
  pivot_wider(
  names_from = "smp",
  values_from = c("ratio"),
  names_sep = "_"
)
dat_w <- dat_w[,c(1,2,3,8,4,5,6,7)]
dat_w <- dat_w[complete.cases(dat_w),]
nrow(dat_w)/29343
cor_caps1 <- get_upper_tri(cor(dat_w[,-c(1:3)]))
meltNum <- melt(cor_caps1, na.rm = T)
meltNum$value <- as.numeric(meltNum$value)
p <- ggplot(meltNum, aes(x = Var1, y = Var2, fill = value)) + geom_tile()+coord_fixed()+
  geom_text(aes(label = round(as.numeric(value),2)), color = "black")+  
  scale_fill_gradient(low = "white", high = "#B2182B", limit = c(0,1), name = "Pearson\nCorrelation") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.background = element_blank(), panel.border = element_blank())
ggsave("heatmap.correlation_sc_taps.pdf", p,width=10, height=10)
ggsave("heatmap.correlation_sc_taps.png", p,width=10, height=10)



#### 3. heatmap in  selected regions ####
# code in cluster
library(scales)
library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(tidyr)
library(pheatmap)
smps <- gsub(".md.filter.meth.sta.txt.gz","",list.files(path = "align/", pattern="*rev*.md.filter.meth.sta.txt.gz"))
res <- data.frame(row.names = c("chr","idx","mC","aC","ratio","rC","smp"))
window_size <- 100000
for(smp in smps){
  input_file1 <- paste0("align/",smp,".md.filter.meth.sta.txt.gz")
  lowinput <- fread(paste0('zcat ', input_file1), header=TRUE)
  lowinput$ratio <- lowinput$mC/lowinput$aC
  lowinput <- lowinput[lowinput$chr%in%paste0("chr",seq(1,22)),]
  lowinput$idx <- round(lowinput$pos/window_size)*window_size
  lowinput_a <- lowinput %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC), ratio=mean(ratio)) 
  lowinput_a$rC <- lowinput_a$mC/lowinput_a$aC
  lowinput_a$smp <- smp
  res <- rbind(res,lowinput_a)
}
res_wide <- res %>%
  select(chr,idx,rC,smp) %>%
  pivot_wider(values_from = rC, names_from = smp)
res_wide %>%
  as.data.frame() %>%
  filter(chr == "chr2" & idx > 80000000 & idx < 87000000) %>%
  select(!contains("chr") & !contains("idx")) %>%
  t() %>%
  pheatmap::pheatmap(cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
                     color=colorRampPalette(c("navy", "white", "red"))(50),
                     filename = "plots/chr2_80000000_87000000.png", width = 24,height = 12,show_rownames = FALSE)
res_wide %>%
  as.data.frame() %>%
  filter(chr == "chr2" & idx > 80000000 & idx < 87000000) %>%
  select(!contains("chr") & !contains("idx")) %>%
  t() %>%
  pheatmap::pheatmap(cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
                     color=colorRampPalette(c("navy", "white", "red"))(50),
                     filename = "plots/chr2_80000000_87000000.1.pdf", width = 24,height = 12,show_rownames = FALSE)

res_wide %>%
  as.data.frame() %>%
  filter(chr == "chr2" & idx > 80000000 & idx < 87000000) %>%
  select(!contains("chr") & !contains("idx")) %>%
  t() %>%
  pheatmap::pheatmap(cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
                     filename = "plots/chr2_80000000_87000000.2.pdf", width = 24,height = 12,show_rownames = FALSE)


#### 4. methylation in chrhmm ####
dat <- read.table("all_state.meth.txt",col.names = c("smp","V1","mC","aC","rC"))
dat$V1 <- gsub("GapArtf", "assembly gaps and alignment artifacts", dat$V1)
dat$V1 <- gsub("HET", "heterochromatin", dat$V1)
dat$V1 <- gsub("Quies", "quiescent states", dat$V1)
dat$V1 <- gsub("Acet", "acetylations marks", dat$V1)
dat$V1 <- gsub("znf", "ZNF genes", dat$V1)
dat$V1 <- gsub("ReprPC", "polycomb repressed", dat$V1)
dat$V1 <- gsub("TSS", "TSS", dat$V1)
dat$V1 <- gsub("BivProm", "bivalent promoters", dat$V1)
dat$V1 <- gsub("PromF", "flanking promoters", dat$V1)
dat$V1 <- gsub("DNase", "only DNase I hypersensitivity", dat$V1)
dat$V1 <- gsub("EnhA", "strong enhancers", dat$V1)
dat$V1 <- gsub("EnhWk", "weak enhancers", dat$V1)
dat$V1 <- gsub("TxEnh", "transcribed candidate enhancers", dat$V1)
dat$V1 <- gsub("TxEx", "exon", dat$V1)
dat$V1 <- gsub("TxWk", "weak transcription", dat$V1)
dat$V1 <- gsub("Tx", "transcription", dat$V1)
levels <- c("TSS", "flanking promoters","bivalent promoters","only DNase I hypersensitivity","polycomb repressed", 
            "strong enhancers", "weak enhancers",  "transcribed candidate enhancers",
            "assembly gaps and alignment artifacts", "heterochromatin", "quiescent states", 
            "acetylations marks", "ZNF genes",
            "transcription","exon","weak transcription")
dat$V1 <- factor(dat$V1, levels=levels)
p <- dat[!dat$V1%in%c("assembly gaps and alignment artifacts","ZNF genes"),] %>% 
  ggplot(aes(x =V1 , y = as.numeric(rC)*100, fill = V1)) + 
  geom_violin() + 
  geom_jitter(size=0.5, width = 0.1) +
  theme_light() +   
  stat_summary(geom="text", fun=mean,
                                 aes(label=sprintf("%1.2f", ..y..)),
                                 position=position_nudge(y=8), size=5)+
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none")  +
  scale_y_continuous(limits=c(0,100),expand = c(0, 0))+
  #scale_fill_manual(values = c("#2D6A4E","#c7e1a6")) + 
  xlab("") + ylab("level")

ggsave("cd8_t_cell_meth_chrhmm.pdf",p,width = 10,height = 6)

