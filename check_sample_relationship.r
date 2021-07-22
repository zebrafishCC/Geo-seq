###data analysis for H135 Geo-seq data

#contents
#1. basic statistics 
   #gene number
#heatmap all samples before filtering and after filtering
#plot PCA plots according to PCA plots

#6. sample name transversion

#2. normalization

#3. density plot and boxplot

#4. PCA analysis and heatmap clustering

#5. DEGs



#7. optional strip specific genes



###------------------------------------------------

##heatmap of all samples
setwd("/home/chengchen/data/geo-seq")
library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(reshape2)
library(matrixStats)
library(gplots)
source("heatmap_funcV2_mod.R")

# import H135 count data from local file

dat = read.delim("demo.combined.featurecount.txt",row.names = 1)
dat[1:10,1:5]
dim(dat)


##name transform. An annotation file is required.
colname = as.numeric(str_split_fixed(colnames(dat),"\\.",2)[,2])
dat = dat[,order(colname)]
sample = read.csv("h135_sampleInfo_20200511.csv")
sample$ID = gsub("-",".", sample$sampleId)
colnames(dat) = sample$sampleName[match(colnames(dat),sample$ID)]
#dat[1:10,1:5]
#dim(dat)

##datatransform important step and never forget
#dat = log2(dat+1)



##sort sample names according to their original levels(slice position)
sample_levels = function(dat) {
  colnames_classification = gsub("HE13[.]\\d([^.]+)[_]\\d+(-\\d+)?", "\\1",colnames(dat))
  table(colnames_classification)
  colnames_classification = factor(colnames_classification)
  sample_levels = levels(colnames_classification)
  sample_order = gsub("[X|Y](\\d+)","\\1",sample_levels)
  sample_levels = sample_levels[order(as.numeric(sample_order))]
  return(sample_levels)
}


##median plot for top genes 
#usage:median_plot = function(topn, raw)
#median_plot(10000,dat)

topn_raw <- function(topn, raw) {
  raw = log2(raw+1)
  raw = as.matrix(raw)
  gene_Var <-  order(rowVars(raw),decreasing=TRUE)
  raw <-  raw[gene_Var,]
  topn <- as.numeric(topn) 
  topnRaw <- raw[1:topn, ]
  print(paste("get topn ", as.character(topn) ,"raw genes"))
  return(topnRaw)
}


#plot median,not neccessary
median_plot = function(topn, raw) {
  dat= topn_raw(topn,as.matrix(dat))
  all_median = apply(dat1,2,median)
  all_median = melt(all_median)
  len = nrow(all_median)
  ord = c(1:1:len)
  all_median$ord = ord
  pdf("top10000median.pdf")
  plot = ggplot(all_median,aes(x=ord,y=value))+geom_point()
  print(plot)
  dev.off()
}


##density plot is required
#x = grepl("X80", colnames(dat))
#names1 = grep("X80", dat1$names,value = True)
##plot samples according to their levels(slice position)
plot_samples = function(dat,sample) {
  #dat = dat[,grepl(sample, colnames(dat))]
  dat = log2(dat+1)
  long = melt(dat)
  title=paste(sample,"log2Raw",sep=" ")
  print(title)
  plotfile = paste0(sample,"_density.pdf")
  print(plotfile)
  p = ggplot(data=long, aes (x=value,color=variable)) +
    geom_density(alpha=0.1) +
    theme_classic()+labs(x="log2-Rawcounts")+
    theme(legend.position = "bottom",legend.text=element_text(size=4))+
    guides(col=guide_legend(ncol=8, bycol=TRUE,title=""))+
    xlim(-1, 10)+ggtitle(sample)+theme(plot.title = element_text(hjust = 0.5))
    #ylim(0,5) + 
    #scale_x_log10() + 
    #geom_vline(xintercept = 1)+
    #theme(text = element_text(size = font.size))
  ggsave(filename=plotfile,plot=p)
}

plot_all_samples = plot_samples(dat,"allsamples")
#dat1 = dat[,grepl("X80", colnames(dat))]
#plot_all_samples = plot_samples(dat1,"X80")

for (sample in sample_levels(dat)) {
  dat1 = dat[,grepl(sample, colnames(dat))]
  plot_samples(dat1,sample)
}


##heatmap plot
#for (sample in sample_levels(dat)) {
  #dat = dat[,grepl(sample, colnames(dat))]
  #dat = topn_raw(10000,dat)
  #filename = paste0(sample,"_10000genes.pdf")
  #pdf(filename)
  #plot = heatmap_sf(dat)
  #print(plot)
  #dev.off()
#}


###############################
##plot sample-samples
methods <- c("spearman", "pearson")
for (method in methods ) {
  surfix=paste0(10000, "genes") 
  #2. select topnfpkm genes
  dat1 = topn_raw(10000,dat)
  #==============================plot heatmap 1==============
  #######1.pearson 1-cor ss
  title=paste(surfix, method, "1-cor", sep=" ")
  print(title)
  plotfile = paste0(paste(surfix, method, "1-cor", sep="."),".ss.pdf")
  print(plotfile)
  
  cormat <- cal_cor(dat, method)
  ###2.distance, can adjust the plot such as fontsize and others
  clustering_distance_cols_p <- cal_cluster_distance_col(cormat, "1-cor")
  
  ###3.plot heatmap  ss pearson 1-cor 
  #the rownames of annotation must match the names of cormat
  #ann2 <- ann[,c(3,4,6)]
  p<- pheatmap(cormat,  main=title, clustering_distance_rows=clustering_distance_cols_p, 
               clustering_distance_cols=clustering_distance_cols_p, 
               clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=F, 
               show_colnames=F,  fontsize=12) #annotation=ann2,scale='row')
  
  ggsave(filename=plotfile, plot=p)
  #plts <- c(plts, list(p))
}


#plot sample-features
methods <- c("spearman", "pearson")
for (method in methods ) {
  surfix=paste0(10000, "genes") 
  #2. select topnfpkm genes
  dat1 = topn_raw(10000,dat)
  #==============================plot heatmap 1==============
  #######1.pearson 1-cor ss
  title=paste(surfix, method, "1-cor", sep=" ")
  print(title)
  plotfile = paste0(paste(surfix, method, "1-cor", sep="."),".sf.pdf")
  print(plotfile)
  
  cormat <- cal_cor(dat, method)
  ###2.distance, can adjust the plot such as fontsize and others
  clustering_distance_cols_p <- cal_cluster_distance_col(cormat, "1-cor")
  
  ###3.plot heatmap  ss pearson 1-cor 
  #the rownames of annotation must match the names of cormat
  #ann2 <- ann[,c(3,4,6)]
  p<- pheatmap(cormat,  main=title,
               clustering_distance_cols=clustering_distance_cols_p, 
               clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=F, 
               show_colnames=F,  fontsize=12) #annotation=ann2,scale='row')
  
  ggsave(filename=plotfile, plot=p)
  #plts <- c(plts, list(p))
}


##heatmap correlation plot for specific slices
samples = sample_levels(dat)
for (sample in samples) {
  suffix=paste0(10000, "genes") 
  dat1 = dat[,grepl(sample, colnames(dat))]
  #dat1 = dat[,grepl(sample, colnames(dat))]
  dat1 = topn_raw(10000,dat1)
  #==============================plot heatmap 1==============
  #######1.pearson 1-cor ss
  title = paste0(sample,"_",suffix,"_pearson", "_cor", ".pdf")
  print(title)
  plotfile = paste0(sample,"_",suffix,"_pearson", "_cor", ".pdf")
  print(plotfile)
  
  cormat <- cal_cor(dat1, "pearson")
  ###2.distance, can adjust the plot such as fontsize and others
  #clustering_distance_cols_p <- cal_cluster_distance_col(cormat, "1-cor")
  ###3.plot heatmap  ss pearson 1-cor 
  #the rownames of annotation must match the names of cormat
  #ann2 <- ann[,c(3,4,6)]
  p<- pheatmap(cormat,  main=title,cluster_rows=F, cluster_cols=F, show_rownames=T, 
               show_colnames=F,  fontsize=12) #annotation=ann2,scale='row')
  
  ggsave(filename=plotfile, plot=p)
}


##filter reads
samples_dropped = c("HE13.5X80_28","HE13.5Y81_53","HE13.5X320-393")
dat1 = dat[,!(names(dat) %in% samples_dropped)]




####################################
#filter samples using nFeatures and nCounts
##number of features
features = as.data.frame(apply(dat, 2, function(x) count(x>0)))
colnames(features) = "nfeatures"
features$samples = rownames(features)
slice = gsub("HE13[.]\\d([^.]+)[_]\\d+(-\\d+)?", "\\1",rownames(features))
slice = factor(slice)
sample_levels = levels(slice)
sample_order = gsub("[X|Y](\\d+)","\\1",sample_levels)
sample_levels = sample_levels[order(as.numeric(sample_order))]
slice = factor(slice,levels = sample_levels)
features$slice = slice

violinplot = ggplot(features,aes(slice,nfeatures))+
  geom_violin(aes(fill = factor(slice)))+
  geom_jitter(height = 0,width = 0.01)+
  geom_hline(yintercept = 20000,linetype="dotted")+
  geom_text(data=subset(features, nfeatures<20000),aes(label=samples),
            position = position_nudge(y=-1000),  size=3)
ggsave(filename = "H13.5-nfeature.pdf",plot = violinplot)

### ncounts for each sample
counts = as.data.frame(apply(dat, 2, sum))
colnames(counts) = "ncounts"
counts$samples = rownames(counts)
counts$slice = slice

violinplot1 = ggplot(counts,aes(slice,ncounts))+
  geom_violin(aes(fill = factor(slice)))+
  geom_jitter(height = 0,width = 0.01)+
  geom_hline(yintercept = 800000,linetype="dotted")+
  geom_text(data=subset(counts, ncounts<800000),aes(label=samples),
            position = position_nudge(y=-100000),  size=3)
ggsave(filename = "H13.5-ncounts.pdf",plot = violinplot1)

##plot two variables
x = features$nfeatures
y = counts$ncounts
#features_counts = data.frame(samples=counts$samples,features=features$nfeatures,counts=counts$ncounts)
pdf("nfeatures_ncounts.pdf")
p=plot(x,y,col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3 , xlab="" , ylab="") 
print(p)
dev.off()
#model <- lm(y ~ x + I(x^2) + I(x^3))
#myPredict <- predict( model , interval="predict" )
#ix <- sort(x,index.return=T)$ix
#lines(x[ix], myPredict[ix , 1], col=2, lwd=2 )
#polygon(c(rev(x[ix]), x[ix]), c(rev(myPredict[ ix,3]), myPredict[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)


########################################
#filter 3 samples with low nfeatures
samples_dropped = c("HE13.5X80_28","HE13.5Y81_53","HE13.5X320_393")
dat1 = dat[,!(names(dat) %in% samples_dropped)]
methods <- c("spearman", "pearson")
for (method in methods ) {
  surfix=paste0("filtered3samples_",10000, "genes") 
  #2. select topnfpkm genes
  dat2 = topn_raw(10000,dat1)
  #==============================plot heatmap 1==============
  #######1.pearson 1-cor ss
  title=paste(surfix, method, "1-cor", sep=" ")
  print(title)
  plotfile = paste0(paste(surfix, method, "1-cor", sep="."),".ss.pdf")
  print(plotfile)
  
  cormat <- cal_cor(dat2, method)
  ###2.distance, can adjust the plot such as fontsize and others
  clustering_distance_cols_p <- cal_cluster_distance_col(cormat, "1-cor")
  
  ###3.plot heatmap  ss pearson 1-cor 
  #the rownames of annotation must match the names of cormat
  #ann2 <- ann[,c(3,4,6)]
  p<- pheatmap(cormat,  main=title, clustering_distance_rows=clustering_distance_cols_p, 
               clustering_distance_cols=clustering_distance_cols_p, 
               clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=F, 
               show_colnames=F,  fontsize=12) #annotation=ann2,scale='row')
  
  ggsave(filename=plotfile, plot=p)
  #plts <- c(plts, list(p))
}


#all needs to rerun tomorrow
#============================================================
#seurat analysis
features_dropped = c("HE13.5X80_28","HE13.5Y81_53","HE13.5X320_393")
feature1 = features[!(rownames(features) %in% features_dropped),]
geoseq = CreateSeuratObject(counts = dat1, project = "geoseq",meta.data = feature1)
#add metadata

#normalize,find highest variables, scale data
geoseq = NormalizeData(geoseq,normalization.method = "LogNormalize",scale.factor = 10000)
geoseq = FindVariableFeatures(geoseq, selection.method = "vst", nfeatures = 2000)
top10 = head(VariableFeatures(geoseq), 10)
all.genes = rownames(geoseq)
geoseq = ScaleData(geoseq, features = all.genes)

##run pca
geoseq = RunPCA(geoseq,features = VariableFeatures(object = geoseq))
pdf("allsamples.pdf")
p1 = DimPlot(geoseq,reduction = "pca",group.by = 'slice')
print(p1)
dev.off()

pdf("allsamples_split.pdf",height = 7,width = 21)
p1 = DimPlot(geoseq,reduction = "pca",split.by = 'slice')
print(p1)
dev.off()
#x80 = sum(grepl("X80", colnames(dat1)))
#col=heat.colors(x80)
#sub_cells <- WhichCells(geoseq, slice = "X80")
#sub_geoseq <- subset(geoseq, cells = sub_cells)
#DimPlot(sub_geoseq,reduction = "pca")


###each subset pca
pca_plot = function(slice) {
  dat3 = dat1[,grepl(slice, colnames(dat1))]
  colname3 = gsub("HE13[.]\\d[^.]+[_](\\d+-\\d+)?", "\\1",colnames(dat3))
  colnames(dat3) = colname3
  
  geoseq1 = CreateSeuratObject(counts = dat3, project = "geoseq1")
  geoseq1@meta.data$orig.ident = rownames(geoseq1@meta.data)
  #add metadata
  
  #normalize,find highest variables, scale data
  geoseq1 = NormalizeData(geoseq1,normalization.method = "LogNormalize",scale.factor = 10000)
  geoseq1 = FindVariableFeatures(geoseq1, selection.method = "vst", nfeatures = 2000)
  top10 = head(VariableFeatures(geoseq1), 10)
  all.genes = rownames(geoseq1)
  geoseq1 = ScaleData(geoseq1, features = all.genes)
  
  ##run pca
  geoseq1 = RunPCA(geoseq1,features = VariableFeatures(object = geoseq1),npcs = 20)
  #pdf("allsamples.pdf")
  Idents(geoseq1) = "orig.ident"
  filename = paste0(slice,"_pca.pdf")
  pdf(filename)
  p = DimPlot(geoseq1,reduction = "pca")
  print(p)
  dev.off()
}

slices = sample_levels(dat1)
for (slice in slices) {
  pca_plot(slice)
  print(slice)
} 
