#heatmap_func.R
## Author: Fangfang Qu
##this script calculate correlation of samples, 
#plot heatmap plot for sample-sample correlation and sample-genes

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)


#zscore is same as scale(x,center = T,scale = T) 
#use zscore after normalize

#select hvgenes after log-transformation ONLY FOR FPKM NORMALIZED DATA
topn_fpkm <- function(topn, fpkm, log10=T) {
	if (log10) {
	    logExpMat <- log10(fpkm +1)
	} else {
        logExpMat <- log2(fpkm +1)
	}
	gene_Var <-  order(rowVars(logExpMat),decreasing=TRUE)
	fpkm <-  fpkm[gene_Var,]
	topn <- as.numeric(topn) 
	topnfpkm <- fpkm[1:topn, ]
	print(paste("get topn ", as.character(topn) ,"gene fpkm"))
    return(topnfpkm)
}

topn_raw <- function(topn, raw) {
  raw = log2(dat+1)
  raw = as.matrix(raw)
  gene_Var <-  order(rowVars(raw),decreasing=TRUE)
  raw <-  raw[gene_Var,]
  topn <- as.numeric(topn) 
  topnRaw <- raw[1:topn, ]
  print(paste("get topn ", as.character(topn) ,"raw genes"))
  return(topnRaw)
}

#plot median
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


zscore = function(x){
  y=(x-mean(x))/sd(x)
  return(y)
}

#prepare log2ed-scaled expdata for plotting sample-feature
scaled_log_fpkm <- function(fpkm,selectgene=NULL) {
	Exp_data <- log2(fpkm + 1)
	if (!is.null(selectgene)) {
       Exp_data <- Exp_data[selectgene, ]
       print("select given genes")
	}

	Exp_data <- apply(Exp_data,1,function(x) zscore(x))
    Exp_data = t(as.matrix(Exp_data))
    return(Exp_data)
}




##the rownames of annotation must match the colnames of fpkm
#distmethod is "1-cor" or "euc"
#plotfile is the prefix 
#1. cal correlation
#cal cor  spearman
cal_cor <- function(fpkm, method="spearman", selectgene=NULL) {
	#normalize data by log2
	Exp_data <- log2(fpkm + 1)
	if (!is.null(selectgene)) {
       Exp_data <- Exp_data[selectgene, ]
       print("select given genes")
	}
	cormat <- round(cor(Exp_data, method=method),4)
	return(cormat)
}

#save mean_cor to file as a summary meta data
#mean_cor <- apply(cormat, 1, mean)
#write.csv(mean_cor, file=mean_corfile, quote=FALSE)


cal_cluster_distance_col <- function(cormat, distmethod="1-cor") {
	##use dist as 1-cormat
	if (distmethod == "1-cor")  {
		hc <- hclust(as.dist(1 - cormat))
		print(plot(hc))
        clustering_distance_cols = as.dist(1 - cormat)
        print("using 1-cor as distance")
    } else if (distmethod == "euc") {
        clustering_distance_cols = dist(t(as.matrix(cormat)), method = "euclidean")
        print("using euclidean as distance")
    } else {
    	print("please give one of the distmethod 1-cor or euc")
    }

    return(clustering_distance_cols)
}

 ##plot sf (sample vs features)
heatmap_sf <- function(fpkm, annotation=NULL,  selectgene=NULL, 
	method="spearman",  distmethod="euc", breaks=NULL, col=NULL, show_rownames=T, 
			show_colnames=F,fontsize=3) {

    Exp_data <- scaled_log_fpkm(fpkm, selectgene)

    cormat <- cal_cor(fpkm, method, selectgene)
    clustering_distance_cols <- cal_cluster_distance_col(cormat, distmethod)
    #clustering_distance_rows <-  clustering_distance_cols
    #pp<- pheatmap(Exp_data,  clustering_distance_cols=clustering_distance_cols, 
     #    clustering_method='ward.D2', cluster_rows=T, cluster_cols=T,
      #   show_rownames=show_rownames, show_colnames=show_colnames,  annotation=annotation, fontsize=fontsize)#,scale='row')

    pp<- pheatmap(Exp_data,  clustering_distance_cols=clustering_distance_cols, clustering_distance_rows=clustering_distance_rows, 
                  clustering_method='ward.D2', cluster_rows=T, cluster_cols=T,
                   fontsize=fontsize)#,scale='row')
    
   return(pp)
}

 ##plot sf no sample cluster
 #samples are sorted by some order (strip_num)
# color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(299)
# breaks = c(seq(-4,4,length=300)) #会把颜色弄得比较均一 给定 -4 到 4 的情况下

# heatmap_sf_noSampleCluster <- function(fpkm, annotation= NA, selectgene=NULL, 
# 	show_rownames=T, show_colnames=F, fontsize=3,	
# 	color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
# 	breaks=NA, border_color = "grey60",
# 	filename = NA, width = NA, height = NA) {

#     Exp_data <- scaled_log_fpkm(fpkm, selectgene)    
        	
#     pp<- pheatmap(Exp_data,  annotation=annotation,
#          clustering_method='ward.D2', cluster_rows=T, cluster_cols=F,
#          show_rownames=show_rownames, show_colnames=show_colnames,  
#          color=color, breaks=breaks, border_color=border_color,
#           fontsize=fontsize,
#          filename = filename, width = width, height = height)#,scale='row')  

#     return(pp)
# }




#plotfile, (sample vs samples)
heatmap_ss <- function(fpkm, annotation, selectgene=NULL, 
	method="spearman",  distmethod="euc", breaks=NULL, col=NULL, show_rownames=T, 
			show_colnames=F,fontsize=3) {
	
    cormat <- cal_cor(fpkm, method, selectgene)
    ###5.distance
	clustering_distance_cols <- cal_cluster_distance_col(cormat, distmethod)

    clustering_distance_rows <-  clustering_distance_cols

    if (!is.null(col) & !is.null(breaks)) {
			p<- pheatmap(cormat_spear, color=col, breaks=breaks, clustering_distance_rows=clustering_distance_rows, 
				clustering_distance_cols=clustering_distance_cols, 
			clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=show_rownames, 
			show_colnames=show_colnames, annotation=annotation, fontsize=fontsize)#,scale='row')
        } else if (!is.null(col) & is.null(breaks)) {
        	p<- pheatmap(cormat_spear, color=col, clustering_distance_rows=clustering_distance_rows, 
				clustering_distance_cols=clustering_distance_cols,  
		    clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=show_rownames, 
		    show_colnames=show_colnames, annotation=annotation, fontsize=fontsize)#,scale='row')
     
        } else if (is.null(col) & !is.null(breaks)) {
            p<- pheatmap(cormat_spear, breaks=breaks, clustering_distance_rows=clustering_distance_rows, 
				clustering_distance_cols=clustering_distance_cols, 
		    clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=show_rownames, 
		    show_colnames=show_colnames, annotation=annotation, fontsize=fontsize)#,scale='row')
     

        } else {
        	p<- pheatmap(cormat_spear, clustering_distance_rows=clustering_distance_rows, 
				clustering_distance_cols=clustering_distance_cols, 
			clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=show_rownames, 
			show_colnames=show_colnames, annotation=annotation, fontsize=fontsize)#,scale='row')
        }
	
		#ggsave(p, filename=plotfile)
        #print(paste("heatmap plots saved in",  plotfile, sep=" "))
        return(p)
}










heatmap_ss_sf <- function(fpkm, annotation, plotfile, selectgene=NULL, 
	method="spearman",  distmethod="euc", breaks=NULL, col=NULL ) {
	#cal cor  spearman
	#normalize data by log2
	Exp_data <- log2(fpkm + 1)
	if (!is.null(selectgene)) {
       Exp_data <- Exp_data[selectgene, ]
       print("select given genes")
	}
	cormat_spear <- round(cor(Exp_data, method=method),4)

    ###5.distance
	##use dist as 1-cormat_spear
	if (distmethod == "1-cor")  {
		hc <- hclust(as.dist(1 - cormat_spear))
		print(plot(hc))
        clustering_distance_rows = as.dist(1 - cormat_spear)
        print("using 1-cor as distance")
    } else if (distmethod == "euc") {
        clustering_distance_rows = dist(t(as.matrix(cormat_spear)), method = "euclidean")
        print("using euclidean as distance")
    } else {
    	print("please give one of the distmethod 1-cor or euc")
    }

    clustering_distance_cols = clustering_distance_rows

    if (!is.null(col) & !is.null(breaks)) {
			p<- pheatmap(cormat_spear, color=col, breaks=breaks, clustering_distance_rows=clustering_distance_rows, 
				clustering_distance_cols=clustering_distance_cols, 
			clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=T, 
			show_colnames=F, annotation=annotation, fontsize=6)#,scale='row')
    } else if (!is.null(col) & is.null(breaks)) {
        	p<- pheatmap(cormat_spear, color=col, clustering_distance_rows=clustering_distance_rows, 
				clustering_distance_cols=clustering_distance_cols,  
		    clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=T, 
		    show_colnames=F, annotation=annotation, fontsize=6)#,scale='row')
     
    } else if (is.null(col) & !is.null(breaks)) {
            p<- pheatmap(cormat_spear, breaks=breaks, clustering_distance_rows=clustering_distance_rows, 
				clustering_distance_cols=clustering_distance_cols, 
		    clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=T, 
		    show_colnames=F, annotation=annotation, fontsize=6)#,scale='row')
     

    } else {
        	p<- pheatmap(cormat_spear, clustering_distance_rows=clustering_distance_rows, 
				clustering_distance_cols=clustering_distance_cols, 
			clustering_method='ward.D2', cluster_rows=T, cluster_cols=T, show_rownames=T, 
			show_colnames=F, annotation=annotation, fontsize=6)#,scale='row')
    }
	    
	plotfile_ss <- paste0(plotfile, ".SS.pdf")
	plotfile_sf <- paste0(plotfile, ".SF.pdf")
	ggsave(p, filename=plotfile_ss)

 ##plot sf

	Exp_data <- apply(Exp_data,1,function(x) zscore(x))
    Exp_data = t(as.matrix(Exp_data))
    pp<- pheatmap(Exp_data,  clustering_distance_cols=clustering_distance_cols, 
         clustering_method='ward.D2', cluster_rows=T, 
         show_rownames=F, cluster_cols=T, annotation=annotation, fontsize=6)#,scale='row')
  
	ggsave(pp, filename=plotfile_sf)
	print(paste("heatmap plots saved in",  plotfile_sf, sep=" "))

	return(list(p,pp))

}
