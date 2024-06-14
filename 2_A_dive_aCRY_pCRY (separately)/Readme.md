RNAseq aCRY pCRY data dive (separate)
================
Kelterborn
2024-04-19

- [Load System](#load-system)
  - [-R_libraries](#-r_libraries)
  - [-R_folders](#-r_folders)
  - [-R functions](#-r-functions)
- [1. Load aCRY & pCRY](#1-load-acry--pcry)
  - [- Tximeta files](#--tximeta-files)
  - [Run Deseq2](#run-deseq2)
- [2. Pre-Analysis](#2-pre-analysis)
- [3. Results](#3-results)
  - [Colors](#colors)
  - [Plot Counts](#plot-counts)
  - [Make results](#make-results)
  - [Volcanos](#volcanos)
- [Export](#export)

BiocManager::install()

# Load System

## -R_libraries

## -R_folders

## -R functions

``` r
getresults_SK <- function(contrast,nres){
r1 <- results(dds, contrast = contrast)
r1$symbol <- mcols(dds)$symbol
assign(paste("res",nres,sep="."),r1)
r1 <- list(r1)
names(r1) <- nres
r1
}
```

``` r
topgenes_f <- function(res,p=0.005,bM=10,l2FC=1){
a <- subset(res, padj < p & baseMean > bM & abs(log2FoldChange) > l2FC)
if(nrow(a)>0) {
a <- a[order(a$baseMean, decreasing = T),]
  a$rank.bm <- seq(1:length(rownames(a)))
a <- a[order(a$padj, decreasing = F),]
  a$rank.padj <- seq(1:length(rownames(a)))
a <- a[order(abs(a$log2FoldChange), decreasing = T),]
  a$rank.l2FC <- seq(1:length(rownames(a)))
a$rank.sum <- a$rank.l2FC+a$rank.bm+a$rank.padj
  a <- a[order(a$rank.sum),]
  }
a
}
```

``` r
# res <- res_ashr_list[[i]][6]

Vulcano_SK <- function(res,
                       n=NULL,
                       list1=NULL,
                       list2=NULL,
                       list1.col=NULL,
                       list2.col=NULL,
                       ol.col=NULL,
                       list1.n=NULL,
                       list2.n=NULL,
                       xlim=NULL,
                       ylim=NULL,
                       
                       ...)
{
# Name
if(is.null(n))  {
  if(is.list(res)){
      n <- names(res)
  } else {
    n <- "results"
  }
  }
  
# unlist results
if(is.list(res)){
  res <- res[[1]]
}
  
# define list 1
if(is.null(list1) ) {
  list1 <- topgenes_f(res=res,p=0.05,bM = 10,l2FC = 1) %>% rownames()
}
# define list 2
if(is.null(list2) ) {
  list2 <- topgenes_f(res=res, p=0.01, l2FC = 2, bM=10) %>% rownames()
}
# list 1 color
if(is.null(list1.col) ) {
  list1.col <- "royalblue1"
}
# list 2 color
if(is.null(list2.col) ) {
  list2.col <- "royalblue4"
}
# list ol color
if(is.null(ol.col) ) {
  ol.col <- "yellow"
} 

# calculate overlap
ol <- calculate.overlap(list(list1,list2))

# name list 1
if(is.null(list1.n) ) {
  list1.n <- "list1"
}
# name list 2
if(is.null(list2.n) ) {
  list2.n <- "list2"
}

# xlim
if(is.null(xlim) ) {
  xlim <- abs(c(min(res$log2FoldChange, na.rm=TRUE),
              max(res$log2FoldChange, na.rm=TRUE))) %>% max()
}

# ylim
if(is.null(ylim) ) {
  ylim <- max(-log10(res$padj), na.rm = TRUE) + 5
}

# colors
pcol <- "black"
lcol <- colours[1+1]

# add symbol
res$symbol <- mcols(dds_list[[i]])$geneSymbol

# define list groups
res$isl1 <- rownames(res) %in% list1
res$isl2 <- rownames(res) %in% list2
res$isol <- rownames(res) %in% ol$a3

# order results (list ol on top)
res <- res[order(abs(res$log2FoldChange), decreasing = T),]
res <- res[order(res$isl1, decreasing = F),]
res <- res[order(res$isl2, decreasing = F),]
res <- res[order(res$isol, decreasing = F),]


# change outlier shape
shape <- 1
if(!is.null(xlim)){
shape <- ifelse(abs(res$log2FoldChange) > xlim, 6,16)
}
if(!is.null(ylim)){
shape[res$padj < 10^-ylim] <- 6
}
# name shape
names(shape)[shape == 6] <- 'outlier'
names(shape)[shape == 16] <- 'in range'
# move outlier to border
res$log2FoldChange[res$log2FoldChange > xlim] <- xlim
res$log2FoldChange[res$log2FoldChange < -xlim] <- -xlim
res$padj[res$padj < 10^-ylim] <- 10^-ylim


# name colors
keyvals <- ifelse(
  res$isol == TRUE, ol.col,
    ifelse(
  res$isl2 == TRUE, list2.col,
    ifelse(
    res$isl1 == TRUE, list1.col,
    'grey') ) )

keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == ol.col] <- "overlap"
  names(keyvals)[keyvals == list1.col] <- list1.n
  names(keyvals)[keyvals == list2.col] <- list2.n
  names(keyvals)[keyvals == 'grey'] <- 'other'

ev_f <- EnhancedVolcano(res,
    x = 'log2FoldChange',
    y = 'padj',
    lab = res$symbol,
    selectLab  = res$symbol[which(names(keyvals) == list2.n | names(keyvals) ==list1.n | names(keyvals) =="overlap")],
    labSize = 3,
    drawConnectors = TRUE,
    boxedLabels = FALSE,
    widthConnectors = 0.5,
    colConnectors = lcol,
    pointSize = 4.0,
    max.overlaps = 100,
    colCustom = keyvals,
    xlim = c(-xlim,xlim),
    ylim = c(0, ylim),
    ylab = "Padj (-log10)",
    title = n,
    subtitle = paste("list1:",length(list1),", list2: ",length(list2)),
    # sub = "SVF",
    # pCutoff = 10^(-60),
    # FCcutoff = 1,
    shapeCustom =shape,
    # pointSize = c(ifelse(rownames(res_WT_D_vs.WT_BL) %in% rownames(top_WT_BL_vs.pcry_BL), 8, 1)),
    legendLabels=c('Not sig.','|L2F| > 2','p-adj < 0.05',
                    'p-adj & L2F'),
    legendPosition = 'bottom',
    col=c('grey', pcol, pcol, pcol),
    legendLabSize = 8,
    legendIconSize = 2.0,
    axisLabSize = 8,
    titleLabSize = 8,
    subtitleLabSize = 8,
    captionLabSize = 8
    )
ev_f

}
```

# 1. Load aCRY & pCRY

## - Tximeta files

## Run Deseq2

![](Readme_files/figure-gfm/dds_acry_pcry-1.png)<!-- -->![](Readme_files/figure-gfm/dds_acry_pcry-2.png)<!-- -->![](Readme_files/figure-gfm/dds_acry_pcry-3.png)<!-- -->![](Readme_files/figure-gfm/dds_acry_pcry-4.png)<!-- -->![](Readme_files/figure-gfm/dds_acry_pcry-5.png)<!-- -->![](Readme_files/figure-gfm/dds_acry_pcry-6.png)<!-- -->![](Readme_files/figure-gfm/dds_acry_pcry-7.png)<!-- -->![](Readme_files/figure-gfm/dds_acry_pcry-8.png)<!-- -->

# 2. Pre-Analysis

### - Data transformations

``` r
# load(file = paste(dir,"rlog_list.RDS",sep="/"))
msdp1 <- lapply(lapply(dds_list,assay,normalized =TRUE),meanSdPlot)
msdp2 <- lapply(lapply(vsd_list,assay,normalized =TRUE),meanSdPlot)
msdp3 <- lapply(lapply(ntd_list,assay,normalized =TRUE),meanSdPlot)
msdp4 <- lapply(lapply(rld_list,assay,normalized =TRUE),meanSdPlot)
```

<img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-1.png" width="50%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-2.png" width="50%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-3.png" width="50%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-4.png" width="50%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-5.png" width="50%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-6.png" width="50%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-7.png" width="50%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-8.png" width="50%" />

### - Check sample distance

<img src="Readme_files/figure-gfm/pre_sample_dist-1.png" width="100%" /><img src="Readme_files/figure-gfm/pre_sample_dist-2.png" width="100%" />

### - Perform principal component analysis

<img src="Readme_files/figure-gfm/pca-1.png" width="80%" /><img src="Readme_files/figure-gfm/pca-2.png" width="80%" />

###### – Advanced PCA

<img src="Readme_files/figure-gfm/pca_advanced-1.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-2.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-3.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-4.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-5.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-6.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-7.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-8.png" width="80%" />

# 3. Results

## Colors

## Plot Counts

![](Readme_files/figure-gfm/plot_counts-1.png)<!-- -->![](Readme_files/figure-gfm/plot_counts-2.png)<!-- -->![](Readme_files/figure-gfm/plot_counts-3.png)<!-- -->![](Readme_files/figure-gfm/plot_counts-4.png)<!-- -->![](Readme_files/figure-gfm/plot_counts-5.png)<!-- -->![](Readme_files/figure-gfm/plot_counts-6.png)<!-- -->

## Make results

## Volcanos

### Volcanos 2.0

#### top genes

``` r
i <- 1
ev1_d <- Vulcano_SK(res=res_ashr_list[[i]][5],
                    xlim=5,ylim=10,
                    list2.col = "black",list1.col = "grey35")
ev1_b <- Vulcano_SK(res=res_ashr_list[[i]][6],
                    xlim=5,ylim=10,
                    list2.col = "royalblue4",list1.col = "royalblue1")
ev1_r <- Vulcano_SK(res=res_ashr_list[[i]][7],
                    xlim=5,ylim=10,
                    list2.col = "firebrick4",list1.col = "firebrick1")

i <- 2
ev2_d <- Vulcano_SK(res=res_ashr_list[[i]][5],
                    xlim=5,ylim=50,
                    list2.col = "black",list1.col = "grey35")
ev2_b <- Vulcano_SK(res=res_ashr_list[[i]][6],
                    xlim=5,ylim=50,
                    list2.col = "royalblue4",list1.col = "royalblue1")
ev2_r <- Vulcano_SK(res=res_ashr_list[[i]][7],
                    xlim=5,ylim=50,
                    list2.col = "firebrick4",list1.col = "firebrick1")

ev1_d+ev1_b+ev1_r  +
  plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'bottom', axis.title=element_text(size=8))
```

![](Readme_files/figure-gfm/top_volcano2-1.png)<!-- -->

``` r
ev2_d+ev2_b+ev2_r  +
  plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'bottom', axis.title=element_text(size=8))
```

![](Readme_files/figure-gfm/top_volcano2-2.png)<!-- -->

``` r
# ga <- ggarrange(evb,
#                 evd+rremove("ylab"),
#                 evr+rremove("ylab"),
#                 ncol = 3, nrow = 2,common.legend = TRUE, legend = "bottom",
#                 widths = c(1.05,1,1))
# 
# 
# ggsave(ga, file="figures/Vulcano_all.pdf",
#        width = 12,
#        height = 8)
```

#### interaction in D-D, Bl-Bl,R-R

``` r
i <- 1
ev1_d <- Vulcano_SK(res=res_ashr_list[[i]][5],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=10)
ev1_b <- Vulcano_SK(res=res_ashr_list[[i]][6],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=10)
ev1_r <- Vulcano_SK(res=res_ashr_list[[i]][7],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=10)

i <- 2
ev2_d <- Vulcano_SK(res=res_ashr_list[[i]][5],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=50)
ev2_b <- Vulcano_SK(res=res_ashr_list[[i]][6],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=50)
ev2_r <- Vulcano_SK(res=res_ashr_list[[i]][7],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=50)

ev1_d+ev1_b+ev1_r  +
  plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'bottom', axis.title=element_text(size=8))
```

![](Readme_files/figure-gfm/vsvs_volcano2-1.png)<!-- -->

``` r
ev2_d+ev2_b+ev2_r  +
  plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'bottom', axis.title=element_text(size=8))
```

![](Readme_files/figure-gfm/vsvs_volcano2-2.png)<!-- -->

#### interaction in BL-D, R-D

``` r
deg_list[[1]][8:9] %>% names()
```

    ## [1] "deg_acry_BLvD.vs.WT_BLvD" "deg_acry_RvD.vs.WT_RvD"

``` r
i <- 1
ev1.1 <- Vulcano_SK(res=res_ashr_list[[i]][1],
                    list1=top_list[[i]][[8]], list2=top_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5, ylim=50)
ev1.2 <- Vulcano_SK(res=res_ashr_list[[i]][2],
                    list1=top_list[[i]][[8]], list2=top_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=50)
ev1.3 <- Vulcano_SK(res=res_ashr_list[[i]][3],
                    list1=top_list[[i]][[8]], list2=top_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=50)
ev1.4 <- Vulcano_SK(res=res_ashr_list[[i]][4],
                    list1=top_list[[i]][[8]], list2=top_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=50)

i <- 2
ev2.1 <- Vulcano_SK(res=res_ashr_list[[i]][1],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=100)
ev2.2 <- Vulcano_SK(res=res_ashr_list[[i]][2],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=100)
ev2.3 <- Vulcano_SK(res=res_ashr_list[[i]][3],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=100)
ev2.4 <- Vulcano_SK(res=res_ashr_list[[i]][4],
                    list1=deg_list[[i]][[8]], list2=deg_list[[i]][[9]],
                    list1.n="blue", list2.n="red",
                    list1.col = "royalblue1",list2.col = "firebrick1",ol.col="orchid",
                    xlim=5,ylim=100)

ev1.1 + ev1.2 + ev1.3 + ev1.4 +
  plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'bottom', axis.title=element_text(size=8))
```

![](Readme_files/figure-gfm/vsvs2_volcano2-1.png)<!-- -->

``` r
ev2.1 + ev2.2 + ev2.3 + ev2.4 +
  plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'bottom', axis.title=element_text(size=8))
```

![](Readme_files/figure-gfm/vsvs2_volcano2-2.png)<!-- -->

``` r
ev1.1 + ev1.2 + ev2.1 + ev2.2 +
  plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'bottom', axis.title=element_text(size=8))
```

![](Readme_files/figure-gfm/vsvs2_volcano2-3.png)<!-- -->

### dark

``` r
# dark
res <- res_ashr_list[[i]][[5]]
n <- 'WT vs. pCRY in dark'
nf <- 'figures/Volcano_dark.png'
l <- nrow(top_list[[i]][['pCRY_D.vs.WT_D']])
pcol <- "black"
lcol <- colours[1+1]


rownames(res) <- mcols(dds_list[[i]])$id.symbol
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$log2FoldChange),]

shape <- ifelse(abs(res$log2FoldChange) > 5, 6,
                ifelse(res$padj < 1e-100,6,16))
summary(is.na(shape))
# shape[is.na(shape)] <- 2
names(shape)[shape == 6] <- 'outliner'
names(shape)[shape == 16] <- 'in range'

res$log2FoldChange[res$log2FoldChange > 5] <- 5
res$log2FoldChange[res$log2FoldChange < -5] <- -5
res$padj[res$padj < 1e-100] <- 1e-100
summary(res$padj < 1e-100)

keyvals <- ifelse(
    res$ishs == TRUE, pcol,
      'grey')

keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == pcol] <- 'highly significant'
  names(keyvals)[keyvals == 'grey'] <- 'other'
  
evd <- EnhancedVolcano(res,
    x = 'log2FoldChange',
    y = 'padj',
    lab = rownames(res),
    selectLab  = rownames(res)[which(names(keyvals) == 'highly significant')],
    # selectLab = anno_tf$geneSymbol,
    # boxedLabels = TRUE,
    labSize = 3,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.5,
    colConnectors = lcol,
    pointSize = 4.0,
    max.overlaps = 10,
    colCustom = keyvals,
    xlim = c(-5.1, 5.1),
    ylim = c(0, 101),
    ylab = "Padj (-log10)",
    title = n,
    subtitle = paste("DE genes:",l),
    # sub = "SVF",
    pCutoff = 10^(-60),
    FCcutoff = 3,
   shapeCustom =shape,
    # pointSize = c(ifelse(rownames(res_WT_D_vs.WT_BL) %in% rownames(top_WT_BL_vs.pcry_BL), 8, 1)),
    legendLabels=c('Not sig.','|L2F| > 2.5','p-adj < 0.05',
                   'p-adj & L2F'),
    legendPosition = 'bottom',
    col=c('grey', pcol, pcol, pcol)
    )
evd
```

### blue

``` r
# blue
res <- res2_shrink
n <- 'WT vs. pCRY in blue'
nf <- 'figures/Volcano_blue.png'
l <- nrow(top.res2)
pcol <- colours[3]
lcol <- colours[3+1]


rownames(res) <- mcols(dds)$id.symbol
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$log2FoldChange),]

shape <- ifelse(abs(res$log2FoldChange) > 5, 17,
                ifelse(res$padj < 1e-100,17,16))
summary(is.na(shape))
# shape[is.na(shape)] <- 2
names(shape)[shape == 6] <- 'outliner'
names(shape)[shape == 16] <- 'in range'

res$log2FoldChange[res$log2FoldChange > 5] <- 5
res$log2FoldChange[res$log2FoldChange < -5] <- -5
res$padj[res$padj < 1e-100] <- 1e-100
summary(res$padj < 1e-100)

keyvals <- ifelse(
    res$ishs == TRUE, pcol,
      'grey')

keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == pcol] <- 'highly significant'
  names(keyvals)[keyvals == 'grey'] <- 'other'
  
evb <- EnhancedVolcano(res,
    x = 'log2FoldChange',
    y = 'padj',
    lab = rownames(res),
    selectLab  = rownames(res)[which(names(keyvals) == 'highly significant')],
    # selectLab = anno_tf$geneSymbol,
    # boxedLabels = TRUE,
    labSize = 3,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.5,
    colConnectors = lcol,
    pointSize = 4.0,
    max.overlaps = 10,
    colCustom = keyvals,
    xlim = c(-5.1, 5.1),
    ylim = c(0, 101),
    ylab = "Padj (-log10)",
    title = n,
    subtitle = paste("DE genes:",l),
    # sub = "SVF",
    pCutoff = 10^(-50),
    FCcutoff = 2.5,
   shapeCustom =shape,
    # pointSize = c(ifelse(rownames(res_WT_D_vs.WT_BL) %in% rownames(top_WT_BL_vs.pcry_BL), 8, 1)),
    legendLabels=c('Not sig.','|L2F| > 2.5','p-adj < 0.05',
                   'p-adj & L2F'),
    legendPosition = 'bottom',
    col=c('grey', pcol, pcol, pcol)
    )
evb
```

### red

``` r
# red
res <- res3_shrink
n <- 'WT vs. pCRY in red'
nf <- 'figures/Volcano_r.png'
l <- nrow(top.res3)
pcol <- colours[5]
lcol <- colours[5+1]

rownames(res) <- mcols(dds)$id.symbol
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$log2FoldChange),]

shape <- ifelse(abs(res$log2FoldChange) > 5, 6,
                ifelse(res$padj < 1e-100,6,16))
summary(is.na(shape))
# shape[is.na(shape)] <- 2
names(shape)[shape == 6] <- 'outliner'
names(shape)[shape == 16] <- 'in range'

res$log2FoldChange[res$log2FoldChange > 5] <- 5
res$log2FoldChange[res$log2FoldChange < -5] <- -5
res$padj[res$padj < 1e-100] <- 1e-100
summary(res$padj < 1e-100)

keyvals <- ifelse(
    res$ishs == TRUE, pcol,
      'grey')

keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == pcol] <- 'highly significant'
  names(keyvals)[keyvals == 'grey'] <- 'other'
  
evr <- EnhancedVolcano(res,
    x = 'log2FoldChange',
    y = 'padj',
    lab = rownames(res),
    selectLab  = rownames(res)[which(names(keyvals) == 'highly significant')],
    # selectLab = anno_tf$geneSymbol,
    # boxedLabels = TRUE,
    labSize = 3,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.5,
    colConnectors = lcol,
    pointSize = 4.0,
    max.overlaps = 10,
    colCustom = keyvals,
    xlim = c(-5.1, 5.1),
    ylim = c(0, 101),
    ylab = "Padj (-log10)",
    title = n,
    subtitle = paste("DE genes:",l),
    # sub = "SVF",
    pCutoff = 10^(-50),
    FCcutoff = 2.5,
   shapeCustom =shape,
    # pointSize = c(ifelse(rownames(res_WT_D_vs.WT_BL) %in% rownames(top_WT_BL_vs.pcry_BL), 8, 1)),
    legendLabels=c('Not sig.','|L2F| > 2.5','p-adj < 0.05',
                   'p-adj & L2F'),
    legendPosition = 'bottom',
    col=c('grey', pcol, pcol, pcol)
    )
evr

#ggsave(file="figures/Vulcano_rtest.pdf",
#       width = 12,
#       height = 12)
```

### all

``` r
ga <- ggarrange(evb,
                evd+rremove("ylab"),
                evr+rremove("ylab"),
                ncol = 3, nrow = 1,common.legend = TRUE, legend = "bottom",
                widths = c(1.05,1,1))


ggsave(ga, file="figures/Vulcano_all.pdf",
       width = 12,
       height = 8)

# + ggplot2::coord_cartesian(xlim=c(-10, 10))
# ggplot2::coord_cartesian(xlim=c(-5, 5)




# with ggplot
ggplot(as.data.frame(resLFC_WT_BL_vs.pcry_BL),aes(log2FoldChange,-log(padj))) +
  geom_point() +
  geom_point(data = as.data.frame(resLFC_WT_BL_vs.pcry_BL[rownames(top_WT_BL_vs.pcry_BL),]), colour = "blue") + 
  geom_point(data = as.data.frame(resLFC_WT_BL_vs.pcry_BL[rocs,]), colour = "red") +
  geom_text_repel(data = as.data.frame(resLFC_WT_BL_vs.pcry_BL[rocs,]), aes(label=anno[rocs,"id.symbol"]),colour = "red",hjust=-0.5, vjust=-1) +
  xlim(-5,5)
getwd()
ggsave("graphs/Vulcano ROCs-2.pdf",
       width = 10,
       height = 6)

## enhanced vulcano
library(airway)
  library(magrittr)

 library('DESeq2')

  dds <- DESeqDataSet(airway, design = ~ cell + dex)
  dds <- DESeq(dds, betaPrior=FALSE)
  res <- results(dds,
    contrast = c('dex','trt','untrt'))
  res <- lfcShrink(dds,
    contrast = c('dex','trt','untrt'), res=res, type = 'normal')
 EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
```

#### —————————

#### COPY & PASTE from pCRY

#### Make results

### combine L2FC with padj

``` r
res3_shrink
subset(mcols(dds), id.symbol =="SELU1") %>% rownames()
res3_shrink["Cre06.g263550",]

# DataFrame with 1 row and 5 columns
#                baseMean log2FoldChange     lfcSE      pvalue        padj
#               <numeric>      <numeric> <numeric>   <numeric>   <numeric>
# Cre06.g263550   3566.14        2.90397  0.167321 1.72033e-71 1.97322e-68

plot(res3_shrink$log2FoldChange,-log(res3_shrink$padj))
res3_shrink$logpadj <- -log(res3_shrink$padj, base=10)

-log(1.38379e-06,base = 10)

x=-2.5
60/(3^2)

y=-8*x^2 + 50
y

res3_shrink$ishs <- ifelse((-20/3*(res3_shrink$log2FoldChange^2) + 60) < res3_shrink$logpadj,TRUE,FALSE)
summary(res3_shrink$ishs)
subset(res3_shrink,ishs==TRUE & logpadj< 50 & abs(log2FoldChange)<2.5)
subset(res3_shrink,ishs==TRUE & logpadj< 50)

# for log2 -/+ 5 and padj 100 
res1_shrink$ishs <- ifelse((-20/3*(res1_shrink$log2FoldChange^2) + 60) < -log(res1_shrink$padj, base=10),TRUE,FALSE)
res2_shrink$ishs <- ifelse((-20/3*(res2_shrink$log2FoldChange^2) + 60) < -log(res2_shrink$padj, base=10),TRUE,FALSE)
res3_shrink$ishs <- ifelse((-20/3*(res3_shrink$log2FoldChange^2) + 60) < -log(res3_shrink$padj, base=10),TRUE,FALSE)
res4_shrink$ishs <- ifelse((-20/3*(res4_shrink$log2FoldChange^2) + 60) < -log(res4_shrink$padj, base=10),TRUE,FALSE)
res5_shrink$ishs <- ifelse((-20/3*(res5_shrink$log2FoldChange^2) + 60) < -log(res5_shrink$padj, base=10),TRUE,FALSE)
res6_shrink$ishs <- ifelse((-20/3*(res6_shrink$log2FoldChange^2) + 60) < -log(res6_shrink$padj, base=10),TRUE,FALSE)
res7_shrink$ishs <- ifelse((-20/3*(res7_shrink$log2FoldChange^2) + 60) < -log(res7_shrink$padj, base=10),TRUE,FALSE)
```

``` r
outdir
figures <- paste(outdir,"figures",sep="/")
dir.create(figures)
```

### Colours

``` r
Greys <- brewer.pal(9, name="Greys")[c(7,5)]
Paired <- brewer.pal(12, name="Paired")[c(2,1,6,5)]
colours <- c(Greys,Paired)
```

### PCA (Panel B)

``` r
library("PCAtools")
library(ggalt)

vsd <- vst(dds, blind=FALSE)
length(vsd)

dubs <- duplicated(mcols(vsd)$id.symbol)
length(dubs)

unique <- unique(mcols(vsd)$id.symbol)
length(unique)

vsd2 <- vsd[!duplicated(mcols(vsd)$id.symbol)]
rownames(vsd2) <- mcols(vsd2)$id.symbol

# colours
Greys <- brewer.pal(9, name="Greys")[c(7,5)]
Paired <- brewer.pal(12, name="Paired")[c(2,1,6,5)]
colours <- c(Greys,Paired)
group.colors <- colours
group.colors <- c("#525252", "#969696", "#1F78B4", "#A6CEE3", "#E31A1C", "#FB9A99")
names(group.colors) <- levels(colData(dds)$condition)
group.colors      

plotPCA(vsd2, intgroup="condition")
vst <- assay(vsd2)
p <- pca(vst, metadata = colData(dds), removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, colby = "condition")
biplot(p, showLoadings = TRUE,labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

biplot(p,x= "PC1", y="PC3",
       lab = p$metadata$sample,
    colby = 'condition', colkey = group.colors,
    encircle = TRUE,
    encircleFill = TRUE,
    legendPosition = 'right',
    labSize = 4)
ggsave(paste(figures,"2023_12_PCA1-3.png",sep="/"),
  width = 6,
  height = 6)

biplot(p,x= "PC3", y="PC4",
       lab = p$metadata$sample,
    colby = 'condition', colkey = group.colors,
    encircle = TRUE,
    encircleFill = TRUE,
    legendPosition = 'right')
ggsave(paste(figures,"2023_12_PCA3-4.pdf",sep="/"),
  width = 12,
  height = 12)

pairsplot(p,colby = 'condition', colkey = group.colors)
ggsave(paste(figures,"2023_12_pairsplot.pdf",sep="/"),
  width = 12,
  height = 12)


plotloadings(p,
    rangeRetain = 0.01,
    labSize = 2.0,
    title = 'Loadings plot',
    subtitle = 'PC1, PC2, PC3, PC4, PC5',
    caption = 'Top 1% variables',
    drawConnectors = TRUE)
ggsave(paste(figures,"2023_12_loadings.pdf",sep="/"),
  width = 12,
  height = 12)

# Determine optimum number of PCs to retain
elbow <- findElbowPoint(p$variance)
elbow

anno["Cre01.g016600",]

 biplot(p,
    lab = p$metadata$condition,
    colby = 'condition', colkey = group.colors,
    hline = 0, vline = 0,
    encircle = TRUE, encircleFill = TRUE,
    legendLabSize = 16, legendIconSize = 8.0,
    legendPosition = 'right')

pairsplot(p,
    components = getComponents(p, c(1:6)),
    triangle = TRUE, trianglelabSize = 12,
    hline = 0, vline = 0,
    pointSize = 2,
    gridlines.major = FALSE, gridlines.minor = FALSE,
    colby = 'condition', colkey = group.colors,
    title = 'Pairs plot', plotaxes = FALSE,
    margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))
```

### Counts

``` r
# TF genes

goi_tf <- c("CCM1","LCR1", "HY5", "QER7", "QER4", "QER6", "ROC15", "ROC40", "ROC66", "ROC75", "ROC114", "ROC55", "CON1", "CRB1")
anno[c("Cre17.g745697","Cre13.g567250","Cre01.g043550","Cre02.g079550"),"geneSymbol"] <- c("QER4","QER6","QER7","ROC110")
anno_tf <- subset(anno,geneSymbol %in% goi_tf | gene_id %in% c("Cre17.g745697","Cre13.g567250","Cre01.g043550","Cre02.g079550"))
anno_tf[,c("gene_id", "geneSymbol","id.symbol")]
summary(anno_tf$gene_id %in% rownames(dds))
anno_tf[!anno_tf$gene_id %in% rownames(dds),]
anno_tf <- anno_tf[anno_tf$gene_id %in% rownames(dds),]

# Combined counts table TFs

goi <- anno_tf
l <- nrow(goi)
all_counts <- {}
for (i in 1:l){
  d <-  plotCounts(dds, gene=goi[i,"gene_id"], intgroup="condition", col=col,main=res$SYMBOL[i],returnData=TRUE)
  d$Gene <- rep(goi[i,"geneSymbol"],length(rownames(d)))
  d$sample <- rownames(d)
  rownames(d) <- {}
  all_counts <- bind_rows(all_counts,d)
  }

all_counts$Gene
levels(all_counts$condition)
all_counts$Gene <- factor(all_counts$Gene)
levels(all_counts$Gene)

all_counts$Gene <- factor(all_counts$Gene, levels = c("CON1", "ROC40", "CRB1" ,  "QER4" , "ROC110",  "QER7" ,"ROC114" , "ROC75", "CCM1" , "ROC55" , "LCR1",  "QER6" , "ROC66",    "ROC15")  )
levels(all_counts$Gene)


# colours
Greys <- brewer.pal(9, name="Greys")[c(7,5)]
Paired <- brewer.pal(12, name="Paired")[c(2,1,6,5)]
colours <- c(Greys,Paired)
group.colors <- colours
group.colors <- c("#525252", "#969696", "#1F78B4", "#A6CEE3", "#E31A1C", "#FB9A99")
max_val <- 1.0*max(all_counts$count)

# Plot
gcounts <- ggplot(all_counts, aes(x = Gene, y = count, fill=condition)) +
  geom_boxplot(fatten = 1) +
  scale_fill_manual(values = group.colors) +
  scale_y_continuous(trans = "log2")
gcounts
ggsave(paste(figures,"2023_12_counts_tfs_aio_log2.pdf",sep="/"), plot = gcounts,
  width = 12,
  height = 8)

## Phot genes all-in-one
goi_phot <- c("CHR1", "CHR2", "PCRY1", "ACRY1", "DCRY1", "PHOT1", "UVR8", "HKR1") #"HKR2" "DCRY2", 

anno_phot <- subset(anno,geneSymbol %in% goi_phot, drop = FALSE)
rownames(anno_phot) <- anno_phot$geneSymbol
anno_phot <- anno_phot[goi_phot,]

# Combined counts table TFs

goi <- anno_phot
l <- nrow(goi)
all_counts <- {}
for (i in 1:l){
  d <-  plotCounts(dds, gene=goi[i,"gene_id"], intgroup="condition", col=col,main=res$SYMBOL[i],returnData=TRUE)
  d$Gene <- rep(goi[i,"geneSymbol"],length(rownames(d)))
  d$sample <- rownames(d)
  rownames(d) <- {}
  all_counts <- bind_rows(all_counts,d)
  }

all_counts$Gene
levels(all_counts$condition)
levels(all_counts$Gene)

all_counts$Gene <- factor(all_counts$Gene)
levels(all_counts$Gene)

all_counts$Gene <- factor(all_counts$Gene, levels = c("PHOT1","CHR1","CHR2","HKR1","UVR8","DCRY1","PCRY1","ACRY1")) # "DCRY2", 
levels(all_counts$Gene)

# colours
Greys <- brewer.pal(9, name="Greys")[c(7,5)]
Paired <- brewer.pal(12, name="Paired")[c(2,1,6,5)]
colours <- c(Greys,Paired)
group.colors <- colours
group.colors <- c("#525252", "#969696", "#1F78B4", "#A6CEE3", "#E31A1C", "#FB9A99")
max_val <- 1.0*max(all_counts$count)

# Plot
gcounts <- ggplot(all_counts, aes(x = Gene, y = count, fill=condition)) +
  geom_boxplot(fatten = 1) +
  scale_fill_manual(values = group.colors) +
  scale_y_continuous(trans = "log2")
gcounts
ggsave(paste(figures,"2023_12_counts_phots_aio_log2.pdf",sep="/"), plot = gcounts,
  width = 12,
  height = 8)
```

### Volcanos

### dark

``` r
# dark
res <- res1_shrink
n <- 'WT vs. pCRY in dark'
nf <- 'figures/Volcano_dark.png'
l <- nrow(top.res1)
pcol <- "black"
lcol <- colours[1+1]


rownames(res) <- mcols(dds)$id.symbol
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$log2FoldChange),]

shape <- ifelse(abs(res$log2FoldChange) > 5, 6,
                ifelse(res$padj < 1e-100,6,16))
summary(is.na(shape))
# shape[is.na(shape)] <- 2
names(shape)[shape == 6] <- 'outliner'
names(shape)[shape == 16] <- 'in range'

res$log2FoldChange[res$log2FoldChange > 5] <- 5
res$log2FoldChange[res$log2FoldChange < -5] <- -5
res$padj[res$padj < 1e-100] <- 1e-100
summary(res$padj < 1e-100)

keyvals <- ifelse(
    res$ishs == TRUE, pcol,
      'grey')

keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == pcol] <- 'highly significant'
  names(keyvals)[keyvals == 'grey'] <- 'other'
  
evd <- EnhancedVolcano(res,
    x = 'log2FoldChange',
    y = 'padj',
    lab = rownames(res),
    selectLab  = rownames(res)[which(names(keyvals) == 'highly significant')],
    # selectLab = anno_tf$geneSymbol,
    # boxedLabels = TRUE,
    labSize = 3,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.5,
    colConnectors = lcol,
    pointSize = 4.0,
    max.overlaps = 10,
    colCustom = keyvals,
    xlim = c(-5.1, 5.1),
    ylim = c(0, 101),
    ylab = "Padj (-log10)",
    title = n,
    subtitle = paste("DE genes:",l),
    # sub = "SVF",
    pCutoff = 10^(-60),
    FCcutoff = 3,
   shapeCustom =shape,
    # pointSize = c(ifelse(rownames(res_WT_D_vs.WT_BL) %in% rownames(top_WT_BL_vs.pcry_BL), 8, 1)),
    legendLabels=c('Not sig.','|L2F| > 2.5','p-adj < 0.05',
                   'p-adj & L2F'),
    legendPosition = 'bottom',
    col=c('grey', pcol, pcol, pcol)
    )
evd
```

### blue

``` r
# blue
res <- res2_shrink
n <- 'WT vs. pCRY in blue'
nf <- 'figures/Volcano_blue.png'
l <- nrow(top.res2)
pcol <- colours[3]
lcol <- colours[3+1]


rownames(res) <- mcols(dds)$id.symbol
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$log2FoldChange),]

shape <- ifelse(abs(res$log2FoldChange) > 5, 17,
                ifelse(res$padj < 1e-100,17,16))
summary(is.na(shape))
# shape[is.na(shape)] <- 2
names(shape)[shape == 6] <- 'outliner'
names(shape)[shape == 16] <- 'in range'

res$log2FoldChange[res$log2FoldChange > 5] <- 5
res$log2FoldChange[res$log2FoldChange < -5] <- -5
res$padj[res$padj < 1e-100] <- 1e-100
summary(res$padj < 1e-100)

keyvals <- ifelse(
    res$ishs == TRUE, pcol,
      'grey')

keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == pcol] <- 'highly significant'
  names(keyvals)[keyvals == 'grey'] <- 'other'
  
evb <- EnhancedVolcano(res,
    x = 'log2FoldChange',
    y = 'padj',
    lab = rownames(res),
    selectLab  = rownames(res)[which(names(keyvals) == 'highly significant')],
    # selectLab = anno_tf$geneSymbol,
    # boxedLabels = TRUE,
    labSize = 3,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.5,
    colConnectors = lcol,
    pointSize = 4.0,
    max.overlaps = 10,
    colCustom = keyvals,
    xlim = c(-5.1, 5.1),
    ylim = c(0, 101),
    ylab = "Padj (-log10)",
    title = n,
    subtitle = paste("DE genes:",l),
    # sub = "SVF",
    pCutoff = 10^(-50),
    FCcutoff = 2.5,
   shapeCustom =shape,
    # pointSize = c(ifelse(rownames(res_WT_D_vs.WT_BL) %in% rownames(top_WT_BL_vs.pcry_BL), 8, 1)),
    legendLabels=c('Not sig.','|L2F| > 2.5','p-adj < 0.05',
                   'p-adj & L2F'),
    legendPosition = 'bottom',
    col=c('grey', pcol, pcol, pcol)
    )
evb
```

### red

``` r
# red
res <- res3_shrink
n <- 'WT vs. pCRY in red'
nf <- 'figures/Volcano_r.png'
l <- nrow(top.res3)
pcol <- colours[5]
lcol <- colours[5+1]

rownames(res) <- mcols(dds)$id.symbol
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$log2FoldChange),]

shape <- ifelse(abs(res$log2FoldChange) > 5, 6,
                ifelse(res$padj < 1e-100,6,16))
summary(is.na(shape))
# shape[is.na(shape)] <- 2
names(shape)[shape == 6] <- 'outliner'
names(shape)[shape == 16] <- 'in range'

res$log2FoldChange[res$log2FoldChange > 5] <- 5
res$log2FoldChange[res$log2FoldChange < -5] <- -5
res$padj[res$padj < 1e-100] <- 1e-100
summary(res$padj < 1e-100)

keyvals <- ifelse(
    res$ishs == TRUE, pcol,
      'grey')

keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == pcol] <- 'highly significant'
  names(keyvals)[keyvals == 'grey'] <- 'other'
  
evr <- EnhancedVolcano(res,
    x = 'log2FoldChange',
    y = 'padj',
    lab = rownames(res),
    selectLab  = rownames(res)[which(names(keyvals) == 'highly significant')],
    # selectLab = anno_tf$geneSymbol,
    # boxedLabels = TRUE,
    labSize = 3,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.5,
    colConnectors = lcol,
    pointSize = 4.0,
    max.overlaps = 10,
    colCustom = keyvals,
    xlim = c(-5.1, 5.1),
    ylim = c(0, 101),
    ylab = "Padj (-log10)",
    title = n,
    subtitle = paste("DE genes:",l),
    # sub = "SVF",
    pCutoff = 10^(-50),
    FCcutoff = 2.5,
   shapeCustom =shape,
    # pointSize = c(ifelse(rownames(res_WT_D_vs.WT_BL) %in% rownames(top_WT_BL_vs.pcry_BL), 8, 1)),
    legendLabels=c('Not sig.','|L2F| > 2.5','p-adj < 0.05',
                   'p-adj & L2F'),
    legendPosition = 'bottom',
    col=c('grey', pcol, pcol, pcol)
    )
evr

#ggsave(file="figures/Vulcano_rtest.pdf",
#       width = 12,
#       height = 12)
```

### all

``` r
ga <- ggarrange(evb,
                evd+rremove("ylab"),
                evr+rremove("ylab"),
                ncol = 3, nrow = 1,common.legend = TRUE, legend = "bottom",
                widths = c(1.05,1,1))


ggsave(ga, file="figures/Vulcano_all.pdf",
       width = 12,
       height = 8)

# + ggplot2::coord_cartesian(xlim=c(-10, 10))
# ggplot2::coord_cartesian(xlim=c(-5, 5)




# with ggplot
ggplot(as.data.frame(resLFC_WT_BL_vs.pcry_BL),aes(log2FoldChange,-log(padj))) +
  geom_point() +
  geom_point(data = as.data.frame(resLFC_WT_BL_vs.pcry_BL[rownames(top_WT_BL_vs.pcry_BL),]), colour = "blue") + 
  geom_point(data = as.data.frame(resLFC_WT_BL_vs.pcry_BL[rocs,]), colour = "red") +
  geom_text_repel(data = as.data.frame(resLFC_WT_BL_vs.pcry_BL[rocs,]), aes(label=anno[rocs,"id.symbol"]),colour = "red",hjust=-0.5, vjust=-1) +
  xlim(-5,5)
getwd()
ggsave("graphs/Vulcano ROCs-2.pdf",
       width = 10,
       height = 6)

## enhanced vulcano
library(airway)
  library(magrittr)

 library('DESeq2')

  dds <- DESeqDataSet(airway, design = ~ cell + dex)
  dds <- DESeq(dds, betaPrior=FALSE)
  res <- results(dds,
    contrast = c('dex','trt','untrt'))
  res <- lfcShrink(dds,
    contrast = c('dex','trt','untrt'), res=res, type = 'normal')
 EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
```

### Heatmap

### plot heatmap

``` r
library("pheatmap")
ntd <- normTransform(dds)
top.all.hm <- head(top.all, n=100)
df <- assay(ntd)[rownames(top.all.hm),]
rownames(df) <- mcols(dds)[rownames(df),"id.symbol"]

# anno_col <- as.data.frame(colData(dds)[,c("condition","treatment","genotype")])
anno_col <- data.frame(condition = dds$condition, treatment = dds$treatment)
rownames(anno_col) <- rownames(colData(dds))

# column order
anno_col <- arrange(anno_col, match(anno_col$condition, levels(anno_col$condition)))
levels(anno_col$condition)
df <- df[,rownames(anno_col)]

annot_colors=list(condition=c(Cancer="#F0978D",Healthy="#63D7DE"))
Greys <- brewer.pal(9, name="Greys")[c(7,5)]
Paired <- brewer.pal(12, name="Paired")[c(2,1,6,5)]
colours <- c(Greys,Paired)

mycolors2 <- list(genotype = c("#525252","#969696"),
     treatment = c("#525252","#1F78B4","#E31A1C"),
     condition = colours)

names(mycolors2$genotype) <- levels(anno_col$genotype)
names(mycolors2$treatment) <- levels(anno_col$treatment)
names(mycolors2$condition) <- levels(anno_col$condition)
mycolors2

xx <- pheatmap(df, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE,
         annotation_col=anno_col,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         annotation_colors=mycolors2
         )

ggsave("figures/2023_12_Heatmap_top.png",plot=xx,
       width = 5,
       height = 15)
```

### Complex Heatmap

``` r
df
genes <- anno_tf

# collapse replicates
colData(dds)
dds_heat <- collapseReplicates(dds, dds$condition,dds$sample)
dds_heat$sample <- factor(dds_heat$sample, levels=c("WT-D1","pcry-D1","WT-BL1","pcry-BL1","WT-R1","pcry-R1"))
colData(dds_heat)
dds_heat$runsCollapsed
colnames(dds_heat)
matchFirstLevel <- dds_heat$sample == levels(dds_heat$sample)[1]
stopifnot(all(rowSums(counts(dds_heat[,matchFirstLevel])) == counts(dds_heat[,1])))

mat <- counts(dds_heat)
length(rownames(counts(dds_heat)))
length(rownames(df))
rownames(mat) <- mcols(dds_heat)$id.symbol
length(rownames(mat))

mat <- mat[df,]

position <- mat
l <- length(rownames(mat))
pos <- data.frame(genes = rownames(mat),pos = 1:l)
rownames(pos) <- rownames(mat)
pos <- pos[genes,2]

# mat <- cbind(res,counts(dds_heat))
# mat <- subset(mat,!SYMBOL2=="" & SYMBOL2 %in% df$SYMBOL2)
# rownames(mat) <- mat$SYMBOL2
# mat <- mat[,c(12:15)]
heat <- t(scale(t(mat)))
Heatmap(heat)

summary(heat)
myCol <- colorRampPalette(c('brown3', 'yellow2', 'darkolivegreen3'))(100)
myBreaks <- seq(-1.5, 1.5, length.out = 100)

ann <- data.frame(condition = colData(dds_heat)[,"condition"],
  stringsAsFactors = FALSE)

colAnn <- HeatmapAnnotation(
  df = ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  # col = colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'))

boxplotCol <- HeatmapAnnotation(
  boxplot = anno_boxplot(
    heat,
    border = FALSE,
    gp = gpar(fill = '#CCCCCC'),
    pch = '.',
    size = unit(2, 'mm'),
    axis = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 12),
      side = 'left')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'col')

boxplotRow <- HeatmapAnnotation(
  boxplot = row_anno_boxplot(
    heat,
    border = FALSE,
    gp = gpar(fill = '#CCCCCC'),
    pch = '.',
    size = unit(2, 'mm'),
    axis = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 12),
      side = 'top')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'row')

genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = pos,
    labels = genes,
    labels_gp = gpar(fontsize = 15, fontface = 'bold'),
    padding = 0.75),
  width = unit(2.0, 'cm') +
    
    max_text_width(
      rownames(heat)[seq(1, nrow(heat), 40)],
      gp = gpar(fontsize = 10,  fontface = 'bold')))

# Clusters
# pamClusters <- cluster::pam(heat, k = 4) # pre-select k = 4 centers
pamClusters <- cluster::pam(heat, k= 8) # pre-select k = 4 centers

pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4','Cluster 5', 'Cluster 6'))


hmap <- Heatmap(heat,
                
                # split the genes / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = TRUE,
                
                name = 'Gene\nZ-\nscore',

                col = colorRamp2(myBreaks, myCol),
                
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(8, 'cm'),
                  legend_height = unit(5.0, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 12, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                
                # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                row_title = 'Statistically significant genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                
                # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2',
                
                # specify top and bottom annotations
                top_annotation = colAnn)
                # bottom_annotation = boxplotCol)


pdf("Heatmap2.pdf", height = 18, width = 6)
draw(hmap + genelabels,
     heatmap_legend_side = 'left',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left')
dev.off()
```

# Export

``` r
library(writexl)

# combine results

res <- res_WT_BL_vs.pcry_BL
res1 <- res_WT_D_vs.pcry_D
res2 <- res_WT_BL_vs.pcry_BL
res3 <- res_WT_R_vs.pcry_R
res4 <- res_WT_D_vs.WT_BL
res5 <- res_WT_D_vs.WT_R

colnames(anno)

res_exp <- data.frame(
  # general
  "gene_id" = rownames(dds),
  "gene_name" = anno[rownames(dds),"geneSymbol"],
  "Alias" = anno[rownames(dds),"prev.symbols"],
  "id.symbol" = anno[rownames(dds),"id.symbol"],
 "Description" = anno[rownames(dds),"Description"],
 "Comments" = anno[rownames(dds),"Comments"],
 "TargetP" = anno[rownames(dds),"TargetP"],
 "Predalgo" = anno[rownames(dds),"Predalgo"],
 "baseMean" = res$baseMean,
 # results 1
 "L2FC.WT_D_vs.pcry_D" = res1$log2FoldChange,
 "pvalue.WT_D_vs.pcry_D" = res1$pvalue,
 "padj.WT_D_vs.pcry_D" = res1$padj,
  # results 2
 "L2FC.WT_BL_vs.pcry_BL" = res2$log2FoldChange,
 "pvalue.WT_BL_vs.pcry_BL" = res2$pvalue,
 "padj.WT_BL_vs.pcry_BL" = res2$padj,
  # results 3
 "L2FC.WT_R_vs.pcry_R" = res3$log2FoldChange,
 "pvalue.WT_R_vs.pcry_R" = res3$pvalue,
 "padj.WT_R_vs.pcry_R" = res3$padj,
  # results 4
 "L2FC.WT_D_vs.WT_BL" = res4$log2FoldChange,
 "pvalue.WT_D_vs.WT_BL" = res4$pvalue,
 "padj.WT_D_vs.WT_BL" = res4$padj,
  # results 5
 "L2FC.WT_D_vs.WT_R" = res5$log2FoldChange,
 "pvalue.WT_D_vs.WT_R" = res5$pvalue,
 "padj.WT_D_vs.WT_R" = res5$padj,
  counts(dds, normalized=TRUE))
res_exp

outdir <- "C:/Users/kelterbs/OneDrive - Charité - Universitätsmedizin Berlin/NGS (RNA-Seq & CRISPR)/Chlamy RNA-Seq P3044"

write_xlsx(data.frame(res_exp),
           paste(outdir,"2023_08 P3044 results.xlsx",sep="/"))

# Counts only
counts <- round(counts(dds, normalized=TRUE))
colnames(counts)
colData(dds)$sample.name <- paste(colData(dds)$condition,colData(dds)$repetition,sep="_")
colnames(counts) == colData(dds)$sample.id
colnames(counts) <- colData(dds)$sample.name

write.csv(counts, file.path(outdir,"Counts.csv"))

Exp_design <- data.frame(names = colData(dds)$sample.name,
                         genotype = colData(dds)$genotype,
                         treatment = colData(dds)$treatment,
                         condition = colData(dds)$condition) %>% t()

Exp_design <- matrix(colData(dds)$sample.name,
                     colData(dds)$genotype,
                     colData(dds)$treatment,
                    colData(dds)$condition)

rownames(Exp_design)[1] <- {1}


write.table(Exp_design,file.path(outdir,"Exp_design.csv"), sep=",",  col.names=FALSE)
```

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/intel/oneapi/mkl/2024.0/lib/libmkl_rt.so.2;  LAPACK version 3.10.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8    
    ##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
    ##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] VennDiagram_1.7.3           futile.logger_1.4.3        
    ##  [3] vsn_3.72.0                  ggpubr_0.6.0               
    ##  [5] viridis_0.6.5               viridisLite_0.4.2          
    ##  [7] knitr_1.47                  kableExtra_1.4.0           
    ##  [9] ape_5.8                     biomaRt_2.60.0             
    ## [11] writexl_1.5.0               pheatmap_1.0.12            
    ## [13] EnhancedVolcano_1.22.0      DESeq2_1.44.0              
    ## [15] SummarizedExperiment_1.34.0 Biobase_2.64.0             
    ## [17] MatrixGenerics_1.16.0       matrixStats_1.3.0          
    ## [19] GenomicRanges_1.56.0        GenomeInfoDb_1.40.1        
    ## [21] IRanges_2.38.0              S4Vectors_0.42.0           
    ## [23] AnnotationHub_3.12.0        BiocFileCache_2.12.0       
    ## [25] dbplyr_2.5.0                BiocGenerics_0.50.0        
    ## [27] curl_5.2.1                  tximport_1.32.0            
    ## [29] tximeta_1.22.1              lubridate_1.9.3            
    ## [31] forcats_1.0.0               dplyr_1.1.4                
    ## [33] purrr_1.0.2                 readr_2.1.5                
    ## [35] tidyr_1.3.1                 tibble_3.2.1               
    ## [37] tidyverse_2.0.0             plyr_1.8.9                 
    ## [39] data.table_1.15.4           sessioninfo_1.2.2          
    ## [41] RColorBrewer_1.1-3          R.utils_2.12.3             
    ## [43] R.oo_1.26.0                 R.methodsS3_1.8.2          
    ## [45] stringr_1.5.1               PCAtools_2.16.0            
    ## [47] ggrepel_0.9.5               ggplot2_3.5.1              
    ## [49] patchwork_1.2.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] BiocIO_1.14.0             bitops_1.0-7             
    ##   [3] filelock_1.0.3            preprocessCore_1.66.0    
    ##   [5] XML_3.99-0.16.1           lifecycle_1.0.4          
    ##   [7] httr2_1.0.1               mixsqp_0.3-54            
    ##   [9] rstatix_0.7.2             MASS_7.3-61              
    ##  [11] lattice_0.22-6            ensembldb_2.28.0         
    ##  [13] backports_1.5.0           magrittr_2.0.3           
    ##  [15] limma_3.60.2              rmarkdown_2.27           
    ##  [17] yaml_2.3.8                cowplot_1.1.3            
    ##  [19] DBI_1.2.3                 maps_3.4.2               
    ##  [21] abind_1.4-5               zlibbioc_1.50.0          
    ##  [23] AnnotationFilter_1.28.0   RCurl_1.98-1.14          
    ##  [25] rappdirs_0.3.3            GenomeInfoDbData_1.2.12  
    ##  [27] irlba_2.3.5.1             dqrng_0.4.1              
    ##  [29] svglite_2.1.3             DelayedMatrixStats_1.26.0
    ##  [31] codetools_0.2-20          DelayedArray_0.30.1      
    ##  [33] xml2_1.3.6                tidyselect_1.2.1         
    ##  [35] UCSC.utils_1.0.0          farver_2.1.2             
    ##  [37] ScaledMatrix_1.12.0       ash_1.0-15               
    ##  [39] GenomicAlignments_1.40.0  jsonlite_1.8.8           
    ##  [41] systemfonts_1.1.0         tools_4.4.0              
    ##  [43] progress_1.2.3            Rcpp_1.0.12              
    ##  [45] glue_1.7.0                Rttf2pt1_1.3.12          
    ##  [47] gridExtra_2.3             SparseArray_1.4.8        
    ##  [49] xfun_0.44                 withr_3.0.0              
    ##  [51] formatR_1.14              BiocManager_1.30.23      
    ##  [53] fastmap_1.2.0             fansi_1.0.6              
    ##  [55] truncnorm_1.0-9           digest_0.6.35            
    ##  [57] rsvd_1.0.5                timechange_0.3.0         
    ##  [59] R6_2.5.1                  colorspace_2.1-0         
    ##  [61] RSQLite_2.3.7             hexbin_1.28.3            
    ##  [63] utf8_1.2.4                generics_0.1.3           
    ##  [65] rtracklayer_1.64.0        prettyunits_1.2.0        
    ##  [67] httr_1.4.7                S4Arrays_1.4.1           
    ##  [69] pkgconfig_2.0.3           gtable_0.3.5             
    ##  [71] blob_1.2.4                XVector_0.44.0           
    ##  [73] htmltools_0.5.8.1         carData_3.0-5            
    ##  [75] ProtGenerics_1.36.0       scales_1.3.0             
    ##  [77] png_0.1-8                 ashr_2.2-63              
    ##  [79] lambda.r_1.2.4            rstudioapi_0.16.0        
    ##  [81] tzdb_0.4.0                reshape2_1.4.4           
    ##  [83] rjson_0.2.21              nlme_3.1-165             
    ##  [85] cachem_1.1.0              KernSmooth_2.23-24       
    ##  [87] BiocVersion_3.19.1        parallel_4.4.0           
    ##  [89] extrafont_0.19            AnnotationDbi_1.66.0     
    ##  [91] restfulr_0.0.15           pillar_1.9.0             
    ##  [93] vctrs_0.6.5               BiocSingular_1.20.0      
    ##  [95] car_3.1-2                 beachmat_2.20.0          
    ##  [97] extrafontdb_1.0           evaluate_0.24.0          
    ##  [99] invgamma_1.1              GenomicFeatures_1.56.0   
    ## [101] cli_3.6.2                 locfit_1.5-9.9           
    ## [103] compiler_4.4.0            futile.options_1.0.1     
    ## [105] Rsamtools_2.20.0          rlang_1.1.4              
    ## [107] crayon_1.5.2              SQUAREM_2021.1           
    ## [109] ggsignif_0.6.4            labeling_0.4.3           
    ## [111] affy_1.82.0               stringi_1.8.4            
    ## [113] BiocParallel_1.38.0       ggalt_0.4.0              
    ## [115] txdbmaker_1.0.0           munsell_0.5.1            
    ## [117] Biostrings_2.72.1         lazyeval_0.2.2           
    ## [119] proj4_1.0-14              Matrix_1.7-0             
    ## [121] hms_1.1.3                 sparseMatrixStats_1.16.0 
    ## [123] bit64_4.0.5               KEGGREST_1.44.0          
    ## [125] statmod_1.5.0             highr_0.11               
    ## [127] broom_1.0.6               memoise_2.0.1            
    ## [129] affyio_1.74.0             bit_4.0.5
