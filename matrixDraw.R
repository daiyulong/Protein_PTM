#!/usr/bin/env Rscript
# Title     : draw matrix data
# Objective : matrix file and sample group file
# Created by: yuhong
# Created on: 2020/9/16
###
library(reshape2)
library(ggplot2)
library(hyplot)
library(stringr)

ARGS <- commandArgs(TRUE)

if(length(ARGS) < 6){
  print(ARGS)
  cat("Usage: matrixDraw.R input_matrix_data_file input_sample_group_file outputfold ncol4facet ylab scale [TRUE|FALSE delete min value]")
  q()
}

input <- as.character(ARGS[1])
inputgroup <- as.character(ARGS[2])
pn    <- as.character(ARGS[3])
ncol <-as.integer(ARGS[4])
ylab <-as.character(ARGS[5])
scaleMethod<-as.character(ARGS[6])


if(FALSE){
  input='matrix.txt'
  inputgroup='groupFile.txt'
  pn='test'
  ncol=3
  ylab='Intensity'
  scaleMethod='log2'
  delMin='TRUE'
}

if(length(ARGS)==7){
	delMin <- as.character(ARGS[7])
	if(delMin == 'TRUE'){
		delMin=TRUE
	}else{
		delMin=FALSE
	}
}else{
	delMin<-FALSE
}



if(scaleMethod!='none'){
	ylab=Hmisc::capitalize(paste0(scaleMethod," Scaled ",ylab))
}

ggsci="lancet"

# input='08h_Rep_log2.txt'
# inputgroup='group.txt'
# pn='8hLog2'
# ncol=3
# ylab="Log2 Intensity"


xlab = "Sample"

if(!dir.exists(pn)){
  dir.create(pn)
}


if(str_detect('.csv$',input)){
  df <- read.csv(input,check.names = FALSE)
}else{
  df <- read.delim(input,check.names = FALSE)
}

if(str_detect('.csv$',inputgroup)){
  dg<-read.csv(inputgroup,check.names = FALSE)
}else{
  dg <- read.delim(inputgroup,check.names = FALSE)
}
colnames(dg)[1:2]<- c("Sample","Group")

dg<-dg[,1:2]

df<-df[,c(1,which(colnames(df) %in% dg[[1]]))]

bakdf<-df

minValue<- min(df[,2:ncol(df)])
if(delMin){
  for(i in 2:ncol(df)){
    df[[i]][df[[i]]==minValue]<-NA
  }
}
cn <- colnames(df)

width = (NCOL(df)-1)*0.8+2
height = 6




ldf<-melt(df,id.vars=cn[1],measure.vars=cn[2:length(cn)],variable.name="Sample",value.name="Intensity")
names(ldf)<-c("ID","Sample","Intensity")




d<-merge(ldf, dg, by='Sample', all=TRUE)
colnames(d)
d<-d[!is.na(d$Intensity) & d$Intensity!=Inf & d$Intensity!=-Inf,]


if(delMin){
	d<-d[!is.na(d$Intensity),]
}

if(scaleMethod == "log2"){
	if(any(d$Intensity==0,na.rm=TRUE)){
		cat("Use log2(x+1)\n")
		scaleMethod="log2plus";
		d$Intensity<-log2(d$Intensity+1)
	}else{
		d$Intensity<-log2(d$Intensity)
	}
}else if(scaleMethod == "log2plus"){
	d$Intensity<-log2(d$Intensity+1)
}else if(scaleMethod == "log10"){
	if(any(d$Intensity==0)){
		cat("Use log10(x+1)\n")
		scaleMethod="log10plus";
		d$Intensity<-log10(d$Intensity+1)
	}else{
		d$Intensity<-log10(d$Intensity)
	}
}

dataConvert<-function(d,method){
  multi.zscore<-function(x){
    (x-mean(x))/sd(x)
  }
  
  if(method == 'log2'){
    d<-cbind(d[[1]],log2(d[,2:NCOL(d)]))
  }else if(method == "log2plus"){
    d<-cbind(d[[1]],log2(d[,2:NCOL(d)]+1))
  }else if(method == 'log10'){
    d<-cbind(d[[1]],log10(d[,2:NCOL(d)]))
  }else if(method == 'zscore'){
    d<-cbind(d[[1]], apply(d[,2:NCOL(d)], 1, multi.zscore))
  }else if(method == 'log10plus'){
    d<-cbind(d[[1]], log10(d[,2:NCOL(d)]+1))
  }
  return(d)
}

scaledf<-dataConvert(df,scaleMethod)

#plotBoxplotGroup
plotBoxplotGroup2 <- function(df, plotPoint=TRUE, plotMeanPoint=TRUE, plotBox=TRUE, outimage="", title="", xlab="", ylab="", filllabel="",
                             colorlabel="", palette="Set2", width=8, height=6, units='in', dpi=300, ggsci=NA){
  library(ggplot2)
  library(cowplot)
  library(ggsci)
  mycols<- c(ggsci::pal_lancet("lanonc")(9)[1:8], ggsci::pal_jama("default")(7), ggsci::pal_nejm("default")(8), ggsci::pal_jco("default")(10))
  
  colnames(df)<-c("Sample","Group","Value")
  df$Group<-as.factor(df$Group)
  df$Sample<-as.factor(df$Sample)
  GroupNumber<-length(levels(factor(df$Group)))
  SimpleNumber<-length(unique(df$Sample))
  
  if(SimpleNumber>20){
    width = SimpleNumber*0.5+2
    barwidth=0.7
    pointsize=1
  }else{
    barwidth=0.5
    pointsize=2
  }
  if(width>40){
    width=40
  }
  
  df$Value[is.na(df$Value)]<-0
  #library(ggplot2)
  #library(cowplot,quietly = TRUE)
  #library(ggbeeswarm)
  
  
  if(plotPoint){
    p<-ggplot(df,aes(x=Sample, y=Value))+
      geom_boxplot(fill=NA,width=barwidth,outlier.colour=NA)+
      ggtitle(title)+xlab(xlab)+ylab(ylab)+labs(color=colorlabel)+
      #scale_fill_brewer(palette=palette)+
      theme_cowplot()
    p<-p+geom_quasirandom(aes(colour=Group),width=0.2)
    
    if(is.na(ggsci)){
      p<-p+scale_color_brewer(palette=palette)
    }else{
      p<-addGGsci(p,'color', ggsci)
    }
    
  }else{
    p<-ggplot(df,aes(x=Sample,y=Value,fill=Group))+
      geom_boxplot(width=barwidth,outlier.colour="darkred")+
      ggtitle(title)+xlab(xlab)+ylab(ylab)+labs(fill=filllabel)+
      theme_cowplot()
    if(is.na(ggsci)){
      p<-p+scale_fill_manual(values=mycols[1:GroupNumber])
    }else{
      p<-addGGsci(p,'fill',ggsci)
    }
  }
  if(plotMeanPoint){
    p<-p+ stat_summary(fun=mean, colour="#ECECFF", geom="point", shape=18, size=pointsize, show.legend = FALSE)
  }
  
  p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
  print(p)
  
  if(outimage!=""){
    ggsave(outimage,dpi=dpi,width=width,height = height, units = units)
  }
  return(p)
}

#plotViolinGroup
plotViolinGroup2 <- function(df, plotPoint=FALSE, plotMeanPoint=FALSE, plotBox=TRUE, outimage="", title="", xlab="", ylab="", filllabel="",
                            colorlabel="", palette="Set2", width=8, height=6, units='in', dpi=300, ggsci=NA){
  colnames(df)<-c("Sample","Group","Value")
  df$Value[is.na(df$Value)]<-0
  
  library(ggplot2)
  library(cowplot)
  library(ggsci)
  mycols<- c(ggsci::pal_lancet("lanonc")(9), ggsci::pal_jama("default")(7), ggsci::pal_nejm("default")(8), ggsci::pal_jco("default")(10))
  GroupNumber<-length(levels(factor(df$Group)))
  SimpleNumber<-length(unique(df$Sample))
  
  if(SimpleNumber>20){
    width = SimpleNumber*0.5+2
    barwidth=0.7
    pointsize=1
  }else{
    barwidth=0.5
    pointsize=2
  }
  print(SimpleNumber)
  print(barwidth)
  
  if(width>40){
    width=40
  }
  
  if(plotPoint){
    p<-ggplot(df,aes(x=Sample,y=Value))+
      #geom_violin(width=.8,trim=FALSE)+
      #geom_boxplot(fill=NA,width=0.5,outlier.colour=NA)+
      ggtitle(title)+xlab(xlab)+ylab(ylab)+labs(color=colorlabel)+
      #scale_fill_brewer(palette=palette)+
      theme_cowplot()
    
    p<-p+geom_quasirandom(aes(colour=Group),width=0.2)
    if(is.na(ggsci)){
      p<-p+scale_color_brewer(palette=palette)
    }else{
      p<-addGGsci(p,'color', ggsci)
    }
  }else{
    p<-ggplot(df,aes(x=Sample,y=Value,fill=Group))+
      geom_violin(width=.8,trim=FALSE)+
      #geom_boxplot(width=0.5,outlier.colour="darkred")+
      ggtitle(title)+xlab(xlab)+ylab(ylab)+labs(fill=filllabel)+
      theme_cowplot()
    if(is.na(ggsci)){
      p<-p+scale_fill_manual(values=mycols[1:GroupNumber])
    }else{
      p<-addGGsci(p,'fill',ggsci)
    }
  }
  if(plotBox){
    if(plotPoint){
      p<-p+geom_boxplot(fill=NA,width=0.2,outlier.colour=NA)
    }else{
      p<-p+geom_boxplot(fill='#FCFCFC',width=0.2,outlier.colour=NA)
    }
    
  }
  if(plotMeanPoint){
    p<-p+ stat_summary(fun=mean, colour="#ECECFF", geom="point", shape=18, size=pointsize, show.legend = FALSE)
  }
  if(length(df$Group[!duplicated(df$Group)]) == 1){
    p<-p + theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1))
  }else{
    p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
  }
  
  print(p)
  
  if(outimage!=""){
    ggsave(outimage,dpi=dpi,width=width,height = height, units = units)
  }
  return(p)
}

#boxplot 绠辩嚎鍥?
plotBoxplotGroup2(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/boxplot.png"),
                 plotPoint=FALSE, width=width, height=4, ylab=ylab)
plotBoxplotGroup2(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/boxplot.pdf"),
                 plotPoint=FALSE, width=width, height=4, ylab=ylab)

#violin 灏忔彁鐞村浘
plotViolinGroup2(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/violin.png"),plotMeanPoint=TRUE,
                plotPoint=FALSE, width=width, height=4, ylab=ylab)
plotViolinGroup2(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/violin.pdf"),plotMeanPoint=TRUE,
                plotPoint=FALSE, width=width, height=4, ylab=ylab)

plotHistogramFacet <- function(df, outimage="", color="grey50",ncol=2,xlab="",ylab="",
                               width=4, height=3, units='in', dpi=300, ggsci=NA,bins=50, barcolor=NA){
  
  colnames(df)<-c("Sample","Group","Value")
  

  
  sample <- df$Sample[!duplicated(df$Sample)]
  sampleNumb<- length(as.character(unique(sample)))
  width <- width * ncol
  height <- height * ceiling(sampleNumb/ncol)
  
  df$Value[is.na(df$Value)]<-0
  
  library(ggplot2)
  library(cowplot)
  library(ggsci)
  mycols<- c(ggsci::pal_lancet("lanonc")(9), ggsci::pal_jama("default")(7), ggsci::pal_nejm("default")(8), ggsci::pal_jco("default")(10))
  GroupNumber<-length(levels(factor(df$Group)))
  SimpleNumber<-length(df$Sample)
  
  p<-ggplot(df,aes(x=Value,fill=Group))
  p<-p+geom_histogram(bins=bins, colour=barcolor)+
    geom_rug(alpha=.2,colour=color)+
    facet_wrap(~Sample,ncol=ncol)+
    xlab(xlab)+ylab(ylab)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank())
  
  if(is.na(ggsci)){
    p<-p+scale_fill_manual(values=mycols[1:GroupNumber])
  }else{
    p<-addGGsci(p,'fill',ggsci)
  }
  
  
  #p <- p+theme_cowplot()
  
  print(p)
  if(outimage!=""){
    ggsave(outimage,dpi=dpi,width=width,height = height, units = units)
  }
  return(p)
}



if(TRUE){
#histogram 鐩存柟鍥?
plotHistogramFacet(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/histogramFacet.png"),
                   xlab=ylab, ylab="Count",ncol=ncol,width=4,height=3,bins=50,barcolor='white')
plotHistogramFacet(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/histogramFacet.pdf"),
                   xlab=ylab, ylab="Count",ncol=ncol,width=4,height=3,bins=50,barcolor='white')

#density 瀵嗗害鍥?
plotDensityFacet(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/densityFacet.png"),
                 xlab=ylab, ylab="Density", ncol=ncol,width=4,height=3)
plotDensityFacet(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/densityFacet.pdf"),
                 xlab=ylab, ylab="Density", ncol=ncol,width=4,height=3)

#density 
plotDensity(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/density.png"),
            xlab=ylab, ylab="Density", width=8, height=6)
plotDensity(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/density.pdf"),
            xlab=ylab, ylab="Density", width=8, height=6)
}
#Correlation coefficient matrix plot 
#mdf<-df[2:ncol(df)]
#rownames(mdf)<-as.character(df[[1]])
#q()


mdf <- scaledf[2:ncol(scaledf)]
rownames(mdf)<-as.character(scaledf[[1]])
m<-as.matrix(mdf)
#rownames(mdf)<-as.vector(df[,1])
mcor <- cor(m,use='complete.obs')
#library(ggcorrplot)
#ggcorrplot(mcor, method='circle', type='upper',show.diag = FALSE, lab=TRUE)
#ggsave(paste0(pn,"/corr.png"), width=7, height=6, units='in', dpi=300)
#ggsave(paste0(pn,"/corr.pdf"), width=7, height=6, units='in', dpi=300)

SampleNumber<-length(unique(d$Sample))
if(SampleNumber>30){
  tl.cex=0.5
}else if(SampleNumber>20){
  tl.cex=0.6
}else if(SampleNumber>15){
  tl.cex=0.8
}else{
  tl.cex=1
}

library(corrplot)
if(SampleNumber<15){
  palette_1 <- RColorBrewer::brewer.pal(n=11, name = "RdYlGn") 
  palette_2 <- rev(palette_1)
  png(paste0(pn,"/corr.png"), width=7, height=7, units='in', res=300)
  corrplot(mcor,
           method = "circle",
           number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           addCoef.col = "blue",
           number.digits=2,
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  cairo_pdf(paste0(pn,"/corr.pdf"), width=7, height=7)
  corrplot(mcor,
           method = "circle",
           number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           addCoef.col = "blue",
           number.digits=2,
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  
  
}else{
  palette_1 <- RColorBrewer::brewer.pal(n=11, name = "RdYlGn") 
  palette_2 <- rev(palette_1)
  png(paste0(pn,"/corr.png"), width=7, height=7, units='in', res=300)
  corrplot(mcor,
           method = "pie",
           #number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           #addCoef.col = "blue",
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  cairo_pdf(paste0(pn,"/corr.pdf"), width=7, height=7)
  corrplot(mcor,
           method = "pie",
           #number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           #addCoef.col = "blue",
           
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  
}

#ggpairs 鐭╅樀鏁ｇ偣鍥?
#if(NCOL(df)<=20){
# library(ggcor)
# quickcor(m, cor.test = TRUE) +
#   geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
#   geom_mark(data = get_data(type = "upper", show.diag = FALSE), size = 2.5) +
#   geom_abline(slope = -1, intercept = 12)
# quickcor(m, type = "upper") + geom_circle2()
library(GGally, quietly = TRUE)
ggpairs(mdf)
p<-ggpairs(mdf, axisLabels="none")
ggsave(paste0(pn,"/ggpairs.png"), plot=p, width=(NCOL(df)-1)*1.5, height=(NCOL(df)-1)*1.5, units='in', dpi=300,limitsize = FALSE)
ggsave(paste0(pn,"/ggpairs.pdf"), plot=p, width=(NCOL(df)-1)*1.5, height=(NCOL(df)-1)*1.5, units='in', dpi=300,limitsize = FALSE)
#}


# 无scale的相关性图
mdf <- df[2:ncol(df)]
rownames(mdf)<-as.character(df[[1]])
m<-as.matrix(mdf)
outfold<-paste0(pn, '/../images/NonScale')
if(!dir.exists(outfold)){
  dir.create(outfold, recursive = TRUE)
}
#rownames(mdf)<-as.vector(df[,1])
write.csv(m, paste0(outfold,"/cor_input.csv"))
mcor <- cor(m, use='complete.obs')
write.csv(mcor,paste0(outfold,'/cor.csv'))
#library(ggcorrplot)
#ggcorrplot(mcor, method='circle', type='upper',show.diag = FALSE, lab=TRUE)
#ggsave(paste0(pn,"/corr.png"), width=7, height=6, units='in', dpi=300)
#ggsave(paste0(pn,"/corr.pdf"), width=7, height=6, units='in', dpi=300)

SampleNumber<-length(unique(d$Sample))
if(SampleNumber>30){
  tl.cex=0.5
}else if(SampleNumber>20){
  tl.cex=0.6
}else if(SampleNumber>15){
  tl.cex=0.8
}else{
  tl.cex=1
}

library(corrplot)
if(SampleNumber<15){
  palette_1 <- RColorBrewer::brewer.pal(n=11, name = "RdYlGn") 
  palette_2 <- rev(palette_1)
  png(paste0(outfold,"/corr.png"), width=7, height=7, units='in', res=300)
  corrplot(mcor,
           method = "circle",
           number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           addCoef.col = "blue",
           number.digits=2,
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  cairo_pdf(paste0(outfold,"/corr.pdf"), width=7, height=7)
  corrplot(mcor,
           method = "circle",
           number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           addCoef.col = "blue",
           number.digits=2,
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  
  
}else{
  palette_1 <- RColorBrewer::brewer.pal(n=11, name = "RdYlGn") 
  palette_2 <- rev(palette_1)
  png(paste0(outfold,"/corr.png"), width=7, height=7, units='in', res=300)
  corrplot(mcor,
           method = "pie",
           #number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           #addCoef.col = "blue",
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  cairo_pdf(paste0(outfold,"/corr.pdf"), width=7, height=7)
  corrplot(mcor,
           method = "pie",
           #number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           #addCoef.col = "blue",
           
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  
}
library(GGally, quietly = TRUE)
ggpairs(mdf)
p<-ggpairs(mdf, axisLabels="none")
ggsave(paste0(outfold,"/ggpairs.png"), plot=p, width=(NCOL(df)-1)*1.5, height=(NCOL(df)-1)*1.5, units='in', dpi=300,limitsize = FALSE)
ggsave(paste0(outfold,"/ggpairs.pdf"), plot=p, width=(NCOL(df)-1)*1.5, height=(NCOL(df)-1)*1.5, units='in', dpi=300,limitsize = FALSE)







checkValue<-function(x){
  if(sum(x==minValue) > length(x)*0.8){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

#heatmap 鐑浘
library(pheatmap)
annotation_col <- dg
rownames(annotation_col)<-dg[[1]]
annotation_col$Sample<-NULL

mmdf<-bakdf[,2:ncol(bakdf)]
rownames(mmdf)<-as.character(bakdf[[1]])

if(delMin){
  flag <- apply(mmdf, 1, checkValue)
  mdf<-mmdf[flag,]
}else{
  mdf<-mmdf
}
border=NA

#write.csv(mdf, "test.csv")
try(
pheatmap(mdf,
         filename=paste0(pn,"/heatmap.png"),
         width=4,
         height=6,
         scale = 'row',
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         annotation_col=annotation_col,
         #color =colorRampPalette(c("blue","white","red"))(51),
         color =colorRampPalette(c("green","black","red"))(51),
         
         border=border,
         clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
         clustering_distance_cols  = "correlation",#"euclidean",#"correlation",
         clustering_method = "complete",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
         #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
         show_rownames = FALSE)
)


if(NROW(mdf)>500){
  try(
  pheatmap(mdf,
           filename=paste0(pn,"/heatmap.pdf"),
           cellwidth=10,
           cellheight=1,
           scale = 'row',
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           annotation_col=annotation_col,
           #color =colorRampPalette(c("blue","white","red"))(51),
           color =colorRampPalette(c("green","black","red"))(51),

           border=border,
           clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
           clustering_distance_cols  = "euclidean",#"correlation",#"correlation",
           clustering_method = "average",#"complete",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
           #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
           show_rownames = FALSE)
  )
}else{
  try(
  pheatmap(mdf,
           filename=paste0(pn,"/heatmap.pdf"),
           cellwidth=10,
           height=8,
           scale = 'row',
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           annotation_col=annotation_col,
           #color =colorRampPalette(c("blue","white","red"))(51),
           color =colorRampPalette(c("green","black","red"))(51),

           border=border,
           clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
           clustering_distance_cols  = "euclidean",#"correlation",#"correlation",
           clustering_method = "average",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
           #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
           show_rownames = FALSE)
  )
}

# 
# 
# 
# 
# 
# doPCA<-function(mydata,pn,dg,retx=TRUE,center=TRUE,scale=TRUE,width=8,height=6,units="in",dpi=300){
#   #make output fold
# 
#   library(ggrepel)
#   library(cowplot)
#   library(ggsci)
#   if(!dir.exists(pn)) dir.create(pn)
# 
#   #pca
#   mydata.pca<-prcomp(mydata,retx=retx,center=center,scale.=scale)
#   z<-summary(mydata.pca)
#   z<-format(z$importance[2,]*100,digits=2)
#   sd<-mydata.pca$sdev
#   loadings<-mydata.pca$rotation
#   rownames(loadings)<-colnames(mydata)
#   scores<-mydata.pca$x
#   
#   groups<-dg[rownames(mydata),][[2]]
#   d<-data.frame(scores,Samples=rownames(mydata),SD=sd,CP=z,Group=groups)
# 
#   #biplot
#   d$Group<-as.factor(d$Group)
#   png(paste(pn,"/biplot.png",sep=""),res=dpi,width=width,height = height, units = units)
#   biplot(mydata.pca,cex=.1)
#   dev.off()
# 
#   p<-ggplot(d,aes(x=PC1,y=PC2,colour=Group))+
#     geom_point(size=5)+
#     geom_text_repel(label=d$Sample,size=3, colour='black')+
#     xlab(paste("PC1(",z[1],"%)",sep=""))+
#     ylab(paste("PC2(",z[2],"%)",sep=""))+
#     stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.1,color=NA)+
#     scale_color_lancet()+
#     scale_fill_lancet()+
#     guides(fill="none")+
#     theme_bw()
#     
#   ggsave(paste(pn,"/PCA1.2.png",sep=""),plot=p,dpi=dpi,width=width,height = height, units =units)
#   ggsave(paste(pn,"/PCA1.2.pdf",sep=""),plot=p,dpi=dpi,width=width,height = height, units =units)
#   #png(paste("PCA/",pn,".PCA1.2.3.png",sep=""),res=dpi,width=height,height = height, units = units)
#   
#   p<-ggplot(d,aes(x=PC1,y=PC2))+
#     geom_point(aes(colour=Group),size=5)+
#     geom_text_repel(label=d$Sample,size=3, colour='black')+
#     xlab(paste("PC1(",z[1],"%)",sep=""))+
#     ylab(paste("PC2(",z[2],"%)",sep=""))+
#     stat_ellipse(type="norm",geom="polygon",alpha=0.1,color=NA,fill="grey")+
#     scale_color_lancet()+
#     scale_fill_lancet()+
#     guides(fill="none")+
#     theme_bw()
#   
#   ggsave(paste(pn,"/PCA1.2_style3.png",sep=""),plot=p,dpi=dpi,width=width,height = height, units =units)
#   ggsave(paste(pn,"/PCA1.2_style3.pdf",sep=""),plot=p,dpi=dpi,width=width,height = height, units =units)
#   
#   # library(ggpubr)
#   # p<-ggscatter(d, x = "PC1", y = "PC2",
#   #           color = "Group", shape = "Group", label = d$Sample,
#   #           palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#   #           xlab=paste("PC1(",z[1],"%)",sep=""),
#   #           ylab=paste("PC2(",z[2],"%)",sep=""),
#   #           ellipse = TRUE, mean.point = TRUE,
#   #           star.plot = TRUE)
#   # ggsave(paste(pn,"/PCA1.2_style2.png",sep=""),plot=p,dpi=dpi,width=width,height = height, units =units)
#   # ggsave(paste(pn,"/PCA1.2_style2.pdf",sep=""),plot=p,dpi=dpi,width=width,height = height, units =units)
# 
#   #涓夌淮鏁ｇ偣鍥?
#   library(scatterplot3d)
#   png(paste(pn,"/PCA123.png",sep=""),res=dpi,width=width,height = height, units = units)
#   colors<-colorDefine(d$Group)
#   s3d<-scatterplot3d(scores[,1:3],xlab=paste0("PC1(",z[1],"%)"),ylab=paste0("PC2(",z[2],"%)"),zlab=paste0("PC3(",z[3],"%)"),
#                      color=colors,
#                      pch=16,
#                      type="h")
#   text(s3d$xyz.convert(scores[,1:3]),labels=rownames(scores),cex= 0.6)
#   dev.off()
#   # library(gg3D)
#   # theta=45
#   # phi=0
#   # #theta鎺у埗鏃嬭浆锛宲hi鎺у埗鍥剧殑鍊炬枩
#   # df3d<-as.data.frame(scores[,1:3])
#   # df3d$Group<-as.character(dg[[2]])
#   # ggplot(df3d, aes(x=PC1,y=PC2, z=PC3,colour=Group))+
#   #   axes_3D(theta=theta, phi=phi) + stat_3D(theta=theta, phi=phi) +
#   #   axis_labs_3D(theta=theta, phi=phi, size=3, hjust=c(1,1,1.2,1.2,1.2,1.2),
#   #                vjust=c(-.5,-.5,-.2,-.2,1.2,1.2)) +
#   #   labs_3D(theta=theta, phi=phi, hjust=c(1,0,0), vjust=c(1.5,1,-.2),
#   #           labs=c("PC1", "PC2", "PC3")) +theme_void()
#   #install.packages('devtools')
#   #devtools::install_github("AckerDWM/gg3D")
# }
# 
# 
# mydata<-t(mdf)
# rownames(dg)<-dg[[1]]
# doPCA(mydata,pn,dg)

#PCA 涓绘垚鍒嗗垎鏋?
#install:
#library(devtools)
#install_github("vqv/ggbiplot")
#library(ggbiplot)
#
# mydata.pca <- prcomp(mydata, scale. = TRUE)
# #----纰庣煶鍥?
# p<-ggbiplot(mydata.pca, obs.scale = 1,
#          choices=1:2,
#          var.scale = 1,
#          groups = as.character(dg[[2]]),
#          pc.biplot =FALSE,
#          ellipse = TRUE,
#          circle = TRUE) +
#   scale_color_discrete(name = '') +
#   theme(legend.direction = 'horizontal', legend.position = 'top')
#
# #----鏁ｇ偣鍥?
# #      PC1 and PC2
# pcad<-as.data.frame(mdf.pca)
#

