#!/usr/bin/env Rscript

# 差异蛋白绘图
#20220921
# 修改group获取方式

suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(hyplot))

ARGS<-commandArgs(TRUE)

if(length(ARGS)<3){
  cat("Usage: Rscript dDEP.R DEP_input.txt group.txt pn [grouporder]")
  q()
}

inputFile<-as.character(ARGS[1])
groupFile<-as.character(ARGS[2])
pn<-as.character(ARGS[3])

ExpCol=5

#限制最多绘图的个数
maxFigure=300

if(FALSE){
  pn<-'test'
  inputFile<-'DEP.xls'
  groupFile<-'../group.txt'
}

gdd<-read.delim(groupFile)
rownames(gdd)<-gdd[[1]]

data<-read.delim(inputFile, check.names = FALSE)
#group<-gdd[colnames(data)[ExpCol:(NROW(gdd)+ExpCol-1)],2]
group<-gdd[[2]]
if(length(ARGS)==4){
  grouporder <- as.character(ARGS[4])
  groupstring <- strsplit(grouporder, ',')[[1]]
  group<-factor(group, levels=groupstring)
}


rownames(data)<-data[[1]]
colnames(data)[2]<-"Symbol"


dJitterPlot<-function(dd,output,n, width=4,height=4,ylab="MS Intensity"){
  colnames(dd)<-c("Group","Value")
  p<-ggplot(dd,aes(x=Group,y=Value, colour=Group))+
    
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 colour="grey60",
                 width = 0.75, size=1.5)+
    geom_jitter(aes(shape=Group,size=Group),width=.3, #size=3,shape=18,
                alpha=.9)+
    scale_shape_manual(values=c(18,17))+
    scale_size_manual(values=c(3,2.4))+
    xlab("")+theme_classic()+ggtitle(n)+ylab(ylab)+
    scale_color_lancet()+theme(plot.title = element_text(hjust = 0.5))
  ggsave(output, plot=p, width=width,height=height, dpi=300, units='in')
}

dBoxPlot<-function(dd,output,n,width=4,height=4,ylab="MS Intensity"){
  colnames(dd)<-c("Group","Value")
  p<-ggplot(dd,aes(x=Group,y=Value, colour=Group))+
    geom_boxplot(fill=NA, 
                 colour="black",
                 width=.6)+
    geom_jitter(aes(shape=Group,size=Group),width=.28,#size=3,
                #shape=18,
                alpha=.9)+
    scale_shape_manual(values=c(18,17))+
    scale_size_manual(values=c(3,2.4))+
    #stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
    #             colour="black",
    #             width = 0.75, size=1.5)+
    stat_summary(fun = mean, geom = "point",
                 colour="#006400",alpha=0.9,
                 size=3)+
    xlab("")+theme_classic()+ggtitle(n)+ylab(ylab)+
    scale_color_lancet()+theme(plot.title = element_text(hjust = 0.5))
  ggsave(output, plot=p, width=width,height=height, dpi=300, units='in')
}

# dggstat<-function(dd, output, n, width=6, height=4.5){
#   colnames(dd)<-c("Group","Value")
#   p<-ggstatsplot::ggbetweenstats(
#     data = dd,
#     x = Group,
#     y = Value,
#     var.equal=TRUE,
#     title = n
#   )
#   ggsave(output, plot=p, width=width,height=height, dpi=300, units='in')
# }

dScatterSE<-function(dd, output, n, width=4, height=4, ylab="MS Intensity"){
  colnames(dd)<-c("Group","Value")
  gdf<-hyplot::summarySE(dd, measurevar='Value',groupvars=c('Group'))
  p<-ggplot(dd,aes(x=Group,y=Value, colour=Group))+
    geom_jitter(aes(shape=Group,size=Group),width=.28,alpha=.9)+
    geom_errorbar(data=gdf,aes(ymin=Value-se, ymax=Value+se), color="black",width=0.6)+
    geom_point(data=gdf, aes(x=Group, y=Value), size=3,color="#006400", alpha=0.9)+
    scale_shape_manual(values=c(18,17))+
    scale_size_manual(values=c(3,2.4))+
    xlab("")+theme_classic()+ggtitle(n)+ylab("MS Intensity")+
    scale_color_lancet()+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(output, plot=p, width=width,height=height, dpi=300, units='in')
}

dBarSE<-function(dd,output,n, width=4,height=4, ylab="MS Intensity"){
  colnames(dd)<-c("Group","Value")
  Group<-as.factor(dd$Group)
  Numb<- NROW(levels(Group))
  gdf<-hyplot::summarySE(dd, measurevar='Value',groupvars=c('Group'))
  p<-ggplot(gdf, aes(x=Group, y=Value, fill=Group))+
    geom_bar(stat="identity",width=0.6)+
    geom_errorbar(aes(ymin=Value-se, ymax=Value+se), color=pal_lancet("lanonc")(9)[1:Numb],width=0.4,size=2)+
    xlab("")+theme_classic()+ggtitle(n)+ylab(ylab)+
    scale_y_continuous(expand=c(0,0),limits=c(0,max(dd$Value)*1.2))+
    ggtitle(n)+
    scale_fill_lancet()+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(output, plot=p, width=4,height=4, dpi=300, units='in')
}

dCompBar<-function(d,output,n,width=4, height=4, ylab="MS Intensity"){
  library(ggpubr)
  colnames(d)<-c("Group","Value")
  
  lg<- levels(factor(d[[1]]))
  complist=list()
  c=1
  for(i in 1:NROW(lg)){
    for(j in i:NROW(lg)){
      if(i!=j){
        complist[[c]]=c(lg[i],lg[j])
        c=c+1
      }
    }
  }
  
  
  ddd<-hyplot::summarySE(d, measurevar="Value",groupvars=c("Group"))
  p<-ggbarplot(d, x="Group",y="Value",
               fill="Group", palette="lancet",
               width=0.5,lab.size = 3,
               add="mean_se",xlab="", ylab=ylab, error.plot = "upper_errorbar",
               title=n)
  
  #my_comp <-list(c("HFB","HFBL"),c("HFB","HFBM"),c("HFB","NC"),c("HFBL","HFBM"),c("HFBL","NC"),c("HFBM","NC"))
  
  p+stat_compare_means(method = "t.test",comparisons = complist,label = "p.signif",hide.ns=TRUE)+
    stat_compare_means(method = "t.test",label.y=max(ddd$Value+ddd$se)*1.5)
  
  ggsave(output,width=width,height=height, units="in",dpi=300)
  
}

mf<-min(maxFigure,NROW(data))

outfold<-paste0(pn,'/Jitter')
if(!dir.exists(outfold)){
  dir.create(outfold, recursive = TRUE)
}

for(i in 1:mf){
  dd<-data.frame(Group=group, Value=as.numeric(data[i,ExpCol:(NROW(gdd)+ExpCol-1)]))
  accession<-rownames(data)[i]
  if(is.na(data[[2]][i]) | data[[2]][i]==""){
    Title <- accession
  }else{
    Title <- paste0(accession,"(",data[[2]][i], ")")
  }
  dJitterPlot(dd, paste0(outfold,'/',i,".pdf"), Title)
  dJitterPlot(dd, paste0(outfold,'/',i,".png"), Title)
}
cat("Draw Jitters finished\n")

outfold<-paste0(pn,'/Boxplot')
if(!dir.exists(outfold)){
  dir.create(outfold, recursive = TRUE)
}

for(i in 1:mf){
  dd<-data.frame(Group=group, Value=as.numeric(data[i,ExpCol:(NROW(gdd)+ExpCol-1)]))
  accession<-rownames(data)[i]
  if(is.na(data[[2]][i]) | data[[2]][i]==""){
    Title <- accession
  }else{
    Title <- paste0(accession,"(",data[[2]][i], ")")
  }
  #dBoxPlot(dd, paste0(outfold,'/',i,".",accession,".pdf"), Title)
  #dBoxPlot(dd, paste0(outfold,'/',i,".",accession,".png"), Title)
  
  dBoxPlot(dd, paste0(outfold,'/',i,".pdf"), Title)
  dBoxPlot(dd, paste0(outfold,'/',i,".png"), Title)
}
cat("Draw Boxplox finished\n")

#######################################################################
# outfold<-paste0(pn,'/DEP/ggstat')
# if(!dir.exists(outfold)){
#   dir.create(outfold)
# }
# 
# for(i in 1:NROW(data)){
#   dd<-data.frame(Group=gdd[[2]], Value=as.numeric(data[i,1:NROW(gdd)]))
#   accession<-rownames(data)[i]
#   if(is.na(data[[NCOL(data)]][i]) | data[[NCOL(data)]][i]==""){
#     Title <- accession
#   }else{
#     Title <- paste0(accession,"(",data[[NCOL(data)]][i], ")")
#   }
#   dggstat(dd, paste0(outfold,'/',i,".",accession,".pdf"), Title)
#   dggstat(dd, paste0(outfold,'/',i,".",accession,".png"), Title)
# }


outfold<-paste0(pn,'/ScatterSE')
if(!dir.exists(outfold)){
  dir.create(outfold, recursive = TRUE)
}

for(i in 1:mf){
  dd<-data.frame(Group=group, Value=as.numeric(data[i,ExpCol:(NROW(gdd)+ExpCol-1)]))
  accession<-rownames(data)[i]
  if(is.na(data[[2]][i]) | data[[2]][i]==""){
    Title <- accession
  }else{
    Title <- paste0(accession,"(",data[[2]][i], ")")
  }

  dScatterSE(dd, paste0(outfold,'/',i,".pdf"), Title)
  dScatterSE(dd, paste0(outfold,'/',i,".png"), Title)
}
cat("Draw ScatterSE finished\n")


########################################################
# outfold<-paste0(pn,'/BarSE')
# if(!dir.exists(outfold)){
#   dir.create(outfold, recursive = TRUE)
# }
# 
# for(i in 1:NROW(data)){
#   dd<-data.frame(Group=gdd[[2]], Value=as.numeric(data[i,ExpCol:(NROW(gdd)+ExpCol-1)]))
#   accession<-rownames(data)[i]
#   if(is.na(data[[2]][i]) | data[[2]][i]==""){
#     Title <- accession
#   }else{
#     Title <- paste0(accession,"(",data[[2]][i], ")")
#   }
#   #dBarSE(dd, paste0(outfold,'/',i,".",accession,".pdf"), Title)
#   #dBarSE(dd, paste0(outfold,'/',i,".",accession,".png"), Title)
#   
#   dBarSE(dd, paste0(outfold,'/',i,".pdf"), Title)
#   dBarSE(dd, paste0(outfold,'/',i,".png"), Title)
# }


outfold<-paste0(pn,'/BarComp')
if(!dir.exists(outfold)){
  dir.create(outfold, recursive = TRUE)
}

for(i in 1:mf){
  dd<-data.frame(Group=group, Value=as.numeric(data[i,ExpCol:(NROW(gdd)+ExpCol-1)]))
  accession<-rownames(data)[i]
  if(is.na(data[[2]][i]) | data[[2]][i]==""){
    Title <- accession
  }else{
    Title <- paste0(accession,"(",data[[2]][i], ")")
  }

  dCompBar(dd, paste0(outfold,'/',i,".pdf"), Title)
  dCompBar(dd, paste0(outfold,'/',i,".png"), Title)
}
cat("Draw barComp finished\n")

