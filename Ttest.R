#!/usr/bin/env Rscript

ARGS<-commandArgs(TRUE)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(hyplot)
library(stringr)
library(xlsx)
################### two sample test ######################
if(length(ARGS)<6){
  cat("Rscript Ttest.R input outputPrefix groupFile method FoldChangeCutoff pvalueCutoff paired var.equal\n")
  q()
}

#input="08h_Rep.txt"
#pn="08hDEP_analysis"
input <- as.character(ARGS[1])
pn    <- as.character(ARGS[2])
ginput <-as.character(ARGS[3])
method <-as.character(ARGS[4])
FCcutoff<-as.numeric(ARGS[5])
pcutoff <-as.numeric(ARGS[6])

paired=FALSE

if(length(ARGS)>6){
  if(ARGS[7] == 'TRUE'){
    paired=TRUE
  }
}

if(FALSE){
  input="test/matrix.txt"
  pn="test/DEP"
  ginput="group.txt"
  method="log2"#"none"
  FCcutoff=2
  pcutoff=0.05
}

# FC = 1 team / 0 team
#cl<-c(0,0,0,1,1,1) # control:0 , Tumer:1 ,  FC  1/0
#method="log2"
#cl<-as.integer(ARGS[3:length(ARGS)])
#maintitle="Volcano graph"
#out=paste(pn,"/ttest.csv",sep=".")

# if(length(ARGS)==3){
# 	FCcutoff=as.numeric(ARGS[3])
# }else{
# 	FCcutoff=1.5
# }
# if(length(ARGS)==4){
#   pcutoff=as.numeric(ARGS[4])
# }else{
#   pcutoff=0.05
# }

width=13.2
height=10.5
dpi=600
units = "cm"

if(!dir.exists(pn)){
  dir.create(pn)
}

getPostfix<-function(s){
  v<-str_split(s,"\\.")
  return(v[length(v)])
}

loadData<-function(input){
  if(getPostfix(input)=="csv"){
    df<-read.csv(input, row.names=1, stringsAsFactors = FALSE, check.names =FALSE)
  }else if(getPostfix(input)=="xlsx"){
    df<-xlsx::read.xlsx(input, sheetIndex = 1)
    rownames(df)<-df[[1]]
    df<-df[,-1]
  }else{
    df<-read.delim(input,row.names=1,  sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
  }
  return(df)
}

df<-loadData(input)

#df<-read.table(input,sep="\t",row.names = 1,header=T)
gd <- loadData(ginput)
cl <- gd[[2]]

m<-as.matrix(df)
mbak<-as.matrix(df)

if(method == "log2"){
  m<-log2(m)
  transm<-log2(df)
  colnames(transm)<-paste0("Log2", colnames(transm))
  df<-cbind(df,transm)
}else if(method =="log2plus1"){
  m<-log2(m+1)
  transm<-log2(df+1)
  colnames(transm)<-paste0("Log2", colnames(transm))
  df<-cbind(df,transm)
}

cv<-function(x){
  100*sd(x)/mean(x)
}

se<-function(x){
  sd(x)/sqrt(length(x))
}

multi.ttest<-function(x){
  errorflag<-try(t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE),silent=TRUE)
  if('try-error' %in% class(errorflag)) NA
  else t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE)$p.value
}

multi.pairedttest<-function(x){
  errorflag<-try(t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE, paired=TRUE),silent=TRUE)
  if('try-error' %in% class(errorflag)) NA
  else t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE, paired=TRUE)$p.value
}

multi.wilcox.test<-function(x){
  errorflag<-try(wilcox.test(x[which(cl==0)],x[which(cl==1)],exact=FALSE,correct=FALSE),silent=TRUE)
  if('try-error' %in% class(errorflag)) NA
  else wilcox.test(x[which(cl==0)],x[which(cl==1)],exact=FALSE,correct=FALSE)$p.value
}

multi.var.test<-function(x){
  var.test(x[which(cl==0)],x[which(cl==1)])$p.value
}

multi.fc<-function(x){
  mean(x[which(cl==1)]) / mean(x[which(cl==0)])
}

if(paired==TRUE){
  df$ttestPvalue<-apply(m,1,multi.pairedttest)
}else{
  df$ttestPvalue<-apply(m,1,multi.ttest)
}

df$ControlMean<-apply(m[,cl==0], 1, mean)
df$ControlSD  <-apply(m[,cl==0], 1, sd)
df$ControlSE  <-apply(m[,cl==0], 1, se)
df$ControlCV  <-apply(m[,cl==0], 1, cv)

df$TestMean   <-apply(m[,cl==1], 1, mean)
df$TestSD     <-apply(m[,cl==1], 1, sd)
df$TestSE     <-apply(m[,cl==1], 1, se)
df$TestCV     <-apply(m[,cl==1], 1, cv)


#apply(m,1,multi.wilcox.test)->df$wilcoxtestPvalue
#apply(m,1,multi.var.test)->df$FtestPvalue

df$FC <- apply(mbak,1,multi.fc)
df$Log2FC<-log2(df$FC)
df$fdr<-p.adjust(df$ttestPvalue,method="fdr")



mydf<-df

df$threshold<-"NoSig"
df$threshold[which(df$Log2FC>=log2(FCcutoff) & df$ttestPvalue <=pcutoff)]<-"Up"
df$threshold[which(df$Log2FC<=log2(1/FCcutoff) & df$ttestPvalue <=pcutoff)]<-"Down"

outdf <- cbind(rownames(df), df)
colnames(outdf)[1] <- "Accession"

write.table(outdf,file=paste0(pn,"/ttest.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(outdf[outdf$threshold=="Up" | outdf$threshold=="Down",],file=paste0(pn,"/DEP.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(outdf[outdf$threshold=="Up",],file=paste0(pn,"/Up.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(outdf[outdf$threshold=="Down",],file=paste0(pn,"/Down.txt"), sep='\t', row.names=FALSE, quote=FALSE)

#write.xlsx(outdf,file=paste0(pn,"/ttest.xlsx"),row.names=FALSE)
#write.xlsx(outdf[outdf$threshold=="Up" | outdf$threshold=="Down",],file=paste(pn,"/DEP.xlsx",sep="."), row.names=FALSE)
#write.xlsx(outdf[outdf$threshold=="Up",],file=paste(pn,"/Up-regulatedDEP.xlsx",sep="."), row.names=FALSE)
#write.xlsx(outdf[outdf$threshold=="Down",],file=paste(pn,"/Down-regulatedDEP.xlsx",sep="."), row.names=FALSE)

# 
# mydf$threshold<-"NoSig"
# mydf$threshold[which(mydf$FC > FCcutoff & mydf$ttestPvalue <0.05)]<-"Up-P0.05"
# mydf$threshold[which(mydf$FC < 1/FCcutoff & mydf$ttestPvalue <0.05)]<-"Down-P0.05"
# mydf$threshold[which(mydf$FC > FCcutoff & mydf$ttestPvalue <0.01)]<-"Up-P0.01"
# mydf$threshold[which(mydf$FC < 1/FCcutoff & mydf$ttestPvalue <0.01)]<-"Down-P0.01"
# 
# table(mydf$threshold)
# 
# colours<-c("Up-P0.01"="#A50026","Up-P0.05"="#D73027","Down-P0.01"="#006837","Down-P0.05"="#1A9850","NoSig"="#CCCCCC")
# shapes<-c("Up-P0.01"=18,"Up-P0.05"=17,"Down-P0.01"=18,"Down-P0.05"=17,"NoSig"=20)
# sizes<-c("Up-P0.01"=1,"Up-P0.05"=1,"Down-P0.01"=1,"Down-P0.05"=1,"NoSig"=0.7)
# 
# maxValue<-as.integer(max(abs(df$Log2FC))+1)
# maxY<-as.integer(max(-log10(df$ttestPvalue))+1)
# 
# mydf$NegLog10p<- -log10(mydf$ttestPvalue)
# mydf$ID<-rownames(mydf)
# 
# vcd<-mydf[,c("ID","Log2FC","NegLog10p")]
# 
# ###火山图绘制
# plotVC(vcd,outimage=paste0(pn,"/Volcano.png"), width = 6, fccutoff = FCcutoff, pcutoff=pcutoff)
# plotVC(vcd,outimage=paste0(pn,"/Volcano.pdf"), width = 6, fccutoff = FCcutoff, pcutoff=pcutoff)
# 
# 
# #火山图版本2
# #colours<-c(Up="#E41A1C",Down="#4DAF4A",NoSig="#999999")
# ##############　火山图　##########################
# ggplot(data=mydf, aes(x=Log2FC, y=-log10(ttestPvalue), colour=threshold, shape=threshold)) +
#   geom_point(alpha=0.7) +theme_bw()+
#   #scale_colour_manual("threshold: \nFC>=1.5 or FC<=0.6667\np-value<0.05",values=colours)+
#   scale_colour_manual(values=colours)+
#   scale_fill_manual(values=colours)+
#   scale_shape_manual(values=shapes)+
#   scale_size_manual(values=sizes)+
#   #theme(legend.position = "none",panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
#   theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.title = element_blank())+
#   xlim(c(-maxValue, maxValue)) + ylim(c(0, maxY)) +
#   xlab("log2 fold change") + ylab("-log10 p-value")+
#   #labs(title = maintitle)+
#   geom_hline(yintercept=-log10(0.05),colour="grey",size=.3,linetype=2)+
#   geom_hline(yintercept=-log10(0.01),colour="grey",size=.3,linetype=2)+
#   geom_vline(xintercept=log2(FCcutoff),colour="grey",size=.3,linetype=2)+
#   geom_vline(xintercept=-log2(FCcutoff),colour="grey",size=.3,linetype=2)+
#   annotate("text",-Inf,-log10(0.05),label="p=0.05",fontface="italic",hjust=-.2,vjust=1,size=2,colour="grey30")+
#   annotate("text",-Inf,-log10(0.01),label="p=0.01",fontface="italic",hjust=-.2,vjust=0,size=2,colour="grey30")+
#   annotate("text",log2(FCcutoff),Inf,label=paste("FC=",FCcutoff,sep=""),angle=90,hjust=1.2,vjust=1.5,size=2,colour="grey30")+
#   annotate("text",-log2(FCcutoff),Inf,label=paste("FC=",sprintf("%.3f", 1/FCcutoff),sep=""),angle=90,hjust=1.1,vjust=-.8,size=2,colour="grey30")
# ggsave(filename=paste(pn,"/VolcanoR_FC",FCcutoff,"_P0.05.png",sep="."),dpi=dpi,width=width,height = height, units = "cm")
# ggsave(filename=paste(pn,"/VolcanoR_FC",FCcutoff,"_P0.05.pdf",sep="."),dpi=1000,width=width,height = height, units = "cm")
# #quit()
# tmpdf<-as.data.frame(table(df$threshold))
# bardf<-tmpdf[tmpdf$Var1!="NoSig",]
# 
# mycol<-brewer.pal(10,"Paired")
# 
# library(ggsci)
# 
# p<-ggplot(bardf,aes(x=Var1,y=Freq,fill=Var1))+
#   geom_bar(width=.5,colour="black",stat="identity")+theme_bw()+
#   geom_text(aes(label=Freq),stat="identity",vjust=-.3,size=3,colour="black")+
#   theme(panel.border = element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x=element_blank(),
#         panel.grid.minor.x = element_blank(),legend.position = "none",axis.line = element_line(colour = "black",size=1))+
#   #scale_fill_manual(values=c(Up="#A50026",Down="#006837"))+
#   scale_y_continuous(expand=c(0,0),limits=c(0,max(bardf$Freq)*1.1))+
#   xlab("")+ylab("Count")+
#   scale_fill_npg()
# ggsave(filename=paste(pn,"/updown.bar.png",sep="."),plot=p,dpi=dpi,width=4,height=4,units="in")
# 
# mythreshold<-as.factor(mydf$threshold)
# 
# 
# write.table(summary(mythreshold),paste(pn,"/summary.txt",sep="."),sep="\t",col.names=FALSE)
# 
# statd<-as.data.frame(summary(mythreshold))
# statd$Term<-rownames(statd)
# statd<-statd[,c(2,1)]
# colnames(statd)<-c("Type","Count")
# drawd<-statd[statd$Type!='NoSig',]
# drawd$updown<-'DOWN'
# drawd$updown[drawd$Type=='Up-P0.01' | drawd$Type=='Up-P0.05'] <-'UP'
# drawd$ptype<-"p0.05"
# drawd$ptype[drawd$Type=='Up-P0.01' | drawd$Type=='Down-P0.01'] <-"p0.01"
# library(plyr)
# 
# drawd<-arrange(drawd,updown, ptype)
# drawd<-ddply(drawd, "updown", transform,label_y=cumsum(Count))
# drawd
# drawd$ptype<-factor(drawd$ptype,levels =c("p0.05","p0.01"))
# p<-ggplot(drawd, aes(x=updown, y=Count, fill=ptype))+geom_bar(stat="identity",width=0.5)+
#   geom_text(aes(y=label_y, label=Count), vjust=1.5, colour='white')+
#   theme_cowplot()+
#   scale_y_continuous(expand=c(0,0))+guides(fill=guide_legend(title=NULL))+
#   xlab("Significant differentially expressed protein")+scale_fill_npg()
# ggsave(paste0(pn,"/DEP.stackbar.png"), plot=p,width=6, height=6, units='in', dpi=300)
#   
# 
# #plotBarText(drawd[,c(1,2,2,1)],legendPosition = 'none', outimage=paste0(pn,"/DEP.bar.png"), width=4,height=4)
# 
# plotPie(statd, outimage=paste0(pn,"/ProteinType.pie.png"),width=6, height=4.5)

