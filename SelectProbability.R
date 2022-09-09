library(openxlsx)
library(dplyr)
setwd("E:/projects_2022/20220818/infiles")

input <-"All-rawSite.xlsx" 
knninput <- "Imputer-KNN.xlsx"
sampNum <- 18
groupFile <- "groupFile.txt"


if(stringr::str_ends(input,'xlsx')){
  dat<-read.xlsx(input,1)
}else{
  dat<-read.delim(input,check.names = FALSE)
}

if(stringr::str_ends(groupFile,'xlsx')){
  dg<-read.xlsx(groupFile,1)
}else{
  dg<-read.delim(groupFile,check.names = FALSE)
}

if(stringr::str_ends(knninput,'xlsx')){
  datknn<-read.xlsx(knninput,1)
}else{
  datknn<-read.delim(knninput)
}


dat <- dat[,1:(ncol(dat)-sampNum)]
dat$no <- 1:nrow(dat) 

datknn2 <- datknn[,-1]
colnames(datknn2)
colnames(datknn2) <- dg[[1]]
datknn <- cbind(datknn[1],datknn2)
# 缺失值过滤并进行knn填充后的完整数据
dat_merge <- merge(dat,datknn,by = "PTM.CollapseKey")
dat_merge <- dat_merge[order(dat_merge$no),]
dat_merge <- dat_merge[,-grep("no",colnames(dat_merge))]

#dat_merge <- dat_merge[,c(2:6,1,7:ncol(dat_merge))]
#write.xlsx(dat_merge,"data_Impute.xlsx")

# select 

dat_prob <- dat_merge

#proteID <-strsplit(dat_prob$PG.ProteinAccessions, ";")
# 获得修饰位点可信度最大值和唯一蛋白group值
dat_sub_pro <- dat_prob[,(ncol(dat_prob)-2*sampNum+1):(ncol(dat_prob)-sampNum)]
dat_sub_pro[dat_sub_pro=="Filtered"] <- 0


for (i in 1:nrow(dat_prob)) {
  dat_prob$pmax[i] <- max(as.numeric(dat_sub_pro[i,]))
#  dat_prob$pmax[i] <- max(as.numeric(dat_prob[i,(ncol(dat_prob)-2*sampNum):(ncol(dat_prob)-sampNum)]))
#  dat_prob$PAID[i] <- proteID[[i]][1]
}


# select site probability > 0.75, 通过修饰位点可信度最大值pmax，挑选可信修饰位点
dat_prob75 <- filter(dat_prob,pmax > 0.75)
#dup_proteAcce <- dat_prob75[!duplicated(dat_prob75$PAID),]

# 删除掉PAID和pmax列
dat_prob75 <- dat_prob75[,!grepl("PAID|pmax",colnames(dat_prob75))]

win.seq <-data.frame(window= str_sub(dat_prob75$PTM.FlankingRegion,start = 2,end = 14))
dat_prob75 <- add_column(dat_prob75,window=win.seq[[1]],.after = "PTM.CollapseKey")

write.csv(dat_prob75,"ProbabilitySite.csv",row.names = F) # 附件3-可靠修饰位点表

# group.WT <- dg[grep("^WT",dg$Sample),]
# group.C1 <- dg[grep("^C1",dg$Sample),]
# group.M8 <- dg[grep("^M8",dg$Sample),]
# 
# write.table(group.WT,"group_wt.txt",quote = F,row.names = F,sep = "\t")
# write.table(group.C1,"group_c1.txt",quote = F,row.names = F,sep = "\t")
# write.table(group.M8,"group_m8.txt",quote = F,row.names = F,sep = "\t")
# 
# write.table(dat_prob75[,c(1,43:48)],"matrix_wt.txt",quote = F,row.names = F,sep = "\t")
# write.table(dat_prob75[,c(1,49:54)],"matrix_c1.txt",quote = F,row.names = F,sep = "\t")
# write.table(dat_prob75[,c(1,55:60)],"matrix_m8.txt",quote = F,row.names = F,sep = "\t")
# dfwt <- data.frame()
# for (i in length(group.WT[[1]])) {
#   strr <- group.WT[[1]][i]
#   matrix.WT <- dat_prob75[,grep(paste("^",strr,sep = ""),colnames(dat_prob75))]
#   dfwt <- cbind(dfwt,matrix.WT)
# }



# Deduplication of protein groups based on trusted sites.
# 获得可信修饰位点对应的修饰蛋白表
dup_proteAcce$PG.ProteinAccessions <- dup_proteAcce$PAID
dup_proteAcce <- dup_proteAcce[,-c(ncol(dup_proteAcce)-1,ncol(dup_proteAcce))]

write.csv(dup_proteAcce,"Modified_protein.csv",row.names = F)
 




## P-value FC 
pfc_dat <- dat_prob75[,c(6,20:25)]
groupF_dat<-read.delim("groupFile.txt", sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(pfc_dat)[-1] <- groupF_dat[[1]]
write.csv(pfc_dat,"pfc_data.txt",row.names = F)

#################################
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(hyplot)
library(stringr)
library(xlsx)
library(openxlsx)
input <- "pfc_data.txt"
ginput <- "groupFile.txt"
method <- "log2"
FCcutoff <- 2
pcutoff <- 0.05
paired <- FALSE
pn <- "./"

width=13.2
height=10.5
dpi=600
units = "cm"


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
gd <- loadData(ginput)
gd$g <- rep(c(0,1),each=3)  ## modified  control:0 , Tumer:1 ,  FC  1/0

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

# 获得可靠的差异修饰位点表
write.table(outdf,file=paste0(pn,"/ttest.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(outdf[outdf$threshold=="Up" | outdf$threshold=="Down",],file=paste0(pn,"/DEP.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(outdf[outdf$threshold=="Up",],file=paste0(pn,"/Up.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(outdf[outdf$threshold=="Down",],file=paste0(pn,"/Down.txt"), sep='\t', row.names=FALSE, quote=FALSE)

# 获得差异修饰蛋白表
DEP_site <- outdf[outdf$threshold=="Up" | outdf$threshold=="Down",]
DEP_UP_site <- outdf[outdf$threshold=="Up",]
DEP_Down_site <- outdf[outdf$threshold=="Down",]
# 获得差异修饰motif window序列
DEP_site_motif_input <- merge(dat_prob75,DEP_site[,c("Accession","FC","Log2FC","fdr","threshold")],by.x = "PTM.CollapseKey",by.y = "Accession")
DEP_UP_site_motif_input <- merge(dat_prob75,DEP_UP_site[,c("Accession","FC","Log2FC","fdr","threshold")],by.x = "PTM.CollapseKey",by.y = "Accession")
DEP_Down_site_motif_input <- merge(dat_prob75,DEP_Down_site[,c("Accession","FC","Log2FC","fdr","threshold")],by.x = "PTM.CollapseKey",by.y = "Accession")

DEP_site_motif_input <- DEP_site_motif_input[,c(2:6,1,7:ncol(DEP_site_motif_input))]

write.csv(DEP_site_motif_input,"DEP_ProbabilitySite.csv",row.names = F)
# write.table(DEP_UP_site_motif_input,"DEP_UP_site_motif_input.txt",row.names = F,quote = F,sep = "\t")
# write.table(DEP_Down_site_motif_input,"DEP_Down_site_motif_input.txt",row.names = F,quote = F,sep = "\t")

write.table(DEP_site_motif_input[,"window"],"DEP_site_motifID.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(DEP_UP_site_motif_input[,"window"],"DEP_UP_site_motifID.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(DEP_Down_site_motif_input[,"window"],"DEP_Down_site_motifID.txt",row.names = F,col.names = F,quote = F,sep = "\t")

DEP_protein <- merge(dup_proteAcce,DEP_site[,c("Accession","FC","Log2FC","fdr","threshold")],by.x = "PTM.CollapseKey",by.y = "Accession")
DEP_protein <- DEP_protein[,c(2:6,1,7:ncol(DEP_protein))]
DEP_proteinID <- data.frame(PTM.CollapseKey=DEP_protein[,1])
write.csv(DEP_protein,"DEP_Modified_protein.csv",row.names = F)
write.table(DEP_proteinID,"DEP_Modified_protein.txt",row.names = F,col.names = T,sep = "\t" ,quote = F)





