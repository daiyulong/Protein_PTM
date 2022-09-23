# get probability site table and probability modified probability protein table.

ARGV<-commandArgs(TRUE)

para<-as.character(ARGV[1])
input<-as.character(ARGV[2])
knninput<-as.character(ARGV[3])
pn<-as.character(ARGV[4])


if(FALSE){
  para <- "E:/projectResearch/02磷酸化蛋白组学流程/demo_parameter.xlsx"
  input <-"E:/projectResearch/02磷酸化蛋白组学流程/demo_sites.xlsx" 
  knninput <- "E:/projectResearch/02磷酸化蛋白组学流程/missValue_imputation.xlsx"
  pn <- "E:/projectResearch/02磷酸化蛋白组学流程/parameterResult"
}

suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(stringr))

if(stringr::str_ends(input,'xlsx')){
  dat<-read.xlsx(input,1)
}else{
  dat<-read.delim(input,check.names = FALSE)
}

if(stringr::str_ends(para,'xlsx')){
  dg<-read.xlsx(para,1)
}else{
  dg<-read.delim(para,check.names = FALSE)
}

if(stringr::str_ends(knninput,'xlsx')){
  datknn<-read.xlsx(knninput,1)
}else{
  datknn<-read.delim(knninput)
}


if (!dir.exists(paste0(pn,"/1.Info"))) {
  dir.create(paste0(pn,"/1.Info"))
}

sampNum <- nrow(dg)

dat <- dat[,1:(ncol(dat)-sampNum)]
dat[,(ncol(dat)-sampNum+1):ncol(dat)] <- lapply(dat[,(ncol(dat)-sampNum+1):ncol(dat)],FUN = function(y){as.numeric(y)})
dat$no <- 1:nrow(dat) 

# datknn2 <- datknn[,-1]
# colnames(datknn2) <- dg[[1]]
# datknn <- cbind(datknn[1],datknn2)


dat_merge <- merge(dat,datknn,by = "PTM.CollapseKey")
dat_merge <- dat_merge[order(dat_merge$no),]
dat_merge <- dat_merge[,-grep("no",colnames(dat_merge))]
#colnames(dat_merge)
dat_merge <- dat_merge[,c(2:6,1,7:ncol(dat_merge))]
#缺失值过滤并进行knn填充后的完整数据（filter NA + KNN imputation）
write.xlsx(dat_merge,paste0(pn,"/1.Info/AllSites_KnnImputation.xlsx"),overwrite = TRUE)
write.table(dat_merge[,c("PTM.CollapseKey",dg[[1]])],paste0(pn,"/1.Info/SampAnalysis_Input.txt"),
            sep = "\t",quote = F,row.names = F)

# 生成可信位点表和可信位点对应的蛋白表 

dat_prob <- dat_merge
proteID <-strsplit(dat_prob$PG.ProteinAccessions, ";")

dat_sub_pro <- dat_prob[,(ncol(dat_prob)-2*sampNum+1):(ncol(dat_prob)-sampNum)]
dat_sub_pro[dat_sub_pro=="Filtered"] <- 0


for (i in 1:nrow(dat_prob)) {
  dat_prob$pmax[i] <- max(as.numeric(dat_sub_pro[i,]))
  
  dat_prob$PAID[i] <- proteID[[i]][1]
}


# select site probability > 0.75, 通过修饰位点可信度最大值pmax，挑选可信修饰位点
win.seq <-data.frame(window= str_sub(dat_prob$PTM.FlankingRegion,start = 2,end = 14))
dat_prob <- add_column(dat_prob,window=win.seq[[1]],.after = "PTM.CollapseKey")

dat_prob75 <- filter(dat_prob,pmax > 0.75)
dup_proteAcce <- dat_prob75[!duplicated(dat_prob75$PAID),]
dat_prob75 <- dat_prob75[,!grepl("PAID|pmax",colnames(dat_prob75))]

# 附件3-可靠修饰位点表（filter NA + KNN imputation + probability > 0.75）
write.xlsx(dat_prob75,paste0(pn,"/1.Info/ProbabilitySite.xlsx"),overwrite = TRUE)


# 获得可信修饰位点对应的修饰蛋白表
dup_proteAcce$PG.ProteinAccessions <- dup_proteAcce$PAID
dup_proteAcce <- dup_proteAcce[,!grepl("PAID|pmax",colnames(dup_proteAcce))]
write.xlsx(dup_proteAcce[,c(1,4,2,3,5:ncol(dup_proteAcce))],paste0(pn,"/1.Info/ModifiedProtein.xlsx"),overwrite = TRUE)


###############################################################################################
# # 生成T.Test分析的matrix表
# para1 <- read.xlsx("parameter.xlsx",1) # 仅适用于2组间的比较
# para2 <- read.xlsx("parameter.xlsx",2)
# para2[[2]]
# para2[[3]]
# samp <- para1[grepl(paste0(para2[[2]],"|",para2[[3]]),para1$Sample),1]
# 
# ttest_dat <- dat_prob75[,c("PTM.CollapseKey",samp)]
# write.table(ttest_dat,paste0(para2[[2]],"-",para2[[3]],"-matrix.txt"),row.names = F,sep = "\t",quote = F)
# 
# 
# 
# 
# 
# 
# ################################################################################################
# ## P-value & FC 
# pfc_dat <- dat_prob75[,c(6,20:25)]
# groupF_dat<-read.delim("groupFile.txt", sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
# colnames(pfc_dat)[-1] <- groupF_dat[[1]]
# write.csv(pfc_dat,"pfc_data.txt",row.names = F)
# 
# #################################
# library(RColorBrewer)
# library(ggplot2)
# library(cowplot)
# library(hyplot)
# library(stringr)
# library(xlsx)
# library(openxlsx)
# input <- "pfc_data.txt"
# ginput <- "groupFile.txt"
# method <- "log2"
# FCcutoff <- 2
# pcutoff <- 0.05
# paired <- FALSE
# pn <- "./"
# 
# width=13.2
# height=10.5
# dpi=600
# units = "cm"
# 
# 
# getPostfix<-function(s){
#   v<-str_split(s,"\\.")
#   return(v[length(v)])
# }
# 
# loadData<-function(input){
#   if(getPostfix(input)=="csv"){
#     df<-read.csv(input, row.names=1, stringsAsFactors = FALSE, check.names =FALSE)
#   }else if(getPostfix(input)=="xlsx"){
#     df<-xlsx::read.xlsx(input, sheetIndex = 1)
#     rownames(df)<-df[[1]]
#     df<-df[,-1]
#   }else{
#     df<-read.delim(input,row.names=1,  sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
#   }
#   return(df)
# }
# 
# df<-loadData(input)
# gd <- loadData(ginput)
# gd$g <- rep(c(0,1),each=3)  ## modified  control:0 , Tumer:1 ,  FC  1/0
# 
# cl <- gd[[2]]
# m<-as.matrix(df)
# mbak<-as.matrix(df)
# 
# if(method == "log2"){
#   m<-log2(m)
#   transm<-log2(df)
#   colnames(transm)<-paste0("Log2", colnames(transm))
#   df<-cbind(df,transm)
# }else if(method =="log2plus1"){
#   m<-log2(m+1)
#   transm<-log2(df+1)
#   colnames(transm)<-paste0("Log2", colnames(transm))
#   df<-cbind(df,transm)
# }
# 
# cv<-function(x){
#   100*sd(x)/mean(x)
# }
# 
# se<-function(x){
#   sd(x)/sqrt(length(x))
# }
# 
# multi.ttest<-function(x){
#   errorflag<-try(t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE),silent=TRUE)
#   if('try-error' %in% class(errorflag)) NA
#   else t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE)$p.value
# }
# 
# multi.pairedttest<-function(x){
#   errorflag<-try(t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE, paired=TRUE),silent=TRUE)
#   if('try-error' %in% class(errorflag)) NA
#   else t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE, paired=TRUE)$p.value
# }
# 
# multi.wilcox.test<-function(x){
#   errorflag<-try(wilcox.test(x[which(cl==0)],x[which(cl==1)],exact=FALSE,correct=FALSE),silent=TRUE)
#   if('try-error' %in% class(errorflag)) NA
#   else wilcox.test(x[which(cl==0)],x[which(cl==1)],exact=FALSE,correct=FALSE)$p.value
# }
# 
# multi.var.test<-function(x){
#   var.test(x[which(cl==0)],x[which(cl==1)])$p.value
# }
# 
# multi.fc<-function(x){
#   mean(x[which(cl==1)]) / mean(x[which(cl==0)])
# }
# 
# if(paired==TRUE){
#   df$ttestPvalue<-apply(m,1,multi.pairedttest)
# }else{
#   df$ttestPvalue<-apply(m,1,multi.ttest)
# }
# 
# df$ControlMean<-apply(m[,cl==0], 1, mean)
# df$ControlSD  <-apply(m[,cl==0], 1, sd)
# df$ControlSE  <-apply(m[,cl==0], 1, se)
# df$ControlCV  <-apply(m[,cl==0], 1, cv)
# 
# df$TestMean   <-apply(m[,cl==1], 1, mean)
# df$TestSD     <-apply(m[,cl==1], 1, sd)
# df$TestSE     <-apply(m[,cl==1], 1, se)
# df$TestCV     <-apply(m[,cl==1], 1, cv)
# 
# 
# #apply(m,1,multi.wilcox.test)->df$wilcoxtestPvalue
# #apply(m,1,multi.var.test)->df$FtestPvalue
# 
# df$FC <- apply(mbak,1,multi.fc)
# df$Log2FC<-log2(df$FC)
# df$fdr<-p.adjust(df$ttestPvalue,method="fdr")
# 
# 
# mydf<-df
# 
# df$threshold<-"NoSig"
# df$threshold[which(df$Log2FC>=log2(FCcutoff) & df$ttestPvalue <=pcutoff)]<-"Up"
# df$threshold[which(df$Log2FC<=log2(1/FCcutoff) & df$ttestPvalue <=pcutoff)]<-"Down"
# 
# outdf <- cbind(rownames(df), df)
# colnames(outdf)[1] <- "Accession"
# 
# # 获得可靠的差异修饰位点表
# write.table(outdf,file=paste0(pn,"/ttest.txt"), sep='\t', row.names=FALSE, quote=FALSE)
# write.table(outdf[outdf$threshold=="Up" | outdf$threshold=="Down",],file=paste0(pn,"/DEP.txt"), sep='\t', row.names=FALSE, quote=FALSE)
# write.table(outdf[outdf$threshold=="Up",],file=paste0(pn,"/Up.txt"), sep='\t', row.names=FALSE, quote=FALSE)
# write.table(outdf[outdf$threshold=="Down",],file=paste0(pn,"/Down.txt"), sep='\t', row.names=FALSE, quote=FALSE)
# 
# # 获得差异修饰蛋白表
# DEP_site <- outdf[outdf$threshold=="Up" | outdf$threshold=="Down",]
# DEP_UP_site <- outdf[outdf$threshold=="Up",]
# DEP_Down_site <- outdf[outdf$threshold=="Down",]
# # 获得差异修饰motif window序列
# DEP_site_motif_input <- merge(dat_prob75,DEP_site[,c("Accession","FC","Log2FC","fdr","threshold")],by.x = "PTM.CollapseKey",by.y = "Accession")
# DEP_UP_site_motif_input <- merge(dat_prob75,DEP_UP_site[,c("Accession","FC","Log2FC","fdr","threshold")],by.x = "PTM.CollapseKey",by.y = "Accession")
# DEP_Down_site_motif_input <- merge(dat_prob75,DEP_Down_site[,c("Accession","FC","Log2FC","fdr","threshold")],by.x = "PTM.CollapseKey",by.y = "Accession")
# 
# DEP_site_motif_input <- DEP_site_motif_input[,c(2:6,1,7:ncol(DEP_site_motif_input))]
# 
# write.csv(DEP_site_motif_input,"DEP_ProbabilitySite.csv",row.names = F)
# # write.table(DEP_UP_site_motif_input,"DEP_UP_site_motif_input.txt",row.names = F,quote = F,sep = "\t")
# # write.table(DEP_Down_site_motif_input,"DEP_Down_site_motif_input.txt",row.names = F,quote = F,sep = "\t")
# 
# write.table(DEP_site_motif_input[,"window"],"DEP_site_motifID.txt",row.names = F,col.names = F,quote = F,sep = "\t")
# write.table(DEP_UP_site_motif_input[,"window"],"DEP_UP_site_motifID.txt",row.names = F,col.names = F,quote = F,sep = "\t")
# write.table(DEP_Down_site_motif_input[,"window"],"DEP_Down_site_motifID.txt",row.names = F,col.names = F,quote = F,sep = "\t")
# 
# DEP_protein <- merge(dup_proteAcce,DEP_site[,c("Accession","FC","Log2FC","fdr","threshold")],by.x = "PTM.CollapseKey",by.y = "Accession")
# DEP_protein <- DEP_protein[,c(2:6,1,7:ncol(DEP_protein))]
# DEP_proteinID <- data.frame(PTM.CollapseKey=DEP_protein[,1])
# write.csv(DEP_protein,"DEP_Modified_protein.csv",row.names = F)
# write.table(DEP_proteinID,"DEP_Modified_protein.txt",row.names = F,col.names = T,sep = "\t" ,quote = F)





