# Get probability site table and probability modified probability protein table.

# UPDATE Information
# 2022-10-11: 增加MSFragger搜库结果分析

ARGV<-commandArgs(TRUE)

if(length(ARGS)<5){
  cat("Rscript getProbability_v2.R [para] [input] [knninput] [pn] [type]\n")
  q()
}

para<-as.character(ARGV[1])
input<-as.character(ARGV[2])
knninput<-as.character(ARGV[3])
pn<-as.character(ARGV[4])
type <- as.character(ARGV[5])

if(FALSE){
  para <- "E:/projects_2022/BP20220597-1/parameter.xlsx"
  input <-"E:/projects_2022/BP20220597-1/PTM_site.xlsx" 
  knninput <- "E:/projects_2022/BP20220597-1/1.Info/missValue_imputation.xlsx"
  pn <- "E:/projects_2022/BP20220597-1"
  type <- "m"  # s:Spectronaut; m:MSFragger; p:PD
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

if (type == "s") {
  
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
  write.table(dat_prob75[,c("PTM.CollapseKey",dg[[1]])],paste0(pn,"/1.Info/matrix.txt"),
              sep = "\t",quote = F,row.names = F)
  
  # 获得可信修饰位点对应的修饰蛋白表
  dup_proteAcce$PG.ProteinAccessions <- dup_proteAcce$PAID
  dup_proteAcce <- dup_proteAcce[,!grepl("PAID|pmax",colnames(dup_proteAcce))]
  write.xlsx(dup_proteAcce[,c(1,4,2,3,5:ncol(dup_proteAcce))],paste0(pn,"/1.Info/ModifiedProtein.xlsx"),overwrite = TRUE)

}else if (type == "m") {
  dat <- dat[,1:(ncol(dat)-sampNum)]
  #dat[,(ncol(dat)-sampNum+1):ncol(dat)] <- lapply(dat[,(ncol(dat)-sampNum+1):ncol(dat)],FUN = function(y){as.numeric(y)})
  dat$no <- 1:nrow(dat) 
  
  cat(paste0("rawProteinID: ",colnames(dat)[1],"\n"))
  cat(paste0("knnProteinID: ",colnames(datknn)[1],"\n"))
  
  dat_merge <- merge(dat,datknn,by = "Protein.Position")
  dat_merge <- dat_merge[order(dat_merge$no),]
  dat_merge <- dat_merge[,-grep("no",colnames(dat_merge))]

  #缺失值过滤并进行knn填充后的完整数据（filter NA + KNN imputation）
  write.xlsx(dat_merge,paste0(pn,"/1.Info/AllSites_KnnImputation.xlsx"),overwrite = TRUE)
  
  # Notice:MSFragger的结果全部都是可信的修饰位点，故dat_prob75和dat_merge内容相同
  dat_prob75 <- dat_merge
  # 附件3-可靠修饰位点表（filter NA + KNN imputation + probability > 0.75）
  write.xlsx(dat_prob75,paste0(pn,"/1.Info/ProbabilitySite.xlsx"),overwrite = TRUE)
  write.table(dat_prob75[,c("Protein.Position",dg[[1]])],paste0(pn,"/1.Info/matrix.txt"),
              sep = "\t",quote = F,row.names = F)
  
  
  # 获得可信修饰位点对应的修饰蛋白表
  dup_proteAcce <- dat_prob75
  dup_proteAcce <- dup_proteAcce[!duplicated(dup_proteAcce$Protein),]
  write.xlsx(dup_proteAcce[,c(3,1,2,4:ncol(dup_proteAcce))],paste0(pn,"/1.Info/ModifiedProtein.xlsx"),overwrite = TRUE)
  
}else if (type == "p") {
  cat("暂不支持PD搜库结果文件!\n")
  q()
}else if (type != "s"|type != "m"|type != "p"){
  cat("Input file error,please check it.\n")
  q()
}




