# 基于ProbabilitySite.xlsx；ModifiedProtein.xlsx；DEP.xls；DEPProtein.xls共计4个文件
# 1.生成差异可靠修饰位点表 & 差异修饰蛋白表
# 2.获得位点合并的原始蛋白表（基于spectronaut搜库）

ARGS<-commandArgs(TRUE)

library(openxlsx)
library(stringr)
if(length(ARGS)<3){
  cat("Rscript getAllDEP_Result.R [pf] [type] [pn] \n")
  q()
}
pf <- as.character(ARGS[1])
type  <- as.character(ARGS[2])
pn <- as.character(ARGS[3])

if(FALSE){
  pf="PTM_Protein.xlsx"
  type="s"   # s:Spectronaut; m:MSFragger; p:PD
  pn="E:/projectResearch/02PTM_WorkFlow/TEST"
  
}


# 获得比较组信息
L <- c()
d <- read.delim(paste0(pn,"/","comparison.txt"),sep = "\t",header = T,stringsAsFactors = FALSE,check.names = FALSE)

for (i in 1:NROW(d)) {
  g1=as.character(d[i,2])
  g2=as.character(d[i,3])
  stringr::str_trim(g1)
  stringr::str_trim(g2)
  
  aa <- paste0(g1,".vs.",g2)
  L <- append(L,aa)
  
}


# sheet_list <- list()
# sheet_list2 <- list()
# for (comp in L) {
#   allSite <- read.xlsx(paste0(pn,"/1.Info/ProbabilitySite.xlsx"),1)
#   allProtein <- read.xlsx(paste0(pn, "/1.Info/ModifiedProtein.xlsx"),1)
#   
#   datname<- read.delim(paste0(pn,"/",comp,"/DESelection/DEP.xls"),sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
#   depProtein <- read.delim(paste0(pn,"/",comp,"/DESelection/DEPProtein.xls"),sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
#   
  
if (type=="s") {
  
  # 获得位点合并的蛋白表
  df <- read.xlsx(pf,1)
  df1 <- df[df$PTM.SiteAA %in% c("S","T","Y"),]
  dup_id <- unique(df1$PG.ProteinAccessions)
  
  df_all <- data.frame()
  for (i in 1:length(dup_id)) {
    df_block <- subset(df1,PG.ProteinAccessions==dup_id[i])
    df_block$PTM.sites <- str_c(unique(df_block$PTM.CollapseKey),collapse = ";")
    df_all <- rbind(df_all,df_block)
  }
  df_all <- df_all[!duplicated(df_all$PG.ProteinAccessions),]
  write.xlsx(df_all,paste0(pn, "/1.Info/Protein_PTM.xlsx"),overwrite = T)
  
  
  ###差异可靠修饰位点表 & 差异修饰蛋白表
  sheet_list <- list()
  sheet_list2 <- list()
  for (comp in L) {
    allSite <- read.xlsx(paste0(pn,"/1.Info/ProbabilitySite.xlsx"),1)
    allProtein <- read.xlsx(paste0(pn, "/1.Info/ModifiedProtein.xlsx"),1)
    
    datname<- read.delim(paste0(pn,"/",comp,"/DESelection/DEP.xls"),sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
    depProtein <- read.delim(paste0(pn,"/",comp,"/DESelection/DEPProtein.xls"),sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
    
    selectColumn <- c("PTM.CollapseKey","ttestPvalue","FC","Log2FC","fdr","threshold")
    sheet_list[[comp]] <- merge(allSite,datname[,selectColumn],by.x = "PTM.CollapseKey",by.y = "PTM.CollapseKey")
    sheet_list2[[comp]] <- merge(depProtein[1],allProtein,by.x = "PG.ProteinAccessions",by.y ="PG.ProteinAccessions")
  }
  
}else if (type=="m") {
  ###差异可靠修饰位点表 & 差异修饰蛋白表
  sheet_list <- list()
  sheet_list2 <- list()
  for (comp in L) {
    allSite <- read.xlsx(paste0(pn,"/1.Info/ProbabilitySite.xlsx"),1)
    allProtein <- read.xlsx(paste0(pn, "/1.Info/ModifiedProtein.xlsx"),1)
    
    datname<- read.delim(paste0(pn,"/",comp,"/DESelection/DEP.xls"),sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
    depProtein <- read.delim(paste0(pn,"/",comp,"/DESelection/DEPProtein.xls"),sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
    
    selectColumn <- c("Protein.Position","ttestPvalue","FC","Log2FC","fdr","threshold")
    sheet_list[[comp]] <- merge(allSite,datname[,selectColumn],by.x = "Protein.Position",by.y = "Protein.Position")
    sheet_list2[[comp]] <- merge(depProtein[1],allProtein,by.x = "Protein",by.y ="Protein" )
  }
}else if (type == "p") {
  cat("暂不支持PD搜库结果文件!\n")
  q()
}else if (type != "s"|type != "m"|type != "p"){
  cat("Input file error,please check it.\n")
  q()
}

# 生成差异可靠修饰位点表 & 差异修饰蛋白表
write.xlsx(sheet_list,paste0(pn,"/1.Info/DEPProbabilitySite.xlsx"),overwrite = T)
write.xlsx(sheet_list2,paste0(pn,"/1.Info/DEPModifiedProtein.xlsx"),overwrite = T)
#}


