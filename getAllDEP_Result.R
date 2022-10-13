# 基于ProbabilitySite.xlsx；ModifiedProtein.xlsx；DEP.xls；DEPProtein.xls共计4个文件
# 生成差异可靠修饰位点表 & 差异修饰蛋白表

ARGS<-commandArgs(TRUE)

library(openxlsx)

if(length(ARGS)<3){
  cat("Rscript getAllDEP_Result.R [comp] [type] [pn] \n")
  q()
}

comp <- as.character(ARGS[1])
type  <- as.character(ARGS[2])
pn <- as.character(ARGS[3])

if(FALSE){
  comp="case.vs.ctrl"
  type="s"   # s:Spectronaut; m:MSFragger; p:PD
  pn="E:/projectResearch/02PTM_WorkFlow/TEST"
  
}

allSite <- read.xlsx(paste0(pn,"/1.Info/ProbabilitySite.xlsx"),1)
allProtein <- read.xlsx(paste0(pn, "/1.Info/ModifiedProtein.xlsx"),1)

depSite <- read.delim(paste0(pn,"/",comp,"/DESelection/DEP.xls"),sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
depProtein <- read.delim(paste0(pn,"/",comp,"/DESelection/DEPProtein.xls"),sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)

if (type=="s") {

  selectColumn <- c("PTM.CollapseKey","ttestPvalue","FC","Log2FC","fdr","threshold")
  psite <- merge(allSite,depSite[,selectColumn],by.x = "PTM.CollapseKey",by.y = "PTM.CollapseKey")
  pprotein <- merge(depProtein[1],allProtein,by.x = "PG.ProteinAccessions",by.y ="PG.ProteinAccessions" )

}else if (type=="m") {
    
  selectColumn <- c("Protein.Position","ttestPvalue","FC","Log2FC","fdr","threshold")
  psite <- merge(allSite,depSite[,selectColumn],by.x = "Protein.Position",by.y = "Protein.Position")
  pprotein <- merge(depProtein[1],allProtein,by.x = "Protein",by.y ="Protein" )

}else if (type == "p") {
  cat("暂不支持PD搜库结果文件!\n")
  q()
}else if (type != "s"|type != "m"|type != "p"){
  cat("Input file error,please check it.\n")
  q()
}

# 生成差异可靠修饰位点表 & 差异修饰蛋白表
write.xlsx(psite,paste0(pn,"/1.Info/DEPProbabilitySite.xlsx"),overwrite = T)
write.xlsx(pprotein,paste0(pn,"/1.Info/DEPModifiedProtein.xlsx"),overwrite = T)
