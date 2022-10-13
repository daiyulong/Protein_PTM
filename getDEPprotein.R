ARGS<-commandArgs(TRUE)

library(openxlsx)

if(length(ARGS)<3){
  cat("Rscript getDEPMotif.R [comp] [pn] [type]\n")
  q()
}

comp <- as.character(ARGS[1])
pn  <- as.character(ARGS[2])
type <- as.character(ARGS[3])

if(FALSE){
  comp="D100.vs.Control"
  pn="E:/projects_2022/BP20220597-1"
  type="m"   # s:Spectronaut; m:MSFragger; p:PD
}

if (type=="s") {
  
  for (dtype in c("Up","Down","DEP")) {
  
  dUp <- read.delim(paste0(pn,"/",comp,"/DESelection/",dtype,".txt"), stringsAsFactors = FALSE,check.names = FALSE)

  allP <- read.xlsx(paste0(pn, "/1.Info/ModifiedProtein.xlsx"),1)
  allProtein <- read.delim(paste0(pn,"/input.txt"),check.names = F,sep = "\t")
  
  pUp <- merge(dUp[1],allP[,c("PG.ProteinAccessions","PTM.CollapseKey" )],by.x = "Accession",by.y = "PTM.CollapseKey")
  pup2 <- merge(pUp[2],allProtein[,1:2 ],by.x = "PG.ProteinAccessions",by.y ="PG.ProteinAccessions" )
   
  write.table(pup2,paste0(pn,"/",comp,"/DESelection/",dtype,"Protein.xls"),sep = "\t",row.names = F,quote = F)
  }
}else if (type=="m") {
  
  for (dtype in c("Up","Down","DEP")) {
    
    dUp <- read.delim(paste0(pn,"/",comp,"/DESelection/",dtype,".txt"), stringsAsFactors = FALSE,check.names = FALSE)
    
    allP <- read.xlsx(paste0(pn, "/1.Info/ModifiedProtein.xlsx"),1)
    allProtein <- read.delim(paste0(pn,"/input.txt"),check.names = F,sep = "\t")
    
    pUp <- merge(dUp[1],allP[,c("Protein","Protein.Position" )],by.x = "Accession",by.y = "Protein.Position")
    pup2 <- merge(pUp[2],allProtein[,1:2],by.x = "Protein",by.y ="Protein" )
    
    write.table(pup2,paste0(pn,"/",comp,"/DESelection/",dtype,"Protein.xls"),sep = "\t",row.names = F,quote = F)
  }
  
}else if (type == "p") {
  cat("暂不支持PD搜库结果文件!\n")
  q()
}else if (type != "s"|type != "m"|type != "p"){
  cat("Input file error,please check it.\n")
  q()
}
