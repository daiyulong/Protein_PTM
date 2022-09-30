ARGS<-commandArgs(TRUE)

library(openxlsx)

if(length(ARGS)<2){
  cat("Rscript getDEPMotif.R comp pn\n")
  q()
}

comp <- as.character(ARGS[1])
pn  <- as.character(ARGS[2])


if(FALSE){
  comp="h6g5l1l.vs.h6wtl"
  pn="E:/projectResearch/02PTM_WorkFlow/TEST/DEPmotif"

}


dUp <- read.delim(paste0(pn,"/",comp,"/DESelection/Up.txt"), stringsAsFactors = FALSE,check.names = FALSE)
dDown <- read.delim(paste0(pn,"/",comp,"/DESelection/Down.txt"), stringsAsFactors = FALSE,check.names = FALSE)
dDEP <- read.delim(paste0(pn,"/",comp,"/DESelection/DEP.txt"), stringsAsFactors = FALSE,check.names = FALSE)

allP <- read.xlsx(paste0(pn, "/1.Info/ModifiedProtein.xlsx"),1)
allProtein <- read.delim(paste0(pn,"/input.txt"),check.names = F,sep = "\t")

pUp <- merge(dUp[1],allP[,c("PG.ProteinAccessions","PTM.CollapseKey" )],by.x = "Accession",by.y = "PTM.CollapseKey")
pup2 <- merge(pUp[2],allProtein[,1:2 ],by.x = "PG.ProteinAccessions",by.y ="PG.ProteinAccessions" )

pDown <- merge(dDown[1],allP[,c("PG.ProteinAccessions","PTM.CollapseKey" )],by.x = "Accession",by.y = "PTM.CollapseKey")
pdown2 <- merge(pDown[2],allProtein[,1:2],by.x = "PG.ProteinAccessions",by.y ="PG.ProteinAccessions" )

pDEP <- merge(dDEP[1],allP[,c("PG.ProteinAccessions","PTM.CollapseKey" )],by.x = "Accession",by.y = "PTM.CollapseKey")
pdep2 <- merge(pDEP[2],allProtein[,1:2],by.x = "PG.ProteinAccessions",by.y ="PG.ProteinAccessions" )

write.table(pup2,paste0(pn,"/",comp,"/DESelection/UpProtein.xls"),sep = "\t",row.names = F,quote = F)
write.table(pdown2,paste0(pn,"/",comp,"/DESelection/DownProtein.xls"),sep = "\t",row.names = F,quote = F)
write.table(pdep2,paste0(pn,"/",comp,"/DESelection/DEPProtein.xls"),sep = "\t",row.names = F,quote = F)


