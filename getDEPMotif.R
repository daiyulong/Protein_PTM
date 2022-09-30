ARGS<-commandArgs(TRUE)

library(stringr)
library(openxlsx)
################### two sample test ######################
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


if(!dir.exists(paste0(pn,"/DEPMotif"))){
  dir.create(paste0(pn,"/DEPMotif"))
}


d<-read.delim(paste0(pn, "/MotifAnalysis/Motif_All.txt"), stringsAsFactors = FALSE,check.names = FALSE)

dUp <- read.delim(paste0(pn,"/",comp,"/DESelection/Up.txt"), stringsAsFactors = FALSE,check.names = FALSE)
dDown <- read.delim(paste0(pn,"/",comp,"/DESelection/Down.txt"), stringsAsFactors = FALSE,check.names = FALSE)
dDEP <- read.delim(paste0(pn,"/",comp,"/DESelection/DEP.txt"), stringsAsFactors = FALSE,check.names = FALSE)

rownames(d) <- d[[1]]
ddUp <- d[dUp[[1]],]
ddDown <- d[dDown[[1]],]
ddDEP <- d[dDEP[[1]],]

output <- paste0(pn,"/DEPMotif/",comp)
if (!dir.exists(output)) {
  dir.create(output)
}

write.table(ddUp[2],paste0(output,"/UpMotif_input.txt"),row.names = F,col.names = F,quote = F)
write.table(ddDown[2],paste0(output,"/DownMotif_input.txt"),row.names = F,col.names = F,quote = F)
write.table(ddDEP[2],paste0(output,"/DEPMotif_input.txt"),row.names = F,col.names = F,quote = F)



