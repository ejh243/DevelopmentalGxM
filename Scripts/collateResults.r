### sumarize results of gene x methylation correlations

setwd("/mnt/data1/Helen/GExpxDNAm/")

all.output<-list.files(path = "Output", pattern = "Cor_GExp_TM_Matrix")

cor.pearson<-NULL
for(file in all.output[-grep("Rank", all.output)]){
	tmp<-read.table(paste("Output/", file, sep = ""), header = TRUE)
	tmp<-t(tmp)
	cor.pearson<-rbind(tmp, cor.pearson)
}
save(cor.pearson, file = "Output/All_Pearson_Cor_TM.rdata")

cor.rank<-NULL
for(file in all.output[grep("Rank", all.output)]){
	tmp<-read.table(paste("Output/", file, sep = ""), header = TRUE)
	tmp<-t(tmp)
	cor.rank<-rbind(tmp, cor.rank)
}

save(cor.rank, file = "Output/All_Rank_Cor_TM.rdata")
