#!/usr/bin/env Rscript
library(peer)
by_genes<-read.table("gene.expression.qn_ind.txt",header=T) #splicing_qqnorm.PSI.txt
model=PEER()
PEER_setPhenoMean(model,as.matrix(t(by_genes)))
dim(PEER_getPhenoMean(model))    
PEER_setNk(model,30)
PEER_setNmax_iterations(model,1000)
PEER_getNk(model)
PEER_update(model)
factors=PEER_getX(model)
rownames(factors)<-colnames(by_genes)
colnames(factors)<-c(1:30)
Alpha = PEER_getAlpha(model)
factor_relevance <- 1.0 / Alpha
result <- data.frame(Factor = 1:length(factor_relevance), Relevance = factor_relevance) 
write.table(result, "factor_relevance.txt", sep="\t", row.names=FALSE, quote=FALSE)
pdf(paste0(30,".r_covs.pdf"),width=8,height=8)
plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance", main="")
dev.off()
factors
write.table(factors ,"peer.txt",sep="\t",row.names=T,quote =FALSE)
