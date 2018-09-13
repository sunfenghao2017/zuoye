 Args <- commandArgs()
 inputmatrix <- Args[6]
 inputcondition <- Args[7]
 outputpath <- Args[8]
 lgFDl <- as.numeric(Args[9])
 lgFDu <- as.numeric(Args[10])
 padj <- as.numeric(Args[11])
 library(edgeR)
 library(gplots)
rt=read.csv(file =inputmatrix,header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,-1]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data=trunc(data)
data=na.omit(data)
summary(rowSums(data))
data <- data[rowSums(data)>10,]
group=read.csv(file=inputcondition,header=F,check.names=F)[,2]
design <- model.matrix(~group) 
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y) 
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = levels(group))
ordered_tags <- topTags(et, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts
write.table(diff,file=paste(outputpath,"/edgerOut.xls",sep=""),sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>lgFDu | diff$logFC<lgFDl)),]
write.table(diffSig, file=paste(outputpath,"/diffSig.xls",sep=""),sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>lgFDu)),]
write.table(diffUp, file=paste(outputpath,"/up.xls",sep=""),sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<lgFDl)),]
write.table(diffDown, file=paste(outputpath,"/down.xls",sep=""),sep="\t",quote=F)
normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file=paste(outputpath,"/normalizeExp.txt",sep=""),sep="\t",quote=F,col.names=F)
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file=paste(outputpath,"/diffmRNAExp.txt",sep=""),sep="\t",quote=F,col.names=F) 
heatmapData <- newData[rownames(diffSig),]
pdf(file=paste(outputpath,"/vol.pdf",sep="")) 
xMax=max(-log10(allDiff$FDR))+1
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>lgFDu,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<lgFDl,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()
hmExp=log10(heatmapData+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file=paste(outputpath,"/heatmap.pdf",sep=""),width=60,height=90)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()
