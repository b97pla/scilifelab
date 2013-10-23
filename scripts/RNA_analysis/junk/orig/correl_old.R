args = commandArgs(TRUE)
snames = args[c(1:length(args))]

# Find no of rows
name <- snames[1]
fname <- paste("tophat_out_", name, "/cufflinks_out_", name, "/genes.fpkm_tracking", sep="")
temp <- read.table(fname, header=T)

e <-matrix(nrow=nrow(temp), ncol= (length(snames)))
colnames(e) <- rep("NO_NAME_ASSIGNED", length(snames))

idx <- 0 

for (name in snames) {
    idx <- idx+1
    fname <- paste("tophat_out_", name, "/cufflinks_out_", name, "/genes.fpkm_tracking", sep="")
    temp <- read.table(fname, header=T)
    colnames(e)[idx] <- name
    o <- order(temp$tracking_id)
    e[,idx] <- temp$FPKM[o]
    temp.o <- temp[o,]
    print(temp.o[1, ])
} 

# write a more user friendly version file 
d.0 <- data.frame(e)
d.1 <- cbind(Gene_ID=temp[o,5], d.0)
d.2 <- cbind(ENSEMBL_ID=temp[o,1],d.1)


########## cleanup cufflincs fpkm table
end=ncol(d.2)

#uniq genes
uniq=d.2[unique(d.2[,1]),]

#duplicate genes
dup=d.2[duplicated(d.2[,1]),]

#add rowsums
dup[,'suma']=rowSums(dup[,3:end])
#sort on rowsums
dup=dup[with(dup,order(suma)),]

#rm duplicated genses and ceep the genes w highest rowsums
max_dup=dup[!duplicated(dup[,'suma']),1:end]

#oncat uniq and dupl
all_uniq=rbind(uniq,max_dup)

#remoove genes w 0 coverage
keep = rowSums(all_uniq[3:end])>0
all_uniq = all_uniq[keep, ]

e=all_uniq[,3:end]
##########

write.table(d.2, file="fpkm_table.txt", quote=F, row.names=F, sep="\t")

# continue with numeric e
pdf("FPKM_heatmap.pdf")
if (end>20){heatmap(cor(e), symm=T,cexCol=0.4,cexRow=0.4)}else{heatmap(cor(e), symm=T)}
dev.off()

library(MASS)
library(calibrate)
p <- princomp(e)
l <- p$loadings[,1:2]
expl.vars <- (p$sdev^2)/sum(p$sdev^2)
pdf("FPKM_PCAplot.pdf")
plot(l, xlab=paste("PC 1, expl. var.", signif(expl.vars[1],2), sep=""), ylab=paste("PC 2; expl. var.", signif(expl.vars[2]), sep=""))
textxy(l[,1],l[,2],labs=rownames(l),cx=1)
dev.off()

#pdf("spearman_heatmap.pdf")
#heatmap(as.matrix(cor(e, method="spearman")))
#dev.off()

# HTSeq counts

name <- snames[1]
fname <- paste("tophat_out_", name, "/", name, ".counts", sep="")
nlines <- length(count.fields(fname))
temp <- read.table(fname, header=F, nrow = (nlines-5))
h <- matrix(nrow=nrow(temp), ncol=length(snames))
colnames(h) <- rep("NO_NAME_ASSIGNED", length(snames))
rownames(h) <- temp[,1]
idx <- 0

for (name in snames){
    idx <- idx + 1
    fname <- paste("tophat_out_", name, "/", name, ".counts", sep="")
    temp <- read.table(fname, header=F, nrow = (nlines-5))
    colnames(h)[idx]<-name
    h[,idx]<-temp[,2]

}

write.table(h, file="count_table.txt", quote=F, row.names=T, sep="\t")

pdf("count_heatmap.pdf")
heatmap(cor(h), symm=T)
dev.off()

#library(calibrate)
p <- princomp(h)
l <- p$loadings[,1:2]
expl.vars <- (p$sdev^2)/sum(p$sdev^2)
pdf("count_PCAplot.pdf")
plot(l, xlab=paste("PC 1 (expl. var.", signif(expl.vars[1]), ")", sep=""), ylab=paste("PC 2 (expl. var.", signif(expl.vars[2]), ")", sep=""))
textxy(l[,1],l[,2],labs=rownames(l),cx=1)
dev.off()
