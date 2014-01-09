library(pheatmap)
args = commandArgs(TRUE)
snames = args[c(1:length(args))]

###-------------- FPKM table for genes

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

########## cleanup cufflincs fpkm table ##########
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

#concat uniq and dupl
all_uniq=rbind(uniq,max_dup)

#remoove genes w 0 coverage
keep = rowSums(all_uniq[3:end])>0
all_uniq_0cov_rem = all_uniq[keep, ]

########## end cleening table ##########

write.table(all_uniq, file="fpkm_table.txt", quote=F, row.names=F, sep="\t")

all_uniq = all_uniq[,3:end]
all_uniq_0cov_rem = all_uniq_0cov_rem[,3:end]
pdf("FPKM_heatmap.pdf")
if (end>20){pheatmap(cor(all_uniq), symm=T,fontsize = 4)}else{pheatmap(cor(all_uniq), symm=T)}
dev.off()

###-------------- FPKM table for isoforms

# Find no of rows
name <- snames[1]
fname <- paste("tophat_out_", name, "/cufflinks_out_", name, "/isoforms.fpkm_tracking", sep="")
temp <- read.table(fname, header=T)

e <-matrix(nrow=nrow(temp), ncol= (length(snames)))
colnames(e) <- rep("NO_NAME_ASSIGNED", length(snames))

idx <- 0 

for (name in snames) {
    idx <- idx+1
    fname <- paste("tophat_out_", name, "/cufflinks_out_", name, "/isoforms.fpkm_tracking", sep="")
    temp <- read.table(fname, header=T)
    colnames(e)[idx] <- name
    o <- order(temp$tracking_id)
    e[,idx] <- temp$FPKM[o]
    temp.o <- temp[o,]
    print(temp.o[1, ])
} 

# write a more user friendly version file 
d.0 <- data.frame(e)
d.1 <- cbind(Transcript_ID=temp[o,5], d.0)
d.2 <- cbind(ENSEMBL_ID=temp[o,1],d.1)

write.table(d.2, file="isoform_fpkm_table.txt", quote=F, row.names=F, sep="\t")


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


