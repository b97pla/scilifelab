
args <- commandArgs(TRUE)

print(args)

mycolors<-rainbow(length(args))

for(i in 1:length(args)){
	data<-read.table(args[i],header=TRUE)
	xp=data$TOTAL_READS[seq(0,length(data$TOTAL_READS),100)]
	yp=data$EXPECTED_DISTINCT[seq(0,length(data$EXPECTED_DISTINCT),100)]
	plot(xp,yp,col=mycolors[i], type="o", cex=0.3, axes=F,xlab="",ylab="")
	par(new=TRUE)
}
	m=max(xp)
	mv=seq(0,m,ceiling(m/100))
	print(mv)
	axis(1,at=xp,labels=xp)
	axis(2,at=ceiling(mv), labels=ceiling(mv), cex.axis=0.8, las=2)	
	title(xlab="TOTAL READS")
	title(ylab="EXPECTED DISTINCT")

legend("bottomright", args, text.col=mycolors)
