
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





