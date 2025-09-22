library(ape)
library(Rlsd2)
library(BactDating)

seq_len <- 2232025 ##NCCP11945 NC_011035.1 sequence length
resB <- readRDS("big_timed_filtered.rds")
dates_tab <- read.table("metadata/dates.tab")

tr <- resB$inputtree
tr$edge.length <- tr$edge.length / seq_len #Rescale gubbins tree to subs/site

dates <- apply(dates_tab, 1, median)
names(dates) <- rownames(dates_tab)

dfile <- tempfile()
cat(length(dates),"\n",file = dfile)
for (i in 1:length(dates)){
  cat(names(dates)[i], dates[i], "\n", append = T, file = dfile)
}

tfile <- tempfile()
write.tree(tr, tfile)

resL <- lsd2(inputTree=tfile, inputDate=dfile, seqLen = seq_len, estimateRoot="a", confidenceInterval = 100)
print("LSD Result")
print(readLines(resL$outResultFiles[1])[32:43])
print("Bactdating Result")
v=resB$record[,Ntip(resB$tree)+1]
v=sort(v[(1+length(v)/2):length(v)])
vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
cat(sprintf('Root date=%.2f [%.2f;%.2f]\n',vals[1],vals[2],vals[3]))
v=resB$record[,"mu"]
v=v[(1+length(v)/2):length(v)]
v=sort(v)
vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])/2232025
cat(sprintf('%s=%.2e [%.2e;%.2e]\n',"Rate",vals[1],vals[2],vals[3]))
