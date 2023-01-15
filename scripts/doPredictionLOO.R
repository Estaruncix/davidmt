#!/usr/bin/env Rscript
###############

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
#meta.train.file <- as.character(parameters[2])
#meta.vali.file <- as.character(parameters[3])
#boruta.res.file <- as.character(parameters[4])
myfactors <- as.character(parameters[2])
mylevels.ch <- as.character(parameters[3])
study <- as.character(parameters[4])
mynorm <- as.logical(as.character(parameters[5]))
# Este threshold se usará para hacer cálculos sobre las prediction.extra.tables (filtraremos individuos con probabilidad
# por debajo del threshold antes de calcular el pct.of.CA.for) y se almacenarán los resultados de todos los modelos en una única tabla.
# Si se quiere deshabilitar poner a 0.5.
threshold <- as.numeric(as.character(parameters[6]))

source("myfunctions.downstream.R")

library(randomForest)
library(pracma)

res.all <- c()
mypred.all <- c()

assess.RF <- function(pred.test, prob.cut.positive=0.5, positive.level="YES", ALL=NULL, ALL.classes){
	true.test <- ALL.classes
	if(!is.null(ALL)) names(true.test) <- colnames(ALL) else names(true.test) <- names(ALL.classes)
	true.test <- true.test[names(pred.test)]

	pred.test.positive <- names(pred.test)[which(pred.test>prob.cut.positive)]
	pred.test.negative <- names(pred.test)[which(pred.test<=prob.cut.positive)]
	true.test.positive <- names(true.test)[which(true.test==positive.level)]
	true.test.negative <- names(true.test)[which(true.test!=positive.level)]
	VP <- length(intersect(pred.test.positive, true.test.positive))
	FP <- length(intersect(pred.test.positive, true.test.negative))
	VN <- length(intersect(pred.test.negative, true.test.negative))
	FN <- length(intersect(pred.test.negative, true.test.positive))
	res <- c(VP/(VP+FN), FP/(FP+VN), VP, FP, VN, FN)
	names(res) <- c("VPR", "FPR", "VP", "FP", "VN", "FN")
	return(res)
}

boruta.files <- system("ls results.boruta.LOO.*", intern=TRUE)

for(i in 1:length(boruta.files)){

meta.train.file <- paste("meta.train", i, "tsv", sep=".")
meta.vali.file <- paste("meta.vali", i, "tsv", sep=".")
boruta.res.file <- paste("results.boruta.LOO", i, "tsv", sep=".")

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(dat) <- dat[, 1]
rownames(dat) <- gsub(";", ".", rownames(dat))
rownames(dat) <- gsub("_", ".", rownames(dat))
rownames(dat) <- gsub("-", ".", rownames(dat))
rownames(dat) <- gsub("[(]", ".", rownames(dat))
rownames(dat) <- gsub("[)]", ".", rownames(dat))
rownames(dat) <- gsub("[[]", ".", rownames(dat))
rownames(dat) <- gsub("[]]", ".", rownames(dat))
rownames(dat) <- gsub(" ", ".", rownames(dat))
rownames(dat) <- gsub(":", ".", rownames(dat))
rownames(dat) <- gsub("/", ".", rownames(dat))

meta <- read.table(file=meta.train.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

myfactors.l <- strsplit(myfactors, split=",")
myfactor1 <- myfactors.l[[1]][1]
if(length(myfactors.l[[1]])==2) myfactor2 <- myfactors.l[[1]][2] else myfactor2 <- ""

if(myfactor2!="") GROUPS <- paste(meta[, myfactor1], meta[, myfactor2], sep=".") else GROUPS <- as.character(meta[, myfactor1])
names(GROUPS) <- mysamples
mylevels <- unlist(strsplit(mylevels.ch, split=","))
mylevels <- mylevels[which(is.element(mylevels, unique(GROUPS)))]
ind <- which(is.element(GROUPS, mylevels))

dat <- dat[, ind]
GROUPS <- GROUPS[ind]
yes.group <- mylevels[1]
no.group <- mylevels[2]
meta <- meta[ind, ]

dat.new <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(dat.new) <- dat.new[, 1]
rownames(dat.new) <- gsub(";", ".", rownames(dat.new))
rownames(dat.new) <- gsub("_", ".", rownames(dat.new))
rownames(dat.new) <- gsub("-", ".", rownames(dat.new))
rownames(dat.new) <- gsub("[(]", ".", rownames(dat.new))
rownames(dat.new) <- gsub("[)]", ".", rownames(dat.new))
rownames(dat.new) <- gsub("[[]", ".", rownames(dat.new))
rownames(dat.new) <- gsub("[]]", ".", rownames(dat.new))
rownames(dat.new) <- gsub(" ", ".", rownames(dat.new))
rownames(dat.new) <- gsub(":", ".", rownames(dat.new))
rownames(dat.new) <- gsub("/", ".", rownames(dat.new))

rownames.dat.new <- rownames(dat.new)
meta.new <- read.table(file=meta.vali.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta.new) <- meta.new[, 1]
mysamples.new <- intersect(meta.new[, 1], colnames(dat.new))
dat.new <- as.matrix(dat.new[, mysamples.new])
rownames(dat.new) <- rownames.dat.new
meta.new <- meta.new[mysamples.new, ]

if(mynorm==TRUE){
	dat.all <- prop.table(as.matrix(dat), margin=2)*100
	dat.all.new <- prop.table(as.matrix(dat.new), margin=2)*100
} else {
	dat.all <- dat; dat.all.new <- dat.new
}

mytb <- read.table(file=boruta.res.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)
mytb$feature <- gsub(";", ".", mytb$feature)
mytb$feature <- gsub("_", ".", mytb$feature)
mytb$feature <- gsub("-", ".", mytb$feature)
mytb$feature <- gsub("[(]", ".", mytb$feature)
mytb$feature <- gsub("[)]", ".", mytb$feature)
mytb$feature <- gsub("[[]", ".", mytb$feature)
mytb$feature <- gsub("[]]", ".", mytb$feature)
mytb$feature <- gsub(" ", ".", mytb$feature)
mytb$feature <- gsub(":", ".", mytb$feature)
mytb$feature <- gsub("/", ".", mytb$feature)

mytb <- mytb[which(mytb$decision=="Confirmed"), ]

numvars <- nrow(mytb)

myvars <- intersect(mytb$feature[1:numvars], rownames(dat))

dat <- dat.all[myvars, ]
dat.new <- dat.all.new[myvars, ]

# Entrenamos el modelo con todas las muestras viejas y predecimos sobre las nuevas

my.mtry <- round(sqrt(nrow(dat)-1))

training <- as.data.frame(t(dat))
training$y <- as.factor(GROUPS)
testing <- as.data.frame(t(dat.new))

set.seed(777)

RF <- randomForest(y ~ ., data=training, mtry=my.mtry, ntree=1000, keep.forest=TRUE)

mypred <- predict(RF, type="prob", newdata=testing)

myclass <- apply(mypred, MARGIN=1, FUN=function(x){
		ind <- which(x==max(x))
		if(length(ind)>1) return(NA) else return(colnames(mypred)[which(x==max(x))])
	})

res <- as.data.frame(mypred)
res$class.predicted <- myclass
rownames(res) <- res$sample <- meta.new[, 1] 
res <- res[, c("sample", setdiff(colnames(res), "sample"))]

if(myfactor2!="") GROUPS.new <- paste(meta.new[, myfactor1], meta.new[, myfactor2], sep=".") else GROUPS.new <- as.character(meta.new[, myfactor1])
names(GROUPS.new) <- rownames(meta.new)

res$true.class <- GROUPS.new[rownames(res)]

res.all <- rbind(res.all, res)
mypred.all <- rbind(mypred.all, mypred)

print(i)
}

rownames(mypred.all) <- rownames(res.all)

# filtramos por threshold y calculamos pct.of.CA.for (los individuos filtrados no entran en la cuenta del total de grupo)
ind <- which(res.all[, 2]>=threshold|res.all[, 3]>=threshold)

if(length(ind)>0){
res.filt <- res.all[ind, ]; mypred.filt <- mypred.all[ind, ]
n.yes <- length(which(res.filt$class.predicted==yes.group&res.filt$true.class==yes.group))
n.no <- length(which(res.filt$class.predicted!=yes.group&res.filt$true.class!=yes.group))
tot.yes <- length(which(res.filt$true.class==yes.group))
tot.no <- length(which(res.filt$true.class!=yes.group))
pct.of.CA.for.YES <- 100*n.yes/tot.yes
pct.of.CA.for.NO <- 100*n.no/tot.no
abs.of.CA.for.YES <- n.yes
abs.of.CA.for.NO <- n.no
tot.yes.v <- tot.yes
tot.no.v <- tot.no
write.table(res.filt, file=paste(study, "probs.for.threshold", threshold, "tsv", sep="."),
		sep="\t", quote=FALSE, row.names=FALSE)
} else {
pct.of.CA.for.YES <- NA
pct.of.CA.for.NO <- NA
abs.of.CA.for.YES <- NA
abs.of.CA.for.NO <- NA
tot.yes.v <- NA
tot.no.v <- NA
}

file.pdf <- paste("ROC.RF.for.threshold", threshold, study, "pdf", sep=".")

pdf(file=file.pdf)

res <- res.filt; mypred <- mypred.filt

# pintamos curva roc de la predicción
ALL.classes <- res$true.class
ALL.classes <- as.factor(ALL.classes)
names(ALL.classes) <- rownames(res)
# Cortes de discriminación para el cálculo de la curva ROC
cuts <- seq(from=0, to=1, length.out=50)
cuts <- c(cuts, 0.5); cuts <- cuts[order(cuts, decreasing=FALSE)]
TPR.v <- FPR.v <- TP.v <- FP.v <- TN.v <- FN.v <- rep(NA, length(cuts))
names(TPR.v) <- names(FPR.v) <- names(TP.v) <- names(FP.v) <- names(TN.v) <- names(FN.v) <- as.character(cuts)
# evaluamos los RF para cada split y para cada cut
for(mycut in cuts){
	res <- assess.RF(pred.test=mypred[names(ALL.classes), yes.group], prob.cut.positive=mycut,
			positive.level=yes.group, ALL.classes=ALL.classes)
	if(!is.na(res["VPR"])&!is.na(res["FPR"])&res["VPR"]!=Inf&res["FPR"]!=Inf){
		TPR.v[as.character(mycut)] <- res["VPR"]
		FPR.v[as.character(mycut)] <- res["FPR"]
		TP.v[as.character(mycut)] <- res["VP"]
		FP.v[as.character(mycut)] <- res["FP"]
		TN.v[as.character(mycut)] <- res["VN"]
		FN.v[as.character(mycut)] <- res["FN"]
	}
}

cuts <- format(cuts, digits=1)

ind <- which(!is.na(TPR.v))
TPR.v <- TPR.v[ind]; FPR.v <- FPR.v[ind]; TP.v <- TP.v[ind]; FP.v <- FP.v[ind]; TN.v <- TN.v[ind]; FN.v <- FN.v[ind]; cuts <- cuts[ind]

AUC <- abs(trapz(FPR.v, TPR.v))

if(is.nan(AUC)){
	warning("El AUC sale NaN")
	return()
} else {
	AUC <- format(AUC, digits=2)
}

mymain <- paste("AUC:", AUC, "numvars:", numvars, sep=" ")

plot(x=FPR.v, y=TPR.v, type="n", xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), main=mymain)
abline(a=0, b=1, col="red")
ind <- c(which(cuts==0.53), which(cuts==0.63), which(cuts==0.73))
if(length(ind)>0) text(x=FPR.v[ind], y=TPR.v[ind], labels=cuts[ind], col="blue", cex=0.6)
lines(x=FPR.v, y=TPR.v, col="blue")

dev.off()

points.df <- data.frame(FPR=FPR.v, TPR=TPR.v, TP=TP.v, FP=FP.v, TN=TN.v, FN=FN.v, FNR=1-TPR.v, discrimination.threshold=cuts)
write.table(points.df, file=paste(study, "prediction.points.ROC.for.threshold", threshold, "tsv", sep="."),
		sep="\t", quote=FALSE, row.names=FALSE)

pct.true.pos.for.thr <- data.frame(numvars=numvars, abs.of.CA.for.YES=abs.of.CA.for.YES, abs.of.CA.for.NO=abs.of.CA.for.NO,
						tot.yes.v=tot.yes.v, tot.no.v=tot.no.v,
						pct.of.CA.for.YES=pct.of.CA.for.YES, pct.of.CA.for.NO=pct.of.CA.for.NO,
						tot.pct.of.CA=100*(abs.of.CA.for.YES+abs.of.CA.for.NO)/(tot.yes.v+tot.no.v),
						AUC=AUC)

colnames(pct.true.pos.for.thr) <- c("numvars", paste("num.of.CA.for", yes.group, sep="."), paste("num.of.CA.for", no.group, sep="."),
								paste("total", yes.group, sep="."), paste("total", no.group, sep="."),
								paste("pct.of.CA.for", yes.group, sep="."), paste("pct.of.CA.for", no.group, sep="."),
								"total.pct.of.CA", "AUC")

write.table(pct.true.pos.for.thr, file=paste("pct.true.pos.for.threshold", threshold, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)







