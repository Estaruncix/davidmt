#!/usr/bin/env Rscript
###############
# Genera particiones de un metadata para hacer leave one out
# ./doFoldsLOO.R $META "week" "pre,2w"

parameters <- commandArgs(trailingOnly=TRUE)

meta.file <- as.character(parameters[1])
myfactors <- as.character(parameters[2])
mylevels <- as.character(parameters[3])

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]

myfactors.l <- strsplit(myfactors, split=",")
myfactor1 <- myfactors.l[[1]][1]
if(length(myfactors.l[[1]])==2) myfactor2 <- myfactors.l[[1]][2] else myfactor2 <- ""

if(myfactor2!="") GROUPS <- paste(meta[, myfactor1], meta[, myfactor2], sep=".") else GROUPS <- as.character(meta[, myfactor1])
names(GROUPS) <- rownames(meta)

mylevels <- unlist(strsplit(mylevels, split=","))

mylevels <- mylevels[which(is.element(mylevels, unique(GROUPS)))]

ind <- which(is.element(GROUPS, mylevels))

meta <- meta[ind, ]

for(i in 1:nrow(meta)){
	ind.vali <- i
	meta.vali <- meta[ind.vali, ]
	meta.train <- meta[setdiff(1:nrow(meta), ind.vali), ]
	write.table(meta.vali, file=paste("meta.vali", i, "tsv", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
	write.table(meta.train, file=paste("meta.train", i, "tsv", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
}

