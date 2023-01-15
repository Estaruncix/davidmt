#!/usr/bin/env Rscript
###############
# Necesitamos generar varias copias idénticas (tantas como los folds especificados a doFoldsLOO) de la tabla de resultados de boruta.
# Estas copias serán las que entren como input a doPredictionLOO.R que se encargará de evaluar ese modelo concreto.
#
# ./getResultsSingleBorutaLOO.R results.boruta.LOO.tsv

parameters <- commandArgs(trailingOnly=TRUE)

# Archivo 'results.boruta.CV.averaged.tsv' resultado de ejecutar averageCrossValidationProbs.R
resultsBorutaAveraged <- as.character(parameters[1])

files <- system("ls meta.train.*", intern=TRUE)

ncopies <- length(files)

res <- read.table(file=resultsBorutaAveraged, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)

res <- res[which(res$decision=="Confirmed"), ]

for(i in 1:ncopies) write.table(res, file=paste("results.boruta.LOO", i, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)

