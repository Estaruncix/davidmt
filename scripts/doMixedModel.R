#!/scratch/software/Rstats/bin/Rscript
#/usr/bin/env Rscript
###############
# Calcula un p-valor para cada una de las variables de la matriz de entrada (variables no taxonómicas). El modelo para obtener el p-valor es un modelo mixto:
#
# v ~ F + R
#
# donde v es una variable dada la matriz de entrada, F contiene los grupos a comparar (efecto fijo) y R codifica
# las distintas observaciones sobre el mismo individuo (efecto aleatorio).
#
# F y R deben ser nombres de columna definidos en el metadata.
#
# El tratamiento del efecto fijo se hace de dos formas
#	(1) interceptación variable y slope fija para cada individuo
#	(2) interceptación variable y slope variable para cada individuo
#
# La opción (2) será de mayor precisión que la (1) cuando tengamos suficientes observaciones en cada individuo (al menos 3).
# Pero el script devuelve 2 p-valores, uno a partir de (1) (pvalslopefix) y otro a partir de (2) (pvalslopevar).
#
# En cualquiera de los 2 casos, este modelo es una alternativa de mayor precisión comparada con un simple agregado en el caso de tener réplicas
# o distintos tiempos sobre un mismo paciente.
#
# Aplicamos la función de R lme (modelos mixtos con librería nlme).
#
# ./doMixedModel.R diversity.tsv metadata_pareado.csv "diet" "AltaADM,BajaADM" "mixedModel" "ID" "diet + donor"

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
myfactors <- as.character(parameters[3])
mylevels <- as.character(parameters[4])
study <- as.character(parameters[5])
descr.column <- as.character(parameters[6])
myformula <- as.character(parameters[7])

source("myfunctions.downstream.R")

library(nlme)

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
DESCRIPTION <- dat[, descr.column]; names(DESCRIPTION) <- rownames(dat) <- dat[, 1]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

myfactors.l <- strsplit(myfactors, split=",")
myfactor1 <- myfactors.l[[1]][1]
if(length(myfactors.l[[1]])==2) myfactor2 <- myfactors.l[[1]][2] else myfactor2 <- ""

if(myfactor2!="") GROUPS <- paste(meta[, myfactor1], meta[, myfactor2], sep=".") else GROUPS <- as.character(meta[, myfactor1])
names(GROUPS) <- mysamples

mylevels <- unlist(strsplit(mylevels, split=","))

mylevels <- mylevels[which(is.element(mylevels, unique(GROUPS)))]

ind <- which(is.element(GROUPS, mylevels))

dat <- dat[, ind]
GROUPS <- GROUPS[ind]
meta <- meta[ind, ]

# Eliminamos variables a todo ceros
cond <- apply(dat, MARGIN=1, FUN=function(x) if(all(x==0)) return(FALSE) else return(TRUE))
dat <- dat[which(cond==TRUE), ]

tem <- unlist(strsplit(myformula, split=" [+] "))
myrandvar <- tem[2]
myfixvar <- tem[1]
myrand_formula_slopefix <- paste("~1 | ", "factor(", myrandvar, ")", sep="")
myrand_formula_slopevar <- paste("~", myfixvar, " | factor(", myrandvar, ")", sep="")

g1 <- mylevels[1]
g2 <- mylevels[2]

i1 <- which(GROUPS==g1)
i2 <- which(GROUPS==g2)

ind <- c(i1, i2)
tag <- paste(g1, "vs", g2, sep=".")

#####
# Aplicamos lme (modelos mixtos con librería nlme) a la matriz de datos:

pvalslopefix <- rep(NA, nrow(dat)); names(pvalslopefix) <- rownames(dat); adjpvalslopefix <- rep(NA, nrow(dat))
pvalslopevar <- rep(NA, nrow(dat)); names(pvalslopevar) <- rownames(dat); adjpvalslopevar <- rep(NA, nrow(dat))

myoption <- lmeControl(opt='optim')

for(myvar in rownames(dat)){

	tryCatch({

	mydat <- data.frame(group=meta[, myfixvar], stringsAsFactors=FALSE)
	colnames(mydat) <- myfixvar
	mydat[, myvar] <- dat[myvar, ]; mydat[, myrandvar] <- meta[, myrandvar]
	myform.fixed <- as.formula(paste(myvar, " ~ ", myfixvar, sep=""))
	myform.random <- as.formula(myrand_formula_slopefix)
	lme.res <- lme(myform.fixed, data=mydat, random = myform.random, control=myoption)
	pvalslopefix[myvar] <- anova(lme.res)[[4]][2]

	}, error=function(e){})

}

for(myvar in rownames(dat)){

	tryCatch({

	mydat <- data.frame(group=meta[, myfixvar], stringsAsFactors=FALSE)
	colnames(mydat) <- myfixvar
	mydat[, myvar] <- dat[myvar, ]; mydat[, myrandvar] <- meta[, myrandvar]
	myform.fixed <- as.formula(paste(myvar, " ~ ", myfixvar, sep=""))
	myform.random <- as.formula(myrand_formula_slopevar)
	lme.res <- lme(myform.fixed, data=mydat, random = myform.random, control=myoption)
	pvalslopevar[myvar] <- anova(lme.res)[[4]][2]

	}, error=function(e){})

}

#####

df <- data.frame(feature=rownames(dat),
			pvalslopefix=pvalslopefix,
			adjpvalslopefix=adjpvalslopefix,
			pvalslopevar=pvalslopevar,
			adjpvalslopevar=adjpvalslopevar,
			stringsAsFactors=FALSE)
df$DESCRIPTION <- DESCRIPTION[df$feature]

mytb <- differential.tests(MATRIX=dat[, ind], GROUPS=GROUPS[ind], GROUP1=g1, GROUP2=g2, TEST="wilcox.test",
						PAIRED=FALSE, CUTOFF=1, COLNAME="feature", FILE=paste("wilcox.test", tag, study, sep="."), DESCRIPTION=DESCRIPTION)
mytb$feature <- as.character(mytb$feature)
rownames(mytb) <- mytb$feature

rownames(df) <- df$feature
df[, "adjpvalslopefix"] <- p.adjust(df[, "pvalslopefix"], method="fdr")
df[, "adjpvalslopevar"] <- p.adjust(df[, "pvalslopevar"], method="fdr")
myord <- order(df$pvalslopefix, decreasing=FALSE)
df <- df[myord, ]

mytb <- mytb[rownames(df), ]
df <- cbind(df, mytb)
df <- df[, c(setdiff(colnames(df), "DESCRIPTION"), "DESCRIPTION")]

write.table(df, file=paste("mixedmodel", tag, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)


