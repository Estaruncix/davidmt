#!/scratch/software/R_4.2.1/bin/Rscript
#/scratch/software/Rstats/bin/Rscript
#/usr/bin/env Rscript
###############
# Se pueden pasar más de 2 grupos, pero para cada par de grupos posible se hará la transformación ANCOMBC específica de esos 2 grupos.
# Si se quiere tener una transformación ANCOMBC que incluya a más de 2 grupos usar doANCOMBC_MT2G.R.
#
# La matriz de entrada ha de ser una matriz de frecuencias absolutas.
#
# Si el último parámetro es TRUE, se identifican ceros estructurales (casi todo ceros en un grupo)
# y se les asigna por defecto un p-valor y p-valor ajustado de cero. Si es FALSE, se calcula el p-valor de esos structural zeros
# como si fueran cualquier otra variable, obteniéndose un valor muy cercano a cero pero no cero.
# Ellos en el paper intentan vender eso de que los structural zeros deberían ser automáticamente significativos (de ahí lo de p-valor cero),
# pero en la función de R el tratamiento por defecto es calcular su p-valor como si se tratara de cualquier otra variable.
# Parece como si eso de dar un p-valor de cero a los structural zeros fuera una decisión que asumieron en versiones anteriores
# de la herramienta y que ahora se hayan arrepentido, aunque sigan dando la opción de hacerlo por coherencia con sus métodos antiguos.
#
# El penúltimo parámetro, si es un estudio pareado ha de valer "group + donor". Si no, tenemos la opción de pasar solo "group" o "group + whatever" si
# quisiéramos incorporar el factor 'whatever' como un efecto aleatorio a la hora de analizar las diferencias en base a "group".
#
# IMPORTANTE: El nombre 'whatever' ha de ser distinto de 'donor'. El comportamiento por defecto del script para muestras pareadas es incorporando la
#	variable 'donor' como un efecto fijo sin interacción con 'group' (interceptación variable pero slope fija). Si quisiéramos tratar el análisis
#	pareado desde la óptica de los random effects podemos hacerlo duplicando la variable 'donor' en el metadata usando otro nombre 'whatever' y pasando
#	entonces la fórmula "group + whatever".
#
# El manejo de factores aleatorios (modelos mixtos con librerías nlme y lme4) se incorpora en la versión 1.6.4 de ANCOMBC que necesita R versión >= 4.2.0.
#
# ./doANCOMBC.R counts.asv.tsv meta.asv.tsv "GVHDs,week" "no.pre,yes.pre" "sinfilt" FALSE 0 0 TAXONOMY "group" FALSE

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
myfactors <- as.character(parameters[3])
mylevels <- as.character(parameters[4])
study <- as.character(parameters[5])
doFilt <- as.logical(as.character(parameters[6]))
mypct <- as.numeric(as.character(parameters[7]))
myfac <- as.numeric(as.character(parameters[8]))
descr.column <- as.character(parameters[9])
myformula <- as.character(parameters[10])
structzero <- as.logical(as.character(parameters[11]))

if(length(grep("donor", myformula))>0) mypaired <- TRUE else mypaired <- FALSE

if ("myfunctions.downstream.R" %in% list.dirs()){
  source("myfunctions.downstream.R")
  
} else{
  source("/mnt/c/Users/RositaRos/Desktop/analysis/estudio/scripts/myfunctions.downstream.R")
}


library(ANCOMBC)
library(phyloseq)

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
origIDs <- dat[, 1]
dat[, 1] <- gsub(";", ".", dat[, 1])
dat[, 1] <- gsub("_", ".", dat[, 1])
dat[, 1] <- gsub("-", ".", dat[, 1])
dat[, 1] <- gsub("[(]", ".", dat[, 1])
dat[, 1] <- gsub("[)]", ".", dat[, 1])
dat[, 1] <- gsub("[[]", ".", dat[, 1])
dat[, 1] <- gsub("[]]", ".", dat[, 1])
dat[, 1] <- gsub(" ", ".", dat[, 1])
dat[, 1] <- gsub(":", ".", dat[, 1])
dat[, 1] <- gsub("/", ".", dat[, 1])
dat[, 1] <- gsub("[|]", ".", dat[, 1])
names(origIDs) <- dat[, 1]

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

if(mypaired==TRUE){
	meta$donor <- as.character(meta$donor)
	myfix_formula <- myformula
	myrand_formula <- NULL
} else {
	myrandvar <- setdiff(unlist(strsplit(myformula, split=" [+] ")), "group")
	myfix_formula <- "group"
	if(length(myrandvar)>0){
		myrand_formula <- paste("(1 | ", myrandvar, ")", sep="")
	} else {
		myrand_formula <- NULL
	}
}

dat.pct <- prop.table(dat, margin=2)*100

dat.nofilt <- dat

if(doFilt==TRUE){
	mycut <- min(abs(dat.pct[which(dat!=0)]))*myfac
	good <- doFiltering(dat.pct, cutoff=mycut, group=GROUPS, mypct=mypct)
	dat <- dat[good, ]; dat.pct <- dat.pct[good, ]
	mydata <- as.data.frame(dat)
	mydata$ID <- rownames(mydata)
	mydata <- mydata[, c("ID", colnames(dat))]
	write.table(mydata, file=paste("filtered.rel.freqs", counts.file, sep="."), sep="\t", row.names=FALSE, quote=FALSE)
	cond <- apply(dat, MARGIN=2, sum)
	dat <- dat[, which(cond>0)]; dat.pct <- dat.pct[, which(cond>0)]
	meta <- meta[which(cond>0), ]; GROUPS <- GROUPS[which(cond>0)]
} else { good <- rownames(dat) }

mycombn <- combn(mylevels,2)

for(i in 1:ncol(mycombn)){

g1 <- mycombn[1,i]
g2 <- mycombn[2,i]

i1 <- which(GROUPS==g1)
i2 <- which(GROUPS==g2)
if(mypaired==TRUE){
	# Esta ordenación ya no es necesaria (lo era para usar t.test o wilcox.test).
	names(i1) <- meta$donor[i1]
	names(i2) <- meta$donor[i2]
	mydonors <- intersect(names(i1), names(i2))
	i1 <- i1[mydonors]; i2 <- i2[mydonors]
}
ind <- c(i1, i2)
tag <- paste(g1, "vs", g2, sep=".")

meta.ancom <- meta[ind, ]
meta.ancom$group <- GROUPS[ind]
OTU <- otu_table(dat.nofilt[, ind], taxa_are_rows=TRUE)
META <- sample_data(meta.ancom)
physeq <- phyloseq(OTU, META)

#### Running ancombc2 function.
# If the group variable has < 3 categories the multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated.
# res:
#	columns started with ‘p’: p-values. P-values are obtained from two-sided Z-test using the test statistic ‘W’.
#	columns started with ‘q’: adjusted p-values.  Adjusted p-values are obtained by applying ‘p_adj_method’ to ‘p’.
out <- ancombc2(data=physeq,
			fix_formula=myfix_formula, rand_formula=myrand_formula,
              p_adj_method="fdr", lib_cut=0, prv_cut = 0.05,
              group="group", struc_zero=structzero, neg_lb=FALSE, 
              global=FALSE)
              
# Aparece este warning más de 50 veces:
# In lme4::lmer(formula = tformula, data = df, control = lme_control) :
#  reiniciar evaluación premisa interrumpida
#
# ... y los p-valores salen todos a 1 ... algo no va bien y no encuentro información sobre este warning en la red ..

res <- out$res

#### Bias-adjusted abundances
samp_frac <- out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] <- 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn <- log(out$feature_table + 1) 
# Adjust the log observed abundances
dat.trans <- t(t(log_obs_abn) - samp_frac)

# filtramos después de transformar
dat.trans <- dat.trans[intersect(good, rownames(dat.trans)), ]

if(structzero==TRUE){
	structural_zeros <- as.data.frame(out$zero_ind)
	structural_zeros$DESCRIPTION <- DESCRIPTION[rownames(structural_zeros)]
	structural_zeros$ID <- rownames(structural_zeros)
	structural_zeros <- structural_zeros[, c("ID", setdiff(colnames(structural_zeros), c("ID", "DESCRIPTION")), "DESCRIPTION")]
	cond1 <- apply(structural_zeros[, c(2,3)], MARGIN=1, FUN=function(x) if(all(x==TRUE)) return(TRUE) else return(FALSE))
	cond2 <- apply(structural_zeros[, c(2,3)], MARGIN=1, FUN=function(x) if(any(x==TRUE)) return(TRUE) else return(FALSE))
	structural_zeros <- rbind(structural_zeros[which(cond1==TRUE), ], structural_zeros[which(cond2==TRUE),], structural_zeros[which(cond2==FALSE),])
	write.table(structural_zeros, file=paste("structural_zeros", tag, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
}

indpval <- grep("p_", colnames(res))
tem <- grep("p_group", colnames(res))
if(length(indpval)>1) indpval <- tem
indadjpval <- grep("q_", colnames(res))
tem <- grep("q_group", colnames(res))
if(length(indadjpval)>1) indadjpval <- tem

df <- data.frame(feature=rownames(out$feature_table), ANCOMBC2.pval=res[, indpval], ANCOMBC2.adjpval.fdr=res[, indadjpval], stringsAsFactors=FALSE)
df$DESCRIPTION <- DESCRIPTION[df$feature]
if(structzero==TRUE) df <- cbind(df, structural_zeros[df$feature, c(2,3)])

mytb <- differential.tests(MATRIX=dat.pct[, ind], GROUPS=GROUPS[ind], GROUP1=g1, GROUP2=g2, TEST="wilcox.test",
						PAIRED=mypaired, CUTOFF=0.1, COLNAME="feature", FILE=paste("ANCOMBC2.test", tag, study, sep="."), DESCRIPTION=DESCRIPTION)
mytb$feature <- as.character(mytb$feature)
rownames(mytb) <- mytb$feature
mytb <- mytb[, c("mean1","mean2","log2FC")]

# Filtramos tabla de p-valores de ANCOMBC y recalculamos p-valor ajustado
rownames(df) <- df$feature
df <- df[intersect(good, rownames(df)), ]
df[, "ANCOMBC2.adjpval.fdr"] <- p.adjust(df[, "ANCOMBC2.pval"], method="fdr")
myord <- order(df$ANCOMBC2.pval, decreasing=FALSE)
df <- df[myord, ]

mytb <- mytb[rownames(df), ]
df <- cbind(df, mytb)
df <- df[, c(setdiff(colnames(df), "DESCRIPTION"), "DESCRIPTION")]

myfeatures <- as.character(df[which(df$ANCOMBC2.pval<=0.1), "feature"])

# pintamos boxplots con datos transformados (machacaremos la gráfica de differential.tests)
G1 <- names(GROUPS)[which(GROUPS==g1)]
G2 <- names(GROUPS)[which(GROUPS==g2)]
items.group1 <- which(is.element(colnames(dat.trans), G1))
items.group2 <- which(is.element(colnames(dat.trans), G2))

pdf(file=paste("ANCOMBC2.test", tag, study, "pdf", sep="."))
for(sign.row in myfeatures){
    	values <- list(dat.trans[sign.row, items.group1], dat.trans[sign.row, items.group2])
    	names(values) <- c(g1, g2)
    	mymain <- DESCRIPTION[sign.row]
    	beeswarm(values, pch=16, col=c("red","blue"), las=1, xlab="Condition", ylab=sign.row, main=mymain, cex.lab=0.6, cex.main=0.5)
    	graphics::boxplot(values, las=1,  add=TRUE, col=NULL)
    	legend("topright", legend=c(paste("ANCOMBC2.pval:"         , format(df[sign.row, "ANCOMBC2.pval"    ], digits=2)),
                                paste("ANCOMBC2.adjpval.fdr:", format(df[sign.row, "ANCOMBC2.adjpval.fdr"], digits=2))), cex=0.7)
}
dev.off()
	
if(mypaired==TRUE){
	pdf(file=paste("ANCOMBC2.test.paired.with.lines", tag, study, "pdf", sep="."))
	for(myfeature in myfeatures){
		plotPaired(dat=dat.trans, feature=myfeature, group=GROUPS[colnames(dat.trans)],
				group.x=g1, group.y=g2, individual.factor=meta[colnames(dat.trans), "donor"],
				DESCRIPTION=DESCRIPTION, p.val=df[myfeature, "ANCOMBC2.pval"])
	}
	dev.off()
}

df$feature <- origIDs[df$feature]
write.table(df, file=paste("ANCOMBC2.test", tag, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)

dat.trans <- as.data.frame(dat.trans)
dat.trans$ID <- origIDs[rownames(dat.trans)]
dat.trans <- dat.trans[, c("ID", setdiff(colnames(dat.trans), "ID"))]
dat.trans$DESCRIPTION <- DESCRIPTION[rownames(dat.trans)]
write.table(dat.trans, file=paste("ANCOMBC2.normalized", counts.file, sep="."), sep="\t", row.names=FALSE, quote=FALSE)

# generamos input boruta
tem <- data.frame(feature=df$feature, p.value=rep(0, nrow(df)), adj.p.value=rep(0, nrow(df)), DESCRIPTION=df$DESCRIPTION, stringsAsFactors=FALSE)
write.table(tem, file=paste("input.boruta", counts.file, sep="-"), sep="\t", row.names=FALSE, quote=FALSE)

print(tag)

}

