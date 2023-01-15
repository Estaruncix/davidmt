#!/usr/bin/env Rscript
###############
## ./doVolcanoPlot.R wilcox.test.Control.vs.CRC.species.tsv 1 0.01 "p.value" TRUE
## (parámetro 1) tabla con p-valores y log2FoldChanges. La tabla debe contener las columnas 'feature', 'log2FC' y 'DESCRIPTION'.
## (parámetro 2) Cutoff para el log2FoldChange
## (parámetro 3) Cutoff para el p-valor
## (parámetro 4) Nombre de la columna en la tabla que contiene los p-valores
## (parámetro 5) Si TRUE se entiende que en DESCRIPTION hay taxonomía y se asume que el phylum está en segundo lugar.
##			Si FALSE, se entiende que hay descripciones de KEGG y se toman los primeros 20 caracteres de lo que viene
##			después del primer "|", que se supone se corresponde con la primera categoría donde se clasifica el KEGG en cuestión.

parameters <- commandArgs(trailingOnly=TRUE)

mytb <- as.character(parameters[1])
myfccutoff <- as.numeric(as.character(parameters[2]))
mypcutoff <- as.numeric(as.character(parameters[3]))
pvariable <- as.character(parameters[4])
# si TRUE se entiende que las variables son taxonómicas
asTaxa <- as.logical(as.character(parameters[5]))

library(EnhancedVolcano)

MYPALLETE <- c("red", "blue", "green", "gold4", "darkviolet",
               "orange", "cyan3", "brown", "cornflowerblue", "coral", 
               "yellowgreen", "plum", "wheat", "black", "chocolate1", 
               "khaki1", "rosybrown4", "darkseagreen1", "cadetblue1", "palevioletred4")

getVolcanoColors <- function(PROBS, LOGFCS, GROUPS, P_CUTOFF, FC_CUTOFF){

  groups.colors <- rep("gray", length(GROUPS))
  names(groups.colors) <- rep("Not significant", length(GROUPS))
  i <- which(PROBS <= P_CUTOFF)
  j <- which(abs(LOGFCS) >= FC_CUTOFF)
  k <- intersect(i,j)
  # Significants
  if(length(k) > 0)
  {
    unique.groups <- unique(GROUPS[k])
    if(length(unique.groups) > length(MYPALLETE))
      stop("Too many unique groups")
    unique.groups.colors <- MYPALLETE[1:length(unique.groups)]
    names(unique.groups.colors) <- unique.groups
    groups.colors[k] <- unique.groups.colors[GROUPS[k]]
    names(groups.colors)[k] <- GROUPS[k]
  }
  return(list(groups.colors=groups.colors, lk=length(k)))

}

mygroups <- c()

df <- read.table(mytb, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE, stringsAsFactors=FALSE)
tag <- sub(".tsv", "", mytb)

pdf(file=paste("VolcanoPlot", tag, "pdf", sep="."))

if(asTaxa==TRUE){
	df$phylum <- sapply(strsplit(x=df$DESCRIPTION, split=";", fixed=TRUE), '[', 2)
	df$phylum <- sub("p_", "", df$phylum); df$phylum <- gsub("_", "", df$phylum)
	mygroups <- df$phylum
} else {
	df$category <- substr(sapply(strsplit(x=df$DESCRIPTION, split="|", fixed=TRUE), '[', 2), 1, 20)
	df$category[which(is.na(df$category))] <- "NA"
	mygroups <- df$category
}
    
res <- getVolcanoColors(df[,pvariable], df$log2FC, mygroups, mypcutoff, myfccutoff)
mycolors <- res$groups.colors
mycaption <- paste(nrow(df), " total features, ", res$lk, " with ", pvariable," <= ", mypcutoff, " and |log2FC| >= ", myfccutoff, sep="")

print(EnhancedVolcano(df,
                    lab = df$feature,
                    x = "log2FC",
                    y = pvariable,
                    xlim = c(min(df$log2FC, na.rm=TRUE), max(df$log2FC, na.rm=TRUE) + 2),
                    ylim = c(0, max(-log10(df[,pvariable]), na.rm=TRUE)),
                    xlab = "log2 FC",
                    ylab = paste("-log10",pvariable),
                    title = mytb,
                    caption = mycaption,
                    titleLabSize = 14,
                    captionLabSize = 10,
                    legendLabSize = 8,
                    legendIconSize = 2,
                    axisLabSize = 14,
                    pCutoff = mypcutoff,
                    FCcutoff = myfccutoff,
                    legendPosition = "right",
                    colCustom = mycolors,
                    transcriptPointSize = 2,
                    gridlines.minor = FALSE,
                    gridlines.major = FALSE,
                    colAlpha = 1.0))

dev.off()

