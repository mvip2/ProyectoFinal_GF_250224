# Codigo

#script para dada 2 
library(dada2)
direccion_fastq <- "C:/Users/Haus/Documents/fastq_PF/"#esta es la carpeta local donde están las secuencias fastq, es necerio modificar la dirección dependiendo de la ubicación, ya que no pudimos subir las secuencias a github (2.5 GB)
list.files(direccion_fastq) # nos muestra los archivos que están en la carpeta

fnFs <- sort(list.files(direccion_fastq, pattern="_1.fastq", full.names = TRUE)) #seleccionamos los archivos forward con el patron _1.fastq
fnRs <- sort(list.files(direccion_fastq, pattern="_2.fastq", full.names = TRUE)) # lo mismo de arriba pero con los reversos 

nombres <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


pdf("calidad_phred/calidad_phred_forward.pdf",width=13,height = 8)
plotQualityProfile(fnFs[1:3])
dev.off() # la calidad phred es buena, es cercana a 40, mayor a 30

pdf("calidad_phred/calidad_phred_reverse.pdf",width=13,height = 8)
plotQualityProfile(fnRs[1:3])
dev.off() # La calidad phred también 

filtFs <- file.path(direccion_fastq, "filtered", paste0(nombres, "_F_filt.fastq.gz"))
filtRs <- file.path(direccion_fastq, "filtered", paste0(nombres, "_R_filt.fastq.gz"))
names(filtFs) <- nombres
names(filtRs) <- nombres

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), # aqui vamos a filtrar los datos, es decir, vamos a comprimir los archivos fastq y además vamos a cortar las secuencias con calidad phred inferior a 30
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
# los filtrados fueron subidos a la misma carpeta de drive donde subimos los fastq

errF <- learnErrors(filtFs, multithread = FALSE)
#104908560 total bases in 437119 reads from 20 samples will be used for learning the error rates.

errR <- learnErrors(filtRs, multithread = FALSE)
#100056160 total bases in 625351 reads from 42 samples will be used for learning the error rates.

plotErrors(errF, nominalQ = TRUE)
#  log-10 transformation introduced infinite values.
plotErrors(errR, nominalQ = TRUE)
#  log-10 transformation introduced infinite values.

dadaFs <- dada(filtFs, err = errF, multithread = FALSE) # como estamos usando windows, voy a cambiar el valor de multithread por FALSE
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#[1]   96 7370

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

table(nchar(getSequences(seqtab.nochim))) #remover quimeras

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- nombres
head(track)

library("devtools")
devtools::install_github("benjjneb/dada2")
library(Rcpp)
library(dada2)

taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/mvip2/Desktop/260524/silva_nr_v132_train_set.fa.gz", multithread=F) #la dirección de la base de datos está dada de forma local
taxa <- addSpecies(taxa, "C:/Users/mvip2/Desktop/260524/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# más codigo
library(phyloseq)
library(Biostrings)
library(ggplot2)


