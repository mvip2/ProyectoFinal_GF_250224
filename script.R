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

taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Haus/Downloads/RDP_LSU_fixed_train_set_v2.fa.gz", multithread=FALSE) # utilizamos el RDP fungi LSU trainset 11 en lugar de UNITE
taxa <- addSpecies(taxa, "C:/Users/Haus/Downloads/rdp_species_assignment_LSU_v2.fa.gz")
saveRDS(taxa, file="taxa.RDS")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# más codigo
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(data.table)
library(tidyverse)

# a partir de aqui crearemos el objeto phyloseq, utilizamos el tutorial de dada2, sin embargo, no pudimos colocar los metadatos correctos
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
save(ps,file="ps.RDS")

readRDS(file="taxa.RDS")

otu_table(ps)
sample_data(ps)
tax_table(ps)

atlas.prune = prune_taxa(taxa_sums(ps) > 1, ps)

readcount <- data.table(as(sample_data(atlas.prune), "data.frame"),
                        TotalReads = sample_sums(atlas.prune), 
                        keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")
ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")


head(readcount[order(readcount$TotalReads), c("SampleID", "TotalReads")])


otu.rare <- otu_table(atlas.prune)
otu.rare <- as.data.frame(t(otu.rare))
sample_names <- rownames(otu.rare)

library(vegan)

otu.rarecurve = rarecurve(otu.rare, step = 10000)

otu.rarecurve = rarecurve(otu.rare[1:20], step = 10000,label=FALSE)

plot_richness(atlas.prune, x="taxa",measures=c("Chao1","Simpson","Shannon")) + geom_boxplot()
