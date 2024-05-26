# Codigo

#script para dada 2 
library(dada2)
direccion_fastq <- "C:/Users/Haus/Documents/fastq_PF/"
list.files(direccion_fastq)

fnFs <- sort(list.files(direccion_fastq, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(direccion_fastq, pattern="_2.fastq", full.names = TRUE))

nombres_forward <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
nombres_reverso <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

pdf("calidad_phred/calidad_phred_reverse.pdf",width=13,height = 8)
plotQualityProfile(fnFs[1:96])
dev.off()

