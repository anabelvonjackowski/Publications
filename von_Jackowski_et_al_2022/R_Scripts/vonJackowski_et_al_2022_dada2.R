
##############################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## dada2 Analysis
##############################################################################

# Open screen
screen -S pebcao (81795.pebcao)
Ctrl+a+d (detach)
screen -r 81795.pebcao (reattach)
screen -x 81795.pebcao (restore if already attached)

# Start R server session
module load R 
R
require(dada2)
require(ShortRead)
require(ggplot2)
require(gridExtra)

# Setwd and load
setwd("/scratch1/mwietz/PEBCAO")
load("PEBCAO.Rdata")

# List files
path <- "/scratch1/mwietz/PEBCAO/Clipped/"
fns <- list.files(path)
fns

# Sort: forward/reverse reads in same order
fnFs <- sort(list.files(path, pattern="_R1.fastq"))
fnRs <- sort(list.files(path, pattern="_R2.fastq"))

# Provide names of clipped files
sample.names <- sort(read.table(
  "/scratch1/mwietz/PEBCAO/sample_names.txt", 
  h = F, stringsAsFactors = F)$V1)

# Specify full path to fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#################################

# Quality check
QualityProfileFs <- list()
for(i in 1:length(fnFs)) {QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(fnFs[i])}
pdf("QualityProfileForward.pdf")
for(i in 1:length(fnFs)) {  do.call("grid.arrange", QualityProfileFs[[i]])}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(fnRs)) {QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(fnRs[i])}
pdf("QualityProfileReverse.pdf")
for(i in 1:length(fnRs)) {  do.call("grid.arrange", QualityProfileRs[[i]])}
dev.off()
rm(QualityProfileRs)
# expected max length: 400 // min overlap: 30

# mkdir and filenames for filtered fastqs
filt_path <- file.path(path, "../Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

#################################

# Filter: settings based on quality profiles
out <- filterAndTrim(
  fnFs, 
  filtFs, 
  fnRs, 
  filtRs,
  truncLen = c(230, 195),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 6)

# should be retaining >80% (0.8) OK here!
head(out)
summary(out[, 2]/out[, 1])

#################################

# Quality check 
QualityProfileFs.filt <- list()
for(i in 1:length(filtFs)) {
  QualityProfileFs.filt[[i]] <- list()
  QualityProfileFs.filt[[i]][[1]] <- plotQualityProfile(filtFs[i])}
pdf("QualityProfileForwardFiltered.pdf")
for(i in 1:length(filtFs)) {do.call("grid.arrange", QualityProfileFs.filt[[i]])}
dev.off()
rm(QualityProfileFs.filt)

QualityProfileRs.filt <- list()
for(i in 1:length(filtRs)) {
  QualityProfileRs.filt[[i]] <- list()
  QualityProfileRs.filt[[i]][[1]] <- plotQualityProfile(filtRs[i])}
pdf("QualityProfileReverseFiltered.pdf")
for(i in 1:length(filtRs)) {  do.call("grid.arrange", QualityProfileRs.filt[[i]])}
dev.off()
rm(QualityProfileRs.filt)

#################################

# Learn errors 
errF <- learnErrors(filtFs, multithread = 6, randomize = TRUE, verbose = 1, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread = 6, randomize = TRUE, verbose = 1, MAX_CONSIST = 20)
# convergence after 5-7 rounds -- good!

# Plot error profiles
pdf("ErrorProfiles.pdf")
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()
# should be only few points outside the black line - but still ok here

# Dereplication 
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name objects 
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Denoising
dadaFs <- dada(derepFs, err = errF, multithread = 6, pool = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = 6, pool = TRUE)

#################################

# Read merging
mergers <- mergePairs(
  dadaFs, 
  derepFs, 
  dadaRs,
  derepRs,
  minOverlap = 10,
  verbose = TRUE,
  propagateCol = c("birth_fold", "birth_ham"))

# Create sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  # 11451 sequences
saveRDS(seqtab, "/scratch1/mwietz/Particles_BAH/renamed/seqtab.rds")

# Remove chimeras // Identified 7501 bimeras out of 11451 
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 6, verbose = TRUE)

dim(seqtab.nochim)  # 3950 sequences
summary(rowSums(seqtab.nochim)/rowSums(seqtab))

#################################

# Remove singletons and 'junk' sequences
# Show distribution of read length and select for filtering
table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim)))
seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% 
  c(354:412) & colSums(seqtab.nochim) > 1]

# Inspect output
dim(seqtab.nochim2)
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))
summary(rowSums(seqtab.nochim2))

# Get summaries 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim2))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nochim", "tabled")
rownames(track) <- sample.names
track <- data.frame(track)
head(track)

# inspect output -- looks ok!
summary(track$tabled/track$input) 
summary(track$filtered/track$input)
summary(track$denoised/track$filtered)
summary(track$merged/track$denoised)
summary(track$nochim/track$merged)
summary(track$tabled/track$nochim)

#################################

# Assign taxonomy
tax <- assignTaxonomy(
  seqtab.nochim2, 
  "silva_nr_v132_train_set.fa.gz", 
  tryRC = TRUE,
  multithread = 10)

## Summary -- Archaea 122 // Bacteria 3391
summary(tax) 

# Select BACTERIA AND ARCHAEA 
# remove OTUs unclassified on phylum level 
# Remove EUK sequences
table(tax[, 1])   
sum(is.na(tax[, 2]))   # result: 116
tmp.good <- tax[!is.na(tax[, 2]) & tax[, 1] %in% c("Bacteria", "Archaea"),]
tax.good <- tmp.good[-c(grep("Chloroplast", tmp.good[, 4]), grep("Mitochondria", tmp.good[, 5])), ]
seqtab.nochim2.good <- seqtab.nochim2[, rownames(tax.good)]
summary(rowSums(seqtab.nochim2.good))

# Format output
seqtab.nochim2.print <- t(seqtab.nochim2.good)
tax.print <- tax.good
all.equal(rownames(seqtab.nochim2.print), rownames(tax.print))  #TRUE
rownames(seqtab.nochim2.print) <- paste("sq", 1:ncol(seqtab.nochim2.good), sep = "")
rownames(tax.print) <- rownames(seqtab.nochim2.print)

# Write output
write.table(seqtab.nochim2.print, "bac_seqtab.txt", quote = F, sep = "\t")
write.table(tax.print, "bac_tax.txt", sep = "\t", quote = F)
write.table(track, file="dada_summary.txt", sep="\t", row.names=T)
uniquesToFasta(seqtab.nochim2.good, "bac_uniques.fasta")

#################################

save.image("PEBCAO.Rdata")
