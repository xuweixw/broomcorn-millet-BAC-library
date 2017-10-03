source('E:/Bio/SequenceDatabase/R-code/screening.R')
library(Biostrings)

genomicSequence <- readDNAStringSet("E:/Bio/SequenceDatabase/homo sapiens/Homo_sapiens.GRCh38.dna.chromosome.Y.fa/Homo_sapiens.GRCh38.dna.chromosome.Y.fa")
Seq <- genomicSequence[[1]] 

baseSequence <- "GAATTC" # EcoRI: G/AATTC .
n <- matchPattern(baseSequence, Seq)
digestedSite <- start(n) + 1	#

filter <- screening(digestedSite, n=500)	# change library size, e.g. n=384*12 .

cloneSequences <- DNAStringSet(Seq, start=filter$start, end=filter$end)	# extract clone sequences from reference sequence

writeXStringSet(cloneSequences, "Inserted_sequence_of_BAC_clones.fasta")	# save result using fasta format.
