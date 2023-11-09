
# scripts to run composition analysis using as input protein and nucleotide fasta sequences

# Run analysis using the R environment

# upload library
library(seqinr)

############# calculate aminoacid frequencies and pI (theoretical isoelectric point)

# Get the list of all *.txt files in the current directory
prot_files <- list.files(pattern = "_prot.faa")

# Loop through the list of files
for (file in prot_files) {
  # Read the file into a variable
  myProts <- read.fasta(file = file, seqtype = "AA")

  # Get stats from all amino acid sequences
  all_results <- lapply(myProts, AAstat, plot=FALSE)

  # Get amino acid proportions
  all_proportions <- lapply(all_results, function(x) x$Prop)

  table <- do.call(rbind, all_proportions)

  # Print composition of all proteins
  all_composition <- lapply(all_results, function(x) x$Compo)

  table2 <- do.call(rbind, all_composition)

  # Get protein length
  lengths <- sapply(myProts, getLength)

  # Calculate the pI (theoretical isoelectric point) for each protein sequence
  pIs <- sapply(myProts, computePI)

  # Combine the tables
  table_all_combined <- cbind(table,table2,lengths,pIs)

  # Write the output to a CSV file
  write.csv(table_all_combined, file = paste(file, "_aminoacid_stats.csv", sep=""))
}


############# calculate codon frequencies and codon usage bias

# upload library
library(seqinr)

# Get the list of all *.txt files in the current directory
seq_files <- list.files(pattern = "_nucl.fna")

# Loop through the list of files
for (file in seq_files) {
  # Read the file into a variable
  sequences <- read.fasta(file = file, seqtype = "DNA")

  codons_rscu <- lapply(sequences, uco, index="rscu", NA.rscu = 0)

  table_codons_rscu <- do.call(rbind, codons_rscu)

  # Write the output to a CSV file
  write.csv(table_codons_rscu, file = paste(file, "_table_codons_rscu.csv", sep=""))
}



############# calculate Purine loading-index


# upload library
library(seqinr)

seq_files <- list.files(pattern = "_nucl.fna$") # whole CDSs, make case specific $

# Loop through the list of files
for (file in seq_files) {
  # Read the file into a variable
  sequences <- read.fasta(file = file, seqtype = "DNA")
  
  A_count <- sapply(sequences, function(x) sum(toupper(x) == "A"))
  T_count <- sapply(sequences, function(x) sum(toupper(x) == "T"))
  G_count <- sapply(sequences, function(x) sum(toupper(x) == "G"))
  C_count <- sapply(sequences, function(x) sum(toupper(x) == "C"))
  N_length <- sapply(sequences, length)
  
  purines <- A_count+G_count
  pyrimidines <- C_count+T_count
  
  dW <- (A_count - T_count) / N_length *1000
  dS <- (G_count - C_count) / N_length *1000

  PLI_purine_loading_index <- dW + dS
  
  
  purr <- cbind(A_count,T_count,G_count,C_count,purines,pyrimidines,PLI_purine_loading_index)

  # Write the output to a CSV file
  write.csv(purr, file = paste(file, "_purine_loading_index_PLI.csv", sep=""))
}

