

#### PACKAGES ####

if(!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("sangerseqR")) {
  BiocManager::install("sangerseqR")
}

if (!require("DECIPHER")) {
  BiocManager::install("DECIPHER")
}

packages_used <-
  c(
    "tidyverse",
    "ape",
    "janitor",
    "sangerseqR",
    "Biostrings",
    "DECIPHER",
    "msa",
    "pegas"
  )

packages_to_install <-
  packages_used[!packages_used %in% installed.packages()[, 1]]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install,
                   Ncpus = parallel::detectCores() - 1
  )
}

lapply(packages_used,
       require,
       character.only = TRUE
)

#### rename Sequences from BLaST MSA OUTPUT ####

# rename the fasta file downloaded from ncbi blast msa
# assumes you have run wrangle_blast_output.bash to create tidy tsv of information int he seq names

renameBlastFastaSeqs <-
  function(
    inTsvFile = "../output/sanger_curated_ab1_ischnura_luta/blast_rbd_06_E1_500_better.tsv",
    inFastaFile = "../output/sanger_curated_ab1_ischnura_luta/blast_rbd_06_E1_500.fasta",
    outFastaFile = "../output/sanger_curated_ab1_ischnura_luta/blast_rbd_06_E1_500_renamed.fasta"
  ){
    
    # read in tsv with info from seq names
    data_blast_names <-
      read_tsv(inTsvFile) %>%
      clean_names() %>%
      mutate(
        accession = str_remove(accession,
                               "\\..*$"),
        organism = 
          str_replace_all(organism,
                          "_sp._[A-Z].*$",
                          "_sp.") %>%
          str_replace_all(.,
                          "_",
                          "_"),
        country = str_remove(country,
                             ":.*$"),
        country = replace_na(country, 
                             "Unk."),
        seq_name = str_c(
          organism,
          country,
          accession,
          sep = ", "
        ) 
      )
    
    # read in fasta from blast msa
    # blast_seqs <- read.fasta(inFastaFile)
    blast_seqs <- ape::read.dna(inFastaFile, format = "fasta")
    
    # filter out seqs with "Query" in seq name
    
    # blast_seqs <- 
    #   read.fasta(
    #     inFastaFile
    #   )[!sapply(names(blast_seqs), 
    #             function(x) grepl("Query", 
    #                               x))]
    rows_with_query <- grepl("Query", dimnames(blast_seqs)[[1]])
    blast_seqs_filtered <- blast_seqs[!rows_with_query, ]
    
    # replace original fasta seq names with those from tsv
    dimnames(blast_seqs_filtered)[[1]] <- 
      data_blast_names$seq_name %>% 
      head(
        n=length(
          dimnames(blast_seqs_filtered)[[1]]
        )
      )
    
    # write renamed fasta file
    write.dna(blast_seqs_filtered, 
              file = outFastaFile, 
              format = "fasta")
    
    # remove carriage returns
    
    fasta_fixed <- 
      outFastaFile %>%
      map_chr(readr::read_file) %>%
      str_remove_all(.,"\r") %>%
      str_remove_all(.," ") %>%
      paste(collapse = "\n") %>%
      str_to_upper() # %>%
    # str_remove(.," ")
    
    # Write the combined content to a new file
    writeLines(
      fasta_fixed, 
      outFastaFile
    )
  }


#### Process Curated AB1 Files ####
# prior to running this, must have dir with curated AB1 files
# the rev ab1 should NOT be reverse and complemented
processCuratedAB1 <-
  function(
    ab1_dir = "../output/sanger_curated_ab1_ischnura_luta/",
    out_file = "../output/sanger_curated_ab1_ischnura_luta/consensus_sequences.fasta",
    fwd_primer_name = "LCO1490",
    rev_primer_name = "HCO2198"
  ){
    #### READ IN AB1 ####
    
    ab1_files <- 
      list.files(ab1_dir, 
                 pattern = "\\.ab1$", 
                 full.names = TRUE)
    
    
    
    # Initialize an empty list to store AB1 data
    ab1_data_list <- list()
    
    # Loop through each AB1 file and read its data
    for(file_path in ab1_files) {
      file_name <- basename(file_path) # Extracts the file name from the path
      ab1_data_list[[file_name]] <- sangerseqR::read.abif(file_path)
    }
    
    #### MAKE CONSENSUS SEQUENCES ####
    
    # Correctly identify and organize pairs
    file_pairs <- list()
    
    for (file in ab1_files) {
      specimen_id <- sub("(.+)_([A-Z0-9]+[0-9A-Z]+)_(.+)_(.+)\\.ab1$", "\\1", basename(file))
      primer_type <- sub("(.+)_([A-Z0-9]+[0-9A-Z]+)_(.+)_(.+)\\.ab1$", "\\2", basename(file))
      
      if (!exists(specimen_id, 
                  where = file_pairs)) {
        file_pairs[[specimen_id]] <- 
          list(fwd = character(0), 
               rev = character(0))
      }
      
      if (primer_type == fwd_primer_name) {
        file_pairs[[specimen_id]]$fwd <- 
          c(file_pairs[[specimen_id]]$fwd, file)
      } else if (primer_type == rev_primer_name) {
        file_pairs[[specimen_id]]$rev <- 
          c(file_pairs[[specimen_id]]$rev, file)
      }
    }
    
    # Filter out specimens without both fwd and rev sequences
    file_pairs_filtered <- Filter(function(x) length(x$fwd) == 1 && length(x$rev) == 1, file_pairs)
    
    consensus_sequences <- list()
    
    for (specimen in names(file_pairs_filtered)) {
      fwd_file <- file_pairs_filtered[[specimen]]$fwd
      rev_file <- file_pairs_filtered[[specimen]]$rev
      
      # Read sequences
      fwd_seq <- sangerseqR::read.abif(fwd_file[[1]])@data$PBAS.1
      rev_seq <- 
        sangerseqR::read.abif(rev_file[[1]])@data$PBAS.1 %>%
        DNAString() %>%
        reverseComplement() %>% 
        as.character()
      
      # First, combine them into a DNAStringSet
      sequences <- 
        DNAStringSet(
          c(
            fwd_seq, 
            rev_seq
          )
        )
      
      # Perform the alignment
      alignment <- msa(sequences,
                       method = "ClustalOmega",
                       type = "dna",
                       # gapOpening = 100,
                       # gapExtension = 100,
                       verbose = TRUE)
      
      
      # Calculate the consensus sequence     
      alignment_matrix <- as.matrix(alignment)
      consensus_seq <- ""
      # Iterate through each column in the alignment
      for (i in 1:ncol(alignment_matrix)) {
        # Extract the bases at the current column
        bases <- alignment_matrix[, i]
        
        # Check if the bases are valid (A, C, T, or G)
        valid_bases <- bases[bases %in% c("A", "C", "T", "G", "a", "c", "t", "g")]
        
        # If there are valid bases, take the first one (or apply any consensus rule you prefer)
        if (length(valid_bases) > 0) {
          consensus_seq <- paste0(consensus_seq, valid_bases[1])
        }
      }
      
      # consensus_seq <-
      #   msaConsensusSequence(
      #     alignment,
      #     # type = "upperlower",
      #     # ignoreGaps = TRUE
      #   ) %>%
      #   str_remove_all(.,
      #                  "\\?")
      
      # Print the consensus sequence
      consensus_sequences[[specimen]] <- paste(consensus_seq, collapse = "")
      
    } 
    
    # Identify specimens with only one ab1 file
    single_ab1_specimens <- Filter(function(x) (length(x$fwd) == 1 & length(x$rev) == 0) || (length(x$fwd) == 0 & length(x$rev) == 1), file_pairs)
    
    # Loop through each specimen with a single ab1 file
    for(specimen in names(single_ab1_specimens)) {
      # Assuming a single ab1 file could be either fwd or rev
      ab1_file <- ifelse(length(single_ab1_specimens[[specimen]]$fwd) == 1, 
                         single_ab1_specimens[[specimen]]$fwd, 
                         single_ab1_specimens[[specimen]]$rev)
      
      # Read the sequence from the ab1 file
      ab1_seq <- sangerseqR::read.abif(ab1_file[[1]])@data$PBAS.1
      
      # Print the consensus sequence
      consensus_sequences[[specimen]] <- paste(ab1_seq, collapse = "")
      
    }
    
    dna_seqs <- 
      DNAStringSet(unlist(consensus_sequences)) %>%
      msa() %>% 
      DNAStringSet()
    
    # Write to a FASTA file
    writeXStringSet(dna_seqs, 
                    filepath = out_file)
    
  }

