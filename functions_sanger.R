#### PACKAGES ####

bioc_pkgs <- c(
  "sangerseqR",
  "Biostrings",
  "DECIPHER",
  "msa"
)

# CRAN packages
cran_pkgs <- c(
  "tidyverse",
  "ape",
  "janitor",
  "pegas"
)

#github_fallbacks <- c("somePackage" = "owner/somePackage")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)

for (pkg in bioc_pkgs) {
  # If not installed, install via BiocManager
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
  # Load the package after installation attempt
  library(pkg, character.only = TRUE)
}

for (pkg in cran_pkgs) {
  # First, try installing/loading from CRAN using p_load()
  if (!requireNamespace(pkg, quietly = TRUE)) {
    pacman::p_load(char = pkg, install = TRUE, update = FALSE, character.only = TRUE)
  } else {
    # Already installed, just load
    library(pkg, character.only = TRUE)
  }
  
  # If STILL not installed, fallback to GitHub if a repo is known
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% names(github_fallbacks)) {
      message(paste("Attempting GitHub install for", pkg))
      pacman::p_load_gh(github_fallbacks[[pkg]], character.only = TRUE, update = FALSE)
    } else {
      warning(paste("Package", pkg, "not found on CRAN, and no GitHub repo specified."))
    }
  }
}

rm(pkg_info, special_packages, bioc_pkgs, cran_packages, cran_pkgs, pkg)

# if(!require("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }
# 
# if (!require("sangerseqR")) {
#   BiocManager::install("sangerseqR")
# }
# 
# if (!require("DECIPHER")) {
#   BiocManager::install("DECIPHER")
# }
# 
# packages_used <-
#   c(
#     "tidyverse",
#     "ape",
#     "janitor",
#     "sangerseqR",
#     "Biostrings",
#     "DECIPHER",
#     "msa",
#     "pegas"
#   )
# 
# packages_to_install <-
#   packages_used[!packages_used %in% installed.packages()[, 1]]
# 
# if (length(packages_to_install) > 0) {
#   install.packages(packages_to_install,
#                    Ncpus = parallel::detectCores() - 1
#   )
# }
# 
# lapply(packages_used,
#        require,
#        character.only = TRUE
# )


#### FILTER SEQUENCES FROM RENAMED BLAST MSA OUTPUT ####

filterFastaByAccession <- 
  function(
    fasta_file, 
    ingroup_file, 
    outgroup_file, 
    output_file, 
    use_regex = FALSE
  ) {
    # # Check that the ape package is available
    # if (!requireNamespace("ape", quietly = TRUE)) {
    #   stop("Package 'ape' is required. Please install it with install.packages('ape').")
    # }
    
    # Remove trailing . followed by one or more digits
    ingroups <- gsub("\\.[0-9]+$", "", trimws(readLines(ingroup_file)))
    outgroups <- gsub("\\.[0-9]+$", "", trimws(readLines(outgroup_file)))
    
    
    # # Read the accession numbers from the ingroup and outgroup files
    # ingroups <- readLines(ingroup_file)
    # outgroups <- readLines(outgroup_file)
    # 
    # # Remove any empty lines and trim whitespace
    # ingroups <- trimws(ingroups[ingroups != ""])
    # outgroups <- trimws(outgroups[outgroups != ""])
    
    # Combine into a unique set of accession numbers to keep
    keep_accessions <- unique(c(ingroups, outgroups))
    
    # Read the FASTA file (assumes DNA sequences)
    seqs <- ape::read.dna(fasta_file, format = "fasta")
    
    # Filter sequences: check if any accession number appears in the sequence header
    keep <- 
      sapply(rownames(seqs), function(header) {
        any(sapply(keep_accessions, function(acc) {
          if (use_regex) {
            grepl(acc, header)
          } else {
            grepl(acc, header, fixed = TRUE)
          }
        }))
      })
    
    filtered_seqs <- seqs[keep, ]
    
    if (length(filtered_seqs) == 0) {
      warning("No sequences matched the accession numbers provided.")
    } else {
      # Write the filtered sequences to the output FASTA file
      ape::write.dna(filtered_seqs, file = output_file, format = "fasta", nbcol = -1, colsep = "")
    }
    
    return(filtered_seqs)
  }




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
processCuratedSANGER <-
  function(
    in_files = "../output/sanger_curated_ab1_ischnura_luta/^.*\\.(ab1|scf)$",
    out_file = "../output/sanger_curated_ab1_ischnura_luta/consensus_sequences.fasta",
    fwd_primer_name = "LCO1490",
    rev_primer_name = "HCO2198"
  ){
    #### READ IN AB1 ####
    
    # Extract the directory and file pattern using stringr functions
    if (str_detect(in_files, "/")) {
      indir <- str_replace(in_files, "/[^/]*$", "")   # Remove everything after the last '/'
      in_file_pattern <- str_extract(in_files, "[^/]+$")  # Extract everything after the last '/'
    } else {
      indir <- "."
      in_file_pattern <- in_files
    }
    
    # read in files
    ab1_files <- 
      list.files(indir, 
                 pattern = in_file_pattern, 
                 full.names = TRUE)
    
    
    
    # # Initialize an empty list to store AB1 data
    # ab1_data_list <- list()
    # 
    # # Loop through each AB1 file and read its data
    # for(file_path in ab1_files) {
    #   file_name <- basename(file_path) # Extracts the file name from the path
    #   # ab1_data_list[[file_name]] <- sangerseqR::read.abif(file_path)
    #   ab1_data_list[[file_name]] <- sangerseqR::readsangerseq(file_path)
    # }
    
    #### MAKE CONSENSUS SEQUENCES ####
    
    # Correctly identify and organize pairs
    file_pairs <- list()
    
    for (file in ab1_files) {
      
      #need to fix this to be more compatible with errors in naming or need to fix names of files.
      # specimen_id <- sub("(.+)_([A-Z0-9]+[0-9A-Z]+)_(.+)_(.+)\\.ab1$", "\\1", basename(file))
      # primer_type <- sub("(.+)_([A-Z0-9]+[0-9A-Z]+)_(.+)_(.+)\\.ab1$", "\\2", basename(file))
      
      specimen_id <- sub("(.*)[_\\-](.*)[_\\-](\\d{4}-\\d{2}-\\d{2})[_\\-](.+)\\.(ab1|scf)$", "\\1", basename(file)) # specimen id
      primer_type <- sub("(.*)[_\\-](.*)[_\\-](\\d{4}-\\d{2}-\\d{2})[_\\-](.+)\\.(ab1|scf)$", "\\2", basename(file)) # primer
      # sub("(.*)[_\\-](.*)[_\\-](\\d{4}-\\d{2}-\\d{2})[_\\-](.+)\\.(ab1|scf)$", "\\3", basename(file)) # date
      # sub("(.*)[_\\-](.*)[_\\-](\\d{4}-\\d{2}-\\d{2})[_\\-](.+)\\.(ab1|scf)$", "\\4", basename(file)) # well etc
      # sub("(.*)[_\\-](.*)[_\\-](\\d{4}-\\d{2}-\\d{2})[_\\-](.+)\\.(ab1|scf)$", "\\5", basename(file)) # file ext
      
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
    file_pairs_filtered <- 
      Filter(
        function(x) length(x$fwd) == 1 && 
          length(x$rev) == 1, 
        file_pairs
      )
    
    consensus_sequences <- list()
    
    for (specimen in names(file_pairs_filtered)) {
      fwd_file <- file_pairs_filtered[[specimen]]$fwd
      rev_file <- file_pairs_filtered[[specimen]]$rev
      
      # Read sequences
      fwd_seq <- 
        sangerseqR::readsangerseq(fwd_file[[1]]) %>%
        primarySeq() %>%
        as.character()
      
      rev_seq <- 
        # sangerseqR::read.abif(rev_file[[1]])@data$PBAS.1 %>%
        sangerseqR::readsangerseq(rev_file[[1]]) %>%
        primarySeq() %>%
        DNAString() %>%
        # reverseComplement() %>%
        as.character()
      
      rev_seq_revcomp <- 
        # sangerseqR::read.abif(rev_file[[1]])@data$PBAS.1 %>%
        sangerseqR::readsangerseq(rev_file[[1]]) %>%
        primarySeq() %>%
        DNAString() %>%
        reverseComplement() %>%
        as.character()
      
      # Perform pairwise alignments using a local alignment (or adjust type as needed)
      score_rev_seq <- pairwiseAlignment(fwd_seq, rev_seq, type = "local", scoreOnly = TRUE)
      score_rev_sed_revcomp <- pairwiseAlignment(fwd_seq, rev_seq_revcomp, type = "local", scoreOnly = TRUE)
      
      # Choose the orientation that gives the higher alignment score
      if(score_rev_sed_revcomp >= score_rev_seq) {
        rev_seq <- as.character(rev_seq_revcomp)
      } 
      
      # First, combine them into a DNAStringSet
      sequences <- 
        DNAStringSet(
          c(
            fwd_seq, 
            rev_seq
          )
        )
      
      # Perform the alignment
      alignment <- 
        msa(sequences,
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
      ab1_seq <- 
        # sangerseqR::read.abif(ab1_file[[1]])@data$PBAS.1
        sangerseqR::readsangerseq(ab1_file[[1]]) %>%
        primarySeq() %>%
        as.character()
      
      # Print the consensus sequence
      consensus_sequences[[specimen]] <- paste(ab1_seq, collapse = "")
      
    }
    
    #dna_seqs <- 
    #  DNAStringSet(unlist(consensus_sequences)) %>%
    #  msa() %>% 
    #  DNAStringSet()
    
    # Create a DNAStringSet from the consensus sequences
    dna_seqs <- DNAStringSet(unlist(consensus_sequences))
    
    # Only run msa() if there is more than one sequence
    if (length(dna_seqs) > 1) {
      dna_seqs <- msa(dna_seqs) %>% DNAStringSet()
    } else {
      dna_seqs <- dna_seqs %>% DNAStringSet()
    }
    
    
    # Write to a FASTA file
    writeXStringSet(dna_seqs, 
                    filepath = out_file)
    
  }

processCuratedAB1 <-
  function(
    ab1_dir = "../output/sanger_curated_ab1_ischnura_luta/",
    out_file = "../output/sanger_curated_ab1_ischnura_luta/consensus_sequences.fasta",
    fwd_primer_name = "LCO1490",
    rev_primer_name = "HCO2198"
  ){
    processCuratedSANGER(
      in_files = 
        str_c(
          ab1_dir,
          "/.*\\.(ab1|scf|AB1|SCF)",
          sep=""
        ),
      out_file = out_file,
      fwd_primer_name = fwd_primer_name,
      rev_primer_name = rev_primer_name
    )
  }

processCuratedSCF <-
  function(
    scf_dir = "../output/sanger_curated_ab1_ischnura_luta/",
    out_file = "../output/sanger_curated_ab1_ischnura_luta/consensus_sequences.fasta",
    fwd_primer_name = "LCO1490",
    rev_primer_name = "HCO2198"
  ){
    processCuratedAB1(
      in_files = 
        str_c(
          scf_dir,
          "/.*\\.(ab1|scf|AB1|SCF)",
          sep=""
        ),
      out_file = out_file,
      fwd_primer_name = fwd_primer_name,
      rev_primer_name = rev_primer_name
    )
  }
