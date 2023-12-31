#### LIBRARIES ####

packages_used <-
  c(
    "tidyverse",
    "rentrez",
    "ape",
    "purrr",
    "furrr",
    "future",
    "haplotypes",
    "styler"
  )

packages_to_install <-
  packages_used[!packages_used %in% 
                  installed.packages()[, 1]]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install,
                   Ncpus = n_cpu - 1)
}

#below shows list of packages that worked installing and loading
lapply(packages_used,
       require,
       character.only = TRUE
)

# Function to check and install a package if not already installed
check_and_install <- 
  function(package) {
    if (!require(package, 
                 character.only = TRUE)) {
      install.packages(package)
      library(package, 
              character.only = TRUE)
    }
  }

if (!require("restez", character.only = TRUE)) {
  # Check and install 'remotes' for installing from GitHub
  check_and_install("remotes")
  
  # Install 'restez' from GitHub
  remotes::install_github("ropensci/restez", force = TRUE)
}

library(restez)


getGenBankSeqs <- 
  function(
    search_term, 
    out_genbank_file_path = "./ncbi.gb", 
    out_fasta_file_path = "./ncbi.fasta", 
    pattern_filter_out = "", 
    n_cpu = parallel::detectCores()
  ) {
    
    #### USER DEFINED VARIABLES ####
    
    # search_term <- '((txid79456[Organism]) AND ("COI" OR "COX1" OR "COXI" OR "CO1" OR "OXIDASE SUBUNIT 1" OR "mitochondrial genome")  NOT ("edna" OR "environmental" OR "chromosome")) OR ((txid603550[Organism:noexp]) AND HQ566727[Accession])'
    # out_genbank_file_path <- "../prj_rotablue_barcoding/output/ischnura_coi_ncbi.gb"
    # out_fasta_file_path <- "../prj_rotablue_barcoding/output/ischnura_coi_ncbi.fasta"
    # 
    # # regular expression to remove unwanted seqs in query
    # pattern_filter_out <- "ribosomal|RNA|16[Ss]|12[Ss]|cytochrome b|cytb|cyt b"
    # 
    # # number of cpu on your computer, for parallelization of code
    # n_cpu <- parallel::detectCores()
    
    
    #### SETUP QUERY GENBANK ####
    
    # Use the entrez_search function to perform the search and save the UIDs (unique identifiers) of the records:
    search_results <-
      rentrez::entrez_search(
        db="nucleotide",
        term=search_term,
        retmax=5000,  #Note: I set retmax to 5000 here to ensure all 4930 results are retrieved.
        use_history=TRUE
      )
    
    
    # Set up variables to fetch records in batches
    total_records <- 
      search_results$count
    batch_size <- 200 # Fetch 500 records at a time
    batches <-
      seq(
        from = 1,
        to = total_records,
        by = batch_size
      )
    
    # #### QUERY GENGANK IF RUNNING FOR 1ST TIME OR CHANGED `search_term` ####
    #
    # ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # #RUN THIS EACH TIME YOU CHANGE THE `search_term` VARIABLE
    # ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    #
    # # Create an empty vector to store all records
    all_records <- ""
    
    # This will fetch the records in batches of 200 at a time using the web history from your search, and write them all into a text file in GenBank format. If you still encounter issues, consider reducing the batch size even further.
    # only need to run if running for first time or changing the query.  results have been saved to a file, so CEB commented this out
    for(batch in batches){
      records <-
        rentrez::entrez_fetch(
          db="nucleotide",
          web_history=search_results$web_history,
          # WebEnv=search_results$WebEnv,
          # query_key=search_results$query_key,
          retstart=batch,
          retmax=batch_size,
          rettype="gb",
          retmode="text"
        )
      all_records <- paste(all_records, records, sep="")
    }
    
    # Save the records to a file
    write(all_records,
          out_genbank_file_path)
    
    rm(all_records)
    
    #### WRANGLE GENBANK RECORDS ####
    
    # REad in and Split "all_records" into individual records
    genbank_records <- 
      paste(
        readLines(out_genbank_file_path), 
        collapse = "\n"
      ) %>%
      str_split("\n//\n") %>%
      .[[1]]
    
    # define function to parse genbank format
    parseGenbankRecord <- 
      function(gb_record) {
        locus <- 
          as_tibble(
            as.list(
              gb_extract(
                gb_record, 
                "locus"
              )
            )
          )
        
        features <- 
          as_tibble(
            gb_extract(
              gb_record, 
              "features"
            )[[1]]
          )
        
        # definition <- 
        #   as_tibble(
        #     gb_extract(gb_record, 
        #                "definition")
        #   ) %>%
        #   dplyr::rename(definition = value)
        # 
        # 
        # keywords <-
        #   as_tibble(
        #     gb_extract(gb_record,
        #                "keywords")
        #   ) %>%
        #   dplyr::rename(keywords = value)
        
        definition <- 
          as_tibble(
            list(
              definition = 
                gb_extract(
                  gb_record, 
                  "definition"
                )
            )
          )
        
        keywords <- 
          as_tibble(
            list(
              keywords = 
                gb_extract(
                  gb_record, 
                  "keywords"
                )
            )
          ) 
        
        sequence <- 
          as_tibble(
            list(
              sequence = 
                gb_extract(
                  gb_record, 
                  "sequence"
                )
            )
          ) 
        
        bind_cols(
          locus,
          features,
          definition,
          sequence,
          keywords
        )
      }
    
    # # Set up a parallel plan
    # plan(multisession,
    #      workers = n_cpu -1)
    
    # Parse each record into one tibble, this takes a little while
    data_genbank <- 
      # future_map_dfr(
      map_df(
        genbank_records[1:total_records - 1], 
        parseGenbankRecord
      ) %>%
      # remove non target seqs
      # remove non target seqs
      dplyr::filter(
        if_any(
          everything(),
          ~ !str_detect(.,
                        pattern_filter_out) 
        ) 
      ) %>%
      # mild haplotype collapsing prior to alignment
      distinct(
        sequence,
        organism,
        # country,
        .keep_all = TRUE
      ) %>%
      # create names for fasta seqs
      dplyr::mutate(
        organism = str_remove(organism,
                              " [a-zA-Z0-9]*\\-.*$"),
        organism = str_remove(organism,
                              "\\'.*$"),
        organism = str_remove(organism,
                              "cf\\. *"),
        organism = str_remove(organism,
                              " [A-Z0-9_].*$"),
        organism = str_remove(organism,
                              "\\..*$"),
        country = str_remove(country,
                             " *<.*$"),
        country = str_replace(country,
                              "  *",
                              " "),
        location = case_when(
          !is.na(country) ~ str_sub(
            country,
            1,
            20),
          !is.na(lat_lon) ~ lat_lon,
          TRUE ~ "unk"),
        fasta_names = str_c(
          organism,
          accession,
          location,
          sep = " "
        )
      )
    
    #### CREATE A FASTA FILE FROM TIBBLE ####
    
    # # Open a file to write the FASTA sequences
    # write_connection <- 
    #   file(
    #     out_fasta_file_path, 
    #     open = "wt"
    #   )
    # 
    # # Iterate through the rows of the tibble
    # for (i in 1:nrow(data_genbank)) {
    #   # Get the accession and sequence
    #   accession <- data_genbank$fasta_names[i]
    #   sequence <- data_genbank$sequence[i]
    #   
    #   # Print the accession and sequence in FASTA format
    #   cat(paste0(">", accession, "\n", sequence, "\n"), file = write_connection)
    # }
    # 
    # # Close the file connection
    # close(write_connection)
    
    
    # Define a function to format as a FASTA entry
    format_fasta <- 
      function(fasta_names, sequence) {
        paste0(
          ">", 
          fasta_names, 
          "\n", 
          sequence, 
          "\n")
      }
    
    # Write the FASTA sequences to a file
    writeLines(
      data_genbank %>%
        select(fasta_names,
               sequence) %>%
        pmap_chr(function(...) format_fasta(fasta_names = ..1, 
                                            sequence = ..2)), 
      out_fasta_file_path
    )
    
  }