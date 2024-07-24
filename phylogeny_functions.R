#!/usr/env Rscript

#### PACKAGES ####
if (!require("Biostrings")) {
  BiocManager::install("Biostrings")
}

if (!require("seqLogo")) {
  BiocManager::install("seqLogo")
}

if (!require("microRNA")) {
  BiocManager::install("microRNA")
}

if (!require("ggtree")) {
  BiocManager::install("ggtree")
}

if (!require("msa")) {
  BiocManager::install("msa")
}

if (!require("ggimage")) {
  BiocManager::install(
    "ggimage",
    type = "binary"
  )
}

install.packages("R4RNA", repos = "https://bioconductor.org/packages/3.16/bioc")

if (!require("YuLab-SMU/ggmsa")) {
  devtools::install_github(
    "YuLab-SMU/ggmsa",
    type = "binary"
  )
}

if (!require("treedataverse")) {
  BiocManager::install(
    "YuLab-SMU/treedataverse",
    type = "binary",
    force = TRUE
  )
}

if (!require("usedist")) {
  devtools::install_github("kylebittinger/usedist")
}

if (!require("rMSA")) {
  devtools::install_github("mhahsler/rMSA")
}

packages_used <-
  c(
    "tidyverse",
    "BiocManager",
    "janitor",
    "phangorn",
    "styler",
    "ape",
    "msa",
    "tinytex",
    "igraph",
    "ips",
    "seqinr",
    "ggtext",
    "bios2mds",
    "rgl",
    "RCurl",
    "devtools",
    "ggrepel",
    "haplotypes",
    "tidytree",
    "Biostrings",
    "seqLogo",
    "microRNA",
    "ggtree",
    "ggimage",
    "ggmsa",
    "treedataverse",
    "usedist",
    "rMSA"
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

# ## load required library
# library(Biostrings)
# library(microRNA)
# library(ggtree)
# library(ggimage)
# library(ggmsa)
# library(ggrepel)


#### Concatonate Fasta Files ####

# my_seqs_fasta_file = "../data/sequences/eviota_CO1_aligned_trimmed.fasta"
# ncbi_fasta_file = "../output/eviota_coi_ncbi.fasta"
# concat_fasta_file = "../output/eviota_CO1_aligned_trimmed_ncbi.fasta"
# edited_concat_fasta_file = "../output/eviota_CO1_aligned_trimmed_ncbi_aligned_trimmed.fasta"
# no_ambigu_edited_concat_fasta_file = "../output/eviota_CO1_aligned_trimmed_ncbi_aligned_trimmed_no_ambiguities.fasta"
# uniq_concat_fasta_file = "../output/eviota_CO1_aligned_trimmed_ncbi_aligned_trimmed_haps.fasta"
# outgroup_acession_number = "HQ909472"
# output_tree_name = "../output/eviota_CO1_aligned_trimmed_ncbi_aligned_trimmed_haps.nwk"


concatFastas <-
  function(inFilePaths,
           outFilePath){
    # Read and concatenate the content of all FASTA files
    combined_fasta <- 
      inFilePaths %>%
      map_chr(readr::read_file) %>%
      paste(collapse = "\n") %>%
      str_remove_all(.,"\r")
    
    # Write the combined content to a new file
    writeLines(
      combined_fasta, 
      outFilePath
    )
  }

#### UNIQUE HAPLOTYPES ####

# open `concat_fasta_file` with seaview, align, and manually trim
# Manually edit the names of your specimen's seqs
# save this file in seaview with the name assigned to the `edited_concat_fasta_file` variable


uniqueSeqsFastaFile <-
  function(
    inFilePath = "../data/myfasta.fasta",
    outFilePath = "../output/myfasta_haplotypes.fasta"
  ) {
    data_haplotypes_dnabin <-
      haplotypes::read.fas(file = inFilePath) %>%
      # generates haplotype files#
      # selection of indel treatment is rather important here.  Most predictable is to set them as a 5th character state
      haplotype(indels = "5th") %>%
      # some file format contortions...
      as.dna() %>%
      haplotypes::as.DNAbin()
    
    write.dna(data_haplotypes_dnabin,
              file = outFilePath,
              format = "fasta"
    )
  }


#### FILTER GENBANK FILES AFTER ALIGINGING and TRIMMING ####

# open the `edited_concat_fasta_file` in notepad++ or bbedit and remove sequences with the invalid characters from the concat function and 
# save as `no_ambigu_edited_concat_fasta_file`

# filterSimilarSequences(
#   fasta_keep = "../output/eviota_CO1_aligned_trimmed_ncbi_aligned_trimmed_no_ambiguities_phires_only.fasta",
#   fasta_ncbi = "../output/eviota_CO1_aligned_trimmed_ncbi_aligned_trimmed_no_ambiguities_ncbi_only.fasta",
#   outFilePath = "../output/eviota_CO1_aligned_trimmed_ncbi_aligned_trimmed_no_ambiguities_filtered.fasta",
#   threshold = 0.8
# )

filterSimilarSequences <- 
  function(
    fasta_keep, 
    fasta_ncbi, 
    outFilePath,
    threshold
  ) {
    # Load sequences
    my_seqs <- read.fasta(fasta_keep)
    ncbi_seqs <- read.fasta(fasta_ncbi)
    
    # Initialize a list to hold NCBI sequences to keep
    ncbi_to_keep <- list()
    
    # Loop through each NCBI sequence
    for (ncbi_name in names(ncbi_seqs)) {
      ncbi_seq <- ncbi_seqs[[ncbi_name]]
      
      # Check similarity with each of 'my_seqs'
      for (my_name in names(my_seqs)) {
        my_seq <- my_seqs[[my_name]]
        
        # Calculate similarity
        similarity <- sum(my_seq == ncbi_seq) / length(my_seq)
        
        # If similarity exceeds threshold, keep this NCBI sequence
        if (similarity >= threshold) {
          ncbi_to_keep[[ncbi_name]] <- ncbi_seq
          break
        }
      }
    }
    
    # Combine 'my_seqs' with filtered 'ncbi_to_keep'
    combined_seqs <- c(my_seqs, ncbi_to_keep)
    
    # Write to a new FASTA file
    write.fasta(sequences = combined_seqs, names = names(combined_seqs), file.out = outFilePath)
  }


#### ALIGN FASTA ####
alignFastaFile <-
  function(
    inFilePath = "../data/myfasta.fasta",
    outFilePath = "../data/myfasta_aligned.fasta"
  ){
    
    # read in data
    alignment <-
      Biostrings::readDNAStringSet(inFilePath) %>%
      RNA2DNA() %>%
      # multi-sequence alignment.  Muscle is generally the fastest option here.  Other options are ClustalW and ClustalOmega
      msa(
        method = "Muscle",
        type = "dna",
        verbose = TRUE
      ) 
    
    writeXStringSet(DNAStringSet(alignment), 
                    filepath = outFilePath)
    
    alignment %>%
      DNAMultipleAlignment() %>%   # nt format
      msaConvert("phangorn::phyDat") # pd format
  }

#### WRANGLE DATA FROM FASTA INTO TREE ####

fasta2tree <-
  function(
    .data_fasta,
    my_outgroup = "HQ909472_Kraemaria_bryani_isolate_C258_outgroup Kraemeria bryani isolate C258 voucher LACM:T-000093 cytochrome oxidase subunit I (CO1) gene, partial cds; mitochondrial",
    n_bootstraps = 100,
    threshold_bootstraps = 50,
    n_cpu = parallel::detectCores()
  ){
    
    #### ML Model Selection ####
    ## test models, select lowest AIC
    data_modeltest <-
      .data_fasta %>%
      modelTest(
        model = "all",
        multicore = TRUE,
        mc.cores = n_cpu - 1
      )
    
    #### Initial ModelFit ####
    ## initial fit, use best fit model from prev chunk
    data_modeltest_fit <-
      as.pml(data_modeltest)
    # data_modeltest_fit
    
    # automatically read output from previous line and return model in next line
    modeltest_as.pml_bestfit <- 
      data_modeltest_fit$call$tree %>% 
      as.character() %>% 
      gsub("tree_", "", .)
    
    print(modeltest_as.pml_bestfit)
    
    #### Optimize Fit ####
    
    ## MAKE SURE TO CHANGE OPTIONS DEPENDING ON YOUR SELECTED MODEL.
    # ?optim.pml gives guidance on how to set optBf and optQ
    # optInv and optGamma should be set depending on whether your model includes +I and/or +G parameters
    
    evolModelFit_opt <-
      data_modeltest_fit %>%
      optim.pml(
        optBf = TRUE,
        optQ = TRUE,
        optInv = TRUE,
        optGamma = FALSE,
        rearrangement = "NNI",
        control = pml.control(trace = 0)
      )
    evolModelFit_opt
    
    # Dave, is this doing what we want?
    # automatically read output from previous line and return model in next line
    modeltest_optim.pml_bestfit <- 
      evolModelFit_opt %>% 
      capture.output() %>%
      .[1] %>%
      gsub(pattern = "model:\\s",
           replacement = "") %>%
      gsub(pattern = "\\s",
           replacement = "")
    
    modeltest_optim.pml_bestfit
    
    #### and Bootstrap ####
    # bootstrap model
    trees_evolModelFit_opt_bs <-
      bootstrap.pml(
        evolModelFit_opt,
        bs = n_bootstraps,
        optNni = TRUE,
        multicore = TRUE,
        mc.cores = n_cpu - 1
      )
    
    ## plotBS functions
    # type = the type of tree to plot, one of "phylogram", "cladogram", "fan", "unrooted", "radial" or "none". If type is "none" the tree is returned with the bootstrap values assigned to the node labels.
    # method = either "FBP" the classical bootstrap (default) or "TBE" (transfer bootstrap)
    # digits = nteger indicating the number of decimal places.
    # p	= only plot support values higher than this percentage number (default is 0).
    
    tree_evolModelFit_opt_bs <-
      plotBS(
        evolModelFit_opt$tree,
        trees_evolModelFit_opt_bs,
        p = threshold_bootstraps,
        digits = 0,
        type = "phylogram",
        method = "FBP"
      )
    
    #### SPECIFY OUTGROUP & ROOT ####
    # this takes the accession number and makes sure that the name of the outgroup matches naming used in tip labels
    
    
    the_outgroup <-
      tree_evolModelFit_opt_bs$tip.label %>%
      as_tibble() %>%
      filter(str_detect(value,
                        gsub("[_\\s].*",
                             "",
                             my_outgroup))) %>%
      pull()
    
    if(length(the_outgroup) > 0) {
      tree_evolModelFit_opt_bs_outgroup <-
        unroot(tree_evolModelFit_opt_bs) %>%
        root(outgroup = the_outgroup)
      
      tree_evolModelFit_opt_bs_outgroup
      
      # plot(
      #   tree_evolModelFit_opt_bs_outgroup,
      #   main = modeltest_optim.pml_bestfit
      # )
      cat("Outgroup found in the tree. Returning the rerooted tree.\n")
      return(tree_evolModelFit_opt_bs_outgroup)
    } else {
      # If the_outgroup is empty, return the original tree without re-rooting
      cat("Outgroup not found in the tree. Returning the original tree.\n")
      return(tree_evolModelFit_opt_bs)
    }
  }

#### EDIT BRANCH TIP LABELS ####

editTipLabels <-
  function(tree){
    # change tip names, easier to change it inside the tree data structure than outside
    tree$tip.label <-
      tree$tip.label %>%
      # replace _ with spaces
      str_replace_all(
        "_",
        " "
      ) %>%
      # remove from co1 to end of tip label
      str_remove("[\\s_]isolate.*$") %>%
      str_remove("[\\s_]cytochrome.*$") %>%
      str_remove("[\\s_]voucher.*$") %>%
      str_remove("[\\s_]OCF.*$") %>%
      str_remove("[\\s_]LT\\-20.*$")  %>%
      str_remove("[\\s_]KSA.*$") 
    
    return(tree)
  }

#### generate basic GGTREE ####

plotGgTree <- 
  function(
    tree, 
    threshold_bootstraps = 50, 
    tip_color = "red", 
    tip_pattern = "^G",
    show_node_id = FALSE
  ) {
    
    threshold_bootstraps <<- threshold_bootstraps
    
    p <- 
      ggtree(tree) +
      geom_treescale() +
      theme_tree2() +
      xlim(NA, 1.4)
    
    # Conditionally add bootstrap labels
    if (threshold_bootstraps > 0) {
      p <- 
        p +
        geom_nodelab(
          node = "internal",
          hjust = 1,
          vjust = 0.5,
          nudge_x = -.008,
          aes(
            subset = label > threshold_bootstraps,
            label = label
          )
        ) 
    }
    
    if (show_node_id == TRUE) {
      p <- 
        p +
        geom_text(
          aes(label = node), 
          hjust=-.3
        )
    }
    
    # Adding the ability to color tip labels based on pattern
    if (show_node_id == FALSE) {
      p <- 
        p + 
        geom_tiplab(
          aes(
            label = label,
            color = ifelse(grepl(tip_pattern, label), tip_color, "black")
          )
        ) +
        scale_color_identity()
    }
    
    return(p)
  }

#### SAVE TREE ####

saveNewickTree <- 
  function(tree, file_path) {
    write.tree(tree, file = file_path)
  }


