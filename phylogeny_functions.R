#!/usr/bin/env Rscript

#### PACKAGES ####

cran_packages <- c(
  "tidyverse", "BiocManager", "janitor", "phangorn", "styler", "ape",
  "tinytex", "igraph", "ips", "seqinr", "ggtext", "bios2mds", "rgl",
  "RCurl", "devtools", "ggrepel", "haplotypes", "tidytree", "ggimage"
)

special_packages <- list(
  list(pkg = "Biostrings",    source = "Bioc",   repo = "Biostrings",               install_options = list()),
  list(pkg = "seqLogo",       source = "Bioc",   repo = "seqLogo",                  install_options = list()),
  list(pkg = "microRNA",      source = "Bioc",   repo = "microRNA",                 install_options = list()),
  list(pkg = "ggtree",        source = "Bioc",   repo = "ggtree",                   install_options = list()),
  list(pkg = "msa",           source = "Bioc",   repo = "msa",                      install_options = list()),
  list(pkg = "ggimage",       source = "GitHub", repo = "GuangchuangYu/ggimage",    install_options = list()),
  list(pkg = "R4RNA",         source = "Bioc",   repo = "R4RNA",                    install_options = list()),
  list(pkg = "ggmsa",         source = "GitHub", repo = "YuLab-SMU/ggmsa",           install_options = list()),
  list(pkg = "treedataverse", source = "Bioc",   repo = "YuLab-SMU/treedataverse",   install_options = list()),
  list(pkg = "usedist",       source = "GitHub", repo = "kylebittinger/usedist",     install_options = list()),
  list(pkg = "rMSA",          source = "GitHub", repo = "mhahsler/rMSA",              install_options = list())
)

setRepositories(ind=1:2)

# A helper function that:
# 1. Attempts to load a package (quietly, without startup messages).
# 2. If the package is not available, it installs the package using the
#    specified source (CRAN, Bioc, or GitHub) and installation options.
# 3. Finally, it attempts to load the package again, printing a clear message
#    if the package still fails to load.
load_install_pkg <- 
  function(pkg, source = "CRAN", repo = NULL, install_options = list()) {
    if (!suppressPackageStartupMessages(require(pkg, character.only = TRUE))) {
      message("Package '", pkg, "' not found. Attempting installation from ", source, " ...")
      tryCatch({
        if (source == "CRAN") {
          install.packages(pkg, Ncpus = parallel::detectCores() - 1)
        } else if (source == "Bioc") {
          do.call(BiocManager::install, c(list(pkg), install_options))
        } else if (source == "GitHub") {
          if (is.null(repo)) {
            stop("For GitHub installation, 'repo' must be specified for package '", pkg, "'.")
          }
          do.call(devtools::install_github, c(list(repo), install_options))
        } else {
          stop("Unknown installation source '", source, "' for package '", pkg, "'.")
        }
      }, error = function(e) {
        message("Installation of package '", pkg, "' failed: ", e$message)
      })
      
      if (!suppressPackageStartupMessages(require(pkg, character.only = TRUE))) {
        message("Failed to load package '", pkg, "' even after installation. Please check the installation or try again later.")
      } else {
        message("Successfully loaded package '", pkg, "'.")
      }
    } else {
      message("Package '", pkg, "' loaded successfully.")
    }
  }

# Loop through and process each CRAN package.
for (pkg in cran_packages) {
  load_install_pkg(pkg = pkg, source = "CRAN")
}

# Loop through and process each special package.
for (pkg_info in special_packages) {
  load_install_pkg(pkg = pkg_info$pkg,
                   source = pkg_info$source,
                   repo = pkg_info$repo,
                   install_options = pkg_info$install_options)
}

rm(pkg_info, special_packages, bioc_pkgs, cran_packages, cran_pkgs, pkg)


#### Concatenate Fasta Files ####

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
    
    # Ensure unique labels by modifying the row names
    rownames(data_haplotypes_dnabin) <- make.unique(rownames(data_haplotypes_dnabin), sep = "_")
    
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
    model_ml = NULL, # new argument: if evolutionary model is provided, skip model test
    my_outgroup = "HQ909472_Kraemaria_bryani_isolate_C258_outgroup Kraemeria bryani isolate C258 voucher LACM:T-000093 cytochrome oxidase subunit I (CO1) gene, partial cds; mitochondrial",
    n_bootstraps = 100,
    threshold_bootstraps = 50,
    n_cpu = parallel::detectCores()
  ){
    
    #### ML Model Selection & Optimization ####
    if(is.null(model_ml)) {
      ## test models, select lowest AIC
      data_modeltest <-
        .data_fasta %>%
        phangorn::modelTest(
          model = "all",
          multicore = ifelse(n_cpu - 1 == 1, FALSE, TRUE),
          mc.cores = ifelse(n_cpu - 1 == 1, NULL, n_cpu - 1)
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
      
      print(paste("The best evolutionary model is:", modeltest_as.pml_bestfit))
      
    } else {
      
      data_modeltest <-
        .data_fasta %>%
        phangorn::modelTest(
          model = model_ml,
          multicore = ifelse(n_cpu - 1 == 1, FALSE, TRUE),
          mc.cores = ifelse(n_cpu - 1 == 1, NULL, n_cpu - 1)
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
      
      print(paste("The best version of the ", model_ml, 
                  " evolutionary model is:", modeltest_as.pml_bestfit))
      
    }
    
    #### Optimize Fit ####
    
    ## MAKE SURE TO CHANGE OPTIONS DEPENDING ON YOUR SELECTED MODEL.
    # ?optim.pml gives guidance on how to set optBf and optQ
    # optInv and optGamma should be set depending on whether your model includes +I and/or +G parameters
    
    # Extract the base model by removing +I and +G
    base_model <- 
      gsub("\\+I|\\+G.*", "", modeltest_as.pml_bestfit) %>%
      trimws()
    
    # Automatically set optBf and optQ based on the base model
    if(base_model == "JC") {
      auto_optBf <- FALSE
      auto_optQ <- FALSE
    } else if(base_model == "K80") {
      auto_optBf <- FALSE
      auto_optQ <- TRUE
    } else if(base_model == "F81") {
      auto_optBf <- TRUE
      auto_optQ <- FALSE
    } else if(base_model == "SYM") {
      auto_optBf <- FALSE
      auto_optQ <- TRUE
    } else {
      # default to optimizing both if the model is not one of the above
      auto_optBf <- TRUE
      auto_optQ <- TRUE
    }
    
    # Determine if the best model includes +I or +G parameters
    auto_optInv <- grepl("\\+I", modeltest_as.pml_bestfit)
    auto_optGamma <- grepl("\\+G", modeltest_as.pml_bestfit)
    
    # print(
    cat("Optimized evolutionary model parameters:\n",
        "optBf =", auto_optBf, "\n",
        "optQ =", auto_optQ, "\n",
        "optInv =", auto_optInv, "\n",
        "optGamma =", auto_optGamma, "\n")
    # )
    
    evolModelFit_opt <-
      data_modeltest_fit %>%
      optim.pml(
        optBf = auto_optBf,
        optQ = auto_optQ,
        optInv = auto_optInv,
        optGamma = auto_optGamma,
        rearrangement = "NNI",
        control = pml.control(trace = 0)
      )
    evolModelFit_opt
    
    #### and Bootstrap ####
    # bootstrap model
    # Skip NNi if less than 5 tips on tree since that can lead to issues if there are branches with 0 length that are collapsed down to 3 or fewer branches.
    if(length(evolModelFit_opt$data) < 5){
      trees_evolModelFit_opt_bs <-
        bootstrap.pml(
          evolModelFit_opt,
          bs = n_bootstraps,
          optNni =  FALSE,
          multicore = ifelse(n_cpu - 1 == 1, FALSE, TRUE),
          mc.cores = ifelse(n_cpu - 1 == 1, NULL, n_cpu - 1),
        )
    } else {
      trees_evolModelFit_opt_bs <-
        bootstrap.pml(
          evolModelFit_opt,
          bs = n_bootstraps,
          optNni =  TRUE,
          multicore = ifelse(n_cpu - 1 == 1, FALSE, TRUE),
          mc.cores = ifelse(n_cpu - 1 == 1, NULL, n_cpu - 1),
        )
    }
    
    ## plotBS functions
    # type = the type of tree to plot, one of "phylogram", "cladogram", "fan", "unrooted", "radial" or "none". If type is "none" the tree is returned with the bootstrap values assigned to the node labels.
    # method = either "FBP" the classical bootstrap (default) or "TBE" (transfer bootstrap)
    # digits = nteger indicating the number of decimal places.
    # p	= only plot support values higher than this percentage number (default is 0).
    
    tree_evolModelFit_opt_bs <-
      plotBS(
        evolModelFit_opt$tree,
        trees_evolModelFit_opt_bs,
        p = threshold_bootstraps/100,
        digits = 2,
        type = "phylogram",
        method = "FBP"
      )
    
    # if max bs val is less than or equal to 1, then convert to scale of 0-100
    conversion_factor <- 
      if (
        max(
          as.numeric(
            tree_evolModelFit_opt_bs$node.label
          ), 
          na.rm = TRUE) <= 1
      ) 100 else 1
    
    # Convert bootstrap support values to values between 0 and 100
    bootstraps_sig <- 
      (as.numeric(tree_evolModelFit_opt_bs$node.label) * conversion_factor) %>% 
      { ifelse(. < threshold_bootstraps, "", .) } %>%
      { ifelse(is.na(.), "", .) }
    
    #### SPECIFY OUTGROUP & ROOT ####
    # this takes the accession number and makes sure that the name of the outgroup matches naming used in tip labels
    
    if(!is.null(my_outgroup)){
      the_outgroup <-
        tree_evolModelFit_opt_bs$tip.label %>%
        as_tibble() %>%
        filter(str_detect(value,
                          gsub("[_\\s].*",
                               "",
                               my_outgroup))) %>%
        pull()
      
      print(
        paste(
          "The outgroup is: ", the_outgroup
        )
      )
    } else {
      the_outgroup <- character()
    }
    
    
    if(length(the_outgroup) > 0) {
      tree_evolModelFit_opt_bs_outgroup <-
        unroot(tree_evolModelFit_opt_bs) %>%
        root(outgroup = the_outgroup)
      
      tree_evolModelFit_opt_bs_outgroup
      
      cat("Outgroup found in the tree. Returning the rerooted tree.\n")
      return(
        list(
          model_ml = data_modeltest_fit,
          best_model = modeltest_as.pml_bestfit,
          model_optim.pml = evolModelFit_opt,
          tree = tree_evolModelFit_opt_bs_outgroup, 
          bootstraps_sig = bootstraps_sig
        )
      )
    } else {
      # If the_outgroup is empty, return the original tree without re-rooting
      cat("Outgroup not found in the tree. Returning the original tree.\n")
      return(
        list(
          model_ml = data_modeltest_fit,
          best_model = modeltest_as.pml_bestfit,
          model_optim.pml = evolModelFit_opt,
          tree = tree_evolModelFit_opt_bs, 
          bootstraps_sig = bootstraps_sig
        )
      )    
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


