#### INITIALIZE ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### PACKAGES ####

cran_packages <- c(
  "tidyverse", 
  "BiocManager", 
  "janitor", 
  "devtools"
  # "phangorn", 
  # "styler", 
  # "ape",
  # "tinytex", 
  # "igraph", 
  # "ips", 
  # "seqinr", 
  # "ggtext", 
  # "bios2mds", 
  # "rgl",
  # "RCurl", 
  # "ggrepel", 
  # "haplotypes", 
  # "tidytree", 
  # "ggimage"
)

special_packages <- list(
  # list(pkg = "Biostrings",    source = "Bioc",   repo = "Biostrings",               install_options = list()),
  # list(pkg = "seqLogo",       source = "Bioc",   repo = "seqLogo",                  install_options = list()),
  # list(pkg = "microRNA",      source = "Bioc",   repo = "microRNA",                 install_options = list()),
  # list(pkg = "ggtree",        source = "Bioc",   repo = "ggtree",                   install_options = list()),
  # list(pkg = "msa",           source = "Bioc",   repo = "msa",                      install_options = list()),
  # list(pkg = "ggimage",       source = "GitHub", repo = "GuangchuangYu/ggimage",    install_options = list()),
  # list(pkg = "R4RNA",         source = "Bioc",   repo = "R4RNA",                    install_options = list()),
  # list(pkg = "ggmsa",         source = "GitHub", repo = "YuLab-SMU/ggmsa",           install_options = list()),
  # list(pkg = "treedataverse", source = "Bioc",   repo = "YuLab-SMU/treedataverse",   install_options = list()),
  # list(pkg = "usedist",       source = "GitHub", repo = "kylebittinger/usedist",     install_options = list()),
  # list(pkg = "rMSA",          source = "GitHub", repo = "mhahsler/rMSA",              install_options = list())
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
  load_install_pkg(
    pkg = pkg_info$pkg,
    source = pkg_info$source,
    repo = pkg_info$repo,
    install_options = pkg_info$install_options
  )
}

rm(special_packages, cran_packages, pkg, pkg_info)

source("df2fasta.R")
source("fasta2df.R")

#### WRANGLE SANGER SEQS 16S, COI, 28S, H3####

# read in sanger data, only keep isolates with all loci sequenced
loci_all <- c("H3", "COI", "16S", "28S")

data_sanger <- 
  list.files(
    path       = "../output/multilocus_phylogeny",
    pattern    = "ischnura_luta_.*consensus.*\\.fasta$",
    recursive  = TRUE,
    full.names = TRUE
  ) %>%
  set_names(.) %>%              # so .id is just the file name
  map_dfr(
    ~ fasta2df(
      .x
    ),
    .id = "file_path"                        # new column “source” ← basename(.x)
  ) %>%
  clean_names() %>%
  select(-extraction_id) %>%
  group_by(isolate) %>%                  # group by the column of interest
  filter(all(loci_all %in% locus)) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  arrange(organism, isolate, locus)


loci_all %>%
  map(
    ~ data_sanger %>%
      filter(locus == .x) %>%
      df2fasta(
        ., 
        fasta_seq_names = c("organism", "isolate"), 
        fasta_file = str_c(
          "../output/multilocus_phylogeny/16s_coi_28s_h3/ischnura_luta_",
          str_to_lower(.x),
          "_groupisolates.fasta"
        )
      )
  )


#### WRANGLE SANGER SEQS 16S, 28S, H3####

# read in sanger data, only keep isolates with all loci sequenced
loci_3 <- c("H3", "16S", "28S")

data_sanger_3 <- 
  list.files(
    path       = "../output/multilocus_phylogeny",
    pattern    = "ischnura_luta_.*consensus.*\\.fasta$",
    recursive  = TRUE,
    full.names = TRUE
  ) %>%
  set_names(.) %>%              # so .id is just the file name
  map_dfr(
    ~ fasta2df(
      .x
    ),
    .id = "file_path"                        # new column “source” ← basename(.x)
  ) %>%
  clean_names() %>%
  select(-extraction_id) %>%
  filter(locus %in% loci_3) %>%
  group_by(isolate) %>%                  # group by the column of interest
  filter(all(loci_3 %in% locus)) %>%
  filter(n() > 2) %>%
  ungroup() %>%
  arrange(organism, isolate, locus)

loci_3 %>%
  map(
    ~ data_sanger_3 %>%
      filter(locus == .x) %>%
      df2fasta(
        ., 
        fasta_seq_names = c("organism", "isolate"), 
        fasta_file = str_c(
          "../output/multilocus_phylogeny/16s_28s_h3/ischnura_luta_",
          str_to_lower(.x),
          "_groupisolates.fasta"
        )
      )
  )

#### WRANGLE SANGER SEQS 28, H3####

# read in sanger data, only keep isolates with all loci sequenced
loci_nuc <- c("H3", "28S")

data_sanger_nuc <- 
  list.files(
    path       = "../output/multilocus_phylogeny",
    pattern    = "ischnura_luta_.*consensus.*\\.fasta$",
    recursive  = TRUE,
    full.names = TRUE
  ) %>%
  set_names(.) %>%              # so .id is just the file name
  map_dfr(
    ~ fasta2df(
      .x
    ),
    .id = "file_path"                        # new column “source” ← basename(.x)
  ) %>%
  clean_names() %>%
  select(-extraction_id) %>%
  filter(locus %in% loci_nuc) %>%
  group_by(isolate) %>%                  # group by the column of interest
  filter(all(loci_nuc %in% locus)) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  arrange(organism, isolate, locus)


loci_nuc %>%
  map(
    ~ data_sanger_nuc %>%
      filter(locus == .x) %>%
      df2fasta(
        ., 
        fasta_seq_names = c("organism", "isolate"), 
        fasta_file = str_c(
          "../output/multilocus_phylogeny/28s_h3/ischnura_luta_",
          str_to_lower(.x),
          "_groupisolates.fasta"
        )
      )
  )

#### WRANGLE SANGER SEQS 16S, COI, ####

# read in sanger data, only keep isolates with all loci sequenced
loci_mito <- c("COI", "16S")

data_sanger_mito <- 
  list.files(
    path       = "../output/multilocus_phylogeny/",
    pattern    = "ischnura_luta_.*consensus.*\\.fasta$",
    recursive  = TRUE,
    full.names = TRUE
  ) %>%
  set_names(.) %>%              # so .id is just the file name
  map_dfr(
    ~ fasta2df(
      .x
    ),
    .id = "file_path"                        # new column “source” ← basename(.x)
  ) %>%
  clean_names() %>%
  select(-extraction_id) %>%
  filter(locus %in% loci_mito) %>%
  group_by(isolate) %>%                  # group by the column of interest
  filter(all(loci_mito %in% locus)) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  arrange(organism, isolate, locus)


loci_mito %>%
  map(
    ~ data_sanger_mito %>%
      filter(locus == .x) %>%
      df2fasta(
        ., 
        fasta_seq_names = c("organism", "isolate"), 
        fasta_file = str_c(
          "../output/multilocus_phylogeny/16s_coi/ischnura_luta_",
          str_to_lower(.x),
          "_groupisolates.fasta"
        )
      )
  )


