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

#### Wrangle aln Files ####

source("aln2tsv10.R")
source("../../phylogenetic_tools/multilocus_phylogeny_fuctions.R") # moved aln2tsv() to phylogenetic tools

source("df2fasta.R")

wrangled_aln <- 
  list.files(
    path       = "../output/multilocus_phylogeny/",
    pattern    = "\\.aln$",
    recursive  = TRUE,
    full.names = TRUE
  ) %>%
  set_names(.) %>%              # so .id is just the file name
  map_dfr(
    ~ aln2tsv(
      infile      = .x,
      unique      = TRUE,
      dedup       = FALSE,
      rmdups      = TRUE,
      output_dups = FALSE
    ),
    .id = "file_path"                        # new column “source” ← basename(.x)
  ) %>%
  clean_names() %>%
  select(
    file_path:accession,
    organism,
    isolate,
    specimen_voucher,
    gcode,
    sequence,
    everything()
  ) %>%
  filter(
    !str_detect(
      accession, 
      "Query_"
    )
  ) %>%
  mutate(
    locus = 
      str_extract(
        file_path %>% 
          str_to_lower(), 
        '16s|coi|28s|[Hh]3'       # the possible loci
      ) %>%
      str_to_upper(.), 
    .after = file_path
  ) %>%
  as_tibble()


# tsv_coi <- 
#   aln2tsv(
#     infile = "../output/multilocus_phylogeny/blast_rbd_23_E1_5000_coi.aln",
#     outfile = "../output/multilocus_phylogeny/blast_rbd_23_E1_5000_coi_uniq_rmdups.tsv",
#     unique = TRUE,
#     dedup = FALSE,
#     rmdups = TRUE,
#     output_dups = FALSE
#   )
# 
# tsv_16s <- 
#   aln2tsv(
#     infile = "../output/multilocus_phylogeny/blast_rbd_04_E1_5000_16s.aln",
#     outfile = "../output/multilocus_phylogeny/blast_rbd_04_E1_5000_16s_uniq_rmdups.tsv",
#     unique = TRUE,
#     dedup = FALSE,
#     rmdups = TRUE,
#     output_dups = FALSE
#   )
# 
# tsv_28s <- 
#   aln2tsv(
#     infile = "../output/multilocus_phylogeny/blast_rbd_17_E1_5000_28s.aln",
#     outfile = "../output/multilocus_phylogeny/blast_rbd_17_E1_5000_28s_uniq_rmdups.tsv",
#     unique = TRUE,
#     dedup = FALSE,
#     rmdups = TRUE,
#     output_dups = FALSE
#   )
# 
# tsv_h3 <- 
#   aln2tsv(
#     infile = "../output/multilocus_phylogeny/blast_rbd_05_E1_5000_h3.aln",
#     outfile = "../output/multilocus_phylogeny/blast_rbd_05_E1_5000_h3_uniq_rmdups.tsv",
#     unique = TRUE,
#     dedup = FALSE,
#     rmdups = TRUE,
#     output_dups = FALSE
#   )



#### Wrangle TSV Files ####

# wrangled_tsv <- 
#   list.files(
#     '../output', 
#     pattern = '\\.tsv$', 
#     recursive = TRUE, 
#     full.names = TRUE
#   ) %>%
#   # filter out unwanted files
#   str_subset(
#     'deprecated|qpcr', 
#     negate = TRUE
#   ) %>%
#   read_delim(
#     delim = '\t', 
#     id = 'file_path',
#     show_col_types = FALSE
#   ) %>%
#   mutate(
#     locus = 
#       str_extract(
#         file_path %>% 
#           str_to_lower(), 
#         '16s|coi|28s|H3'       # the possible loci
#       ) %>%
#       str_to_upper(.), 
#     .after = file_path
#   ) %>%
#   clean_names() %>%
#   distinct(
#     locus,
#     accession,
#     location,
#     organism,
#     specimen_voucher,
#     isolate,
#     .keep_all = TRUE
#   )

#### GET DUPLICATED ACCESSIONS ####

# duplicated accessions indicate that the sequences belong to a genome accession

accession_dupes <- 
  wrangled_aln %>%                     # <- replace with your object
  group_by(accession) %>%                  # group by the column of interest
  filter(n() > 1) %>%                      # keep groups with >1 row
  ungroup() %>%                                # optional: drop the grouping
  arrange(organism, isolate, accession)

accession_dupes %>%
  filter(
    str_detect(
      organism,
      "Ischnura"
    )
  )

#### GET DUPLICATED SPECIMEN VOUCHERS ####

accession_specvouch <-
  wrangled_aln %>%                     # <- replace with your object
  filter(!is.na(specimen_voucher)) %>%
  group_by(specimen_voucher) %>%                  # group by the column of interest
  filter(n() > 1) %>%                      # keep groups with >1 row
  ungroup() %>%                                # optional: drop the grouping
  arrange(organism, isolate, specimen_voucher) %>%    # keep rows with a voucher
  group_by(specimen_voucher, locus) %>%
  filter(n() == 1) %>%
  group_by(specimen_voucher) %>%                  # group by the column of interest
  filter(n() > 1)

#### GET DUPLICATED ISOLATES ####

isolate_dups <-
  wrangled_aln %>%                     # <- replace with your object
  filter(!is.na(isolate)) %>%
  group_by(isolate) %>%                  # group by the column of interest
  filter(n() > 1) %>%                      # keep groups with >1 row
  ungroup() %>%                                # optional: drop the grouping
  arrange(organism, isolate, locus) %>%
  # filter out 
  filter(isolate > 7) %>%
  filter(isolate != "B")  %>%
  filter(isolate != "v.2019") %>%
  filter(
    str_detect(
      isolate,
      "ishin11M"
    ) == FALSE
  )%>%
  filter(isolate != "fNotCel1") %>%
  filter(isolate != "AG129_WGS33_773") 

isolate_dups %>%
  filter(
    str_detect(
      organism,
      "Ischnura"
    )
  ) 


#### GET DUPLICATED ISOLATES for COI, 16S, H3 & 28S ####

loci_all <- c("H3", "COI", "16S", "28S")

isolate_dups_all_loci <-
  isolate_dups %>%
  group_by(organism, isolate) %>%
  filter(all(loci_all %in% locus)) %>%
  filter(n() > 3) %>% 
  ungroup()  %>%
  arrange(organism, isolate, locus) %>%
  pivot_wider(
    id_cols   = c(organism, isolate, country),
    names_from = locus,
    values_from = c(file_path, gcode, sequence, accession),
    names_sep  = "_"
  ) %>%
  distinct(
    organism,
    across(c(organism, starts_with("sequence_"))),
    .keep_all = TRUE
  ) %>%
  pivot_longer(
    cols      = matches("^(file_path|gcode|sequence|accession)_"),
    names_to  = c(".value", "locus"),
    names_pattern = "^(.*)_(.*)$"
  ) %>%
  arrange(organism, isolate, locus) %>%
  filter(!str_detect(isolate, "PL068"))  

loci_all %>%
  map(
    ~ isolate_dups_all_loci %>%
      filter(locus == .x) %>%
      df2fasta(
        df = ., 
        fasta_seq_names = c("organism", "isolate", "country"), 
        fasta_file = 
          str_c(
            "../output/multilocus_phylogeny/16s_coi_28s_h3/blast_",
            .x,
            "_uniq_rmdups_groupisolates.fa"
          )
        )
  )
        

  
#### GET DUPLICATED ISOLATES for 16S, H3 & 28S ####

loci_three <- c("H3", "16S", "28S")

isolate_dups_3_loci <-
  isolate_dups %>%
  filter(locus %in% loci_three)  %>%
  group_by(organism, isolate) %>%
  filter(all(loci_three %in% locus)) %>%
  filter(n() > 2) %>% 
  ungroup()  %>%
  arrange(organism, isolate, locus) %>%
  pivot_wider(
    id_cols   = c(organism, isolate, country),
    names_from = locus,
    values_from = c(file_path, gcode, sequence, accession),
    names_sep  = "_"
  ) %>%
  distinct(
    organism,
    across(c(organism, starts_with("sequence_"))),
    .keep_all = TRUE
  ) %>%
  pivot_longer(
    cols      = matches("^(file_path|gcode|sequence|accession)_"),
    names_to  = c(".value", "locus"),
    names_pattern = "^(.*)_(.*)$"
  ) %>%
  arrange(organism, isolate, locus) %>%
  filter(!str_detect(isolate, "BYU_PL068")) 


loci_three %>%
  map(
    ~ isolate_dups_3_loci %>%
      filter(locus == .x) %>%
      df2fasta(
        df = ., 
        fasta_seq_names = c("organism", "isolate", "country"), 
        fasta_file = 
          str_c(
            "../output/multilocus_phylogeny/16s_28s_h3/blast_",
            .x,
            "_uniq_rmdups_groupisolates.fa"
          )
      )
  )

#### GET DUPLICATED ISOLATES for H3 & 28S ####

loci_nuc <- c("H3", "28S")

isolate_dups_nuc_loci <-
  isolate_dups %>%
  filter(locus %in% loci_nuc)  %>%
  group_by(organism, isolate) %>%
  filter(all(loci_nuc %in% locus)) %>%
  filter(n() > 1) %>% 
  ungroup() %>%
  arrange(organism, isolate, locus) %>%
  pivot_wider(
    id_cols   = c(organism, isolate, country),
    names_from = locus,
    values_from = c(file_path, gcode, sequence, accession),
    names_sep  = "_"
  ) %>%
  distinct(
    organism,
    across(c(organism, starts_with("sequence_"))),
    .keep_all = TRUE
  ) %>%
  pivot_longer(
    cols      = matches("^(file_path|gcode|sequence|accession)_"),
    names_to  = c(".value", "locus"),
    names_pattern = "^(.*)_(.*)$"
  ) %>%
  arrange(organism, isolate, locus) %>%
  filter(!str_detect(isolate, "BYU_PL068")) %>%
  filter(!str_detect(isolate, "BEA526") )%>%
  filter(!str_detect(isolate, "BEA527"))

loci_nuc %>%
  map(
    ~ isolate_dups_nuc_loci %>%
      filter(locus == .x) %>%
      df2fasta(
        df = ., 
        fasta_seq_names = c("organism", "isolate", "country"), 
        fasta_file = 
          str_c(
            "../output/multilocus_phylogeny/28s_h3/blast_",
            .x,
            "_uniq_rmdups_groupisolates.fa"
          )
      )
  )


#### GET DUPLICATED ISOLATES for COI & 16S ####

loci_mito <- c("COI", "16S")

isolate_dups_mito_loci <-
  isolate_dups %>%
  filter(locus %in% loci_mito)  %>%
  group_by(organism, isolate) %>%
  filter(all(loci_mito %in% locus)) %>%
  filter(n() > 1) %>% 
  ungroup() %>%
  arrange(organism, isolate, locus) %>%
  pivot_wider(
    id_cols   = c(organism, isolate, country),
    names_from = locus,
    values_from = c(file_path, gcode, sequence, accession),
    names_sep  = "_"
  ) %>%
  distinct(
    organism,
    across(c(organism, starts_with("sequence_"))),
    .keep_all = TRUE
  ) %>%
  pivot_longer(
    cols      = matches("^(file_path|gcode|sequence|accession)_"),
    names_to  = c(".value", "locus"),
    names_pattern = "^(.*)_(.*)$"
  ) %>%
  arrange(organism, isolate, locus) %>%
  filter(!str_detect(isolate, "PL068"))

loci_mito %>%
  map(
    ~ isolate_dups_mito_loci %>%
      filter(locus == .x) %>%
      df2fasta(
        df = ., 
        fasta_seq_names = c("organism", "isolate", "country"), 
        fasta_file = 
          str_c(
            "../output/multilocus_phylogeny/16s_coi/blast_",
            .x,
            "_uniq_rmdups_groupisolates.fa"
          )
      )
  )

#### Get caccession#### Get counts for what specimen vouchers are missing ####
# count(wrangled_tsv, 
#       locus, 
#       missing_voucher = is.na(specimen_voucher)) %>%
#   summarise(prop_missing = n[missing_voucher] / sum(n),
#             .by = locus)
# 
# filter(wrangled_tsv, is.na(specimen_voucher))
# 
# filter(wrangled_tsv, is.na(specimen_voucher)) %>%
#   select(accession) %>%
#   distinct
# n_distinct(wrangled_tsv$accession)
# 
# wrangled_tsv %>%
#   filter(!is.na(specimen_voucher)) %>%
#   count(locus, 
#         specimen_voucher) %>%
#   filter(n > 1)
# 
# #### Identify Accession numbers across loci in same sample
# wrangled_tsv %>%
#   select(locus, Accession, Specimen_Voucher) %>%
#   pivot_wider(names_from = locus, 
#               values_from = Accession)
# 
# 
# count(wrangled_tsv, Specimen_Voucher) %>% filter(n > 1)
# 
# 
# wrangled_tsv %>%
#   select(locus, Accession, Specimen_Voucher) %>%
#   filter(n() > 1, .by = Specimen_Voucher) %>%
#   filter(!is.na(Specimen_Voucher)) %>%
#   pivot_wider(names_from = locus, 
#               values_from = Accession) 
# 
# 
# wrangled_tsv %>%
#   select(locus, Accession, Specimen_Voucher) %>%
#   filter(!is.na(Specimen_Voucher)) %>%
#   filter(n() > 1, .by = Specimen_Voucher)  %>%
#   pivot_wider(names_from = locus, 
#               values_from = Accession,
#               values_fn = ~str_c(.x, collapse = '; ')) %>% View
# 
# 
# wrangled_tsv %>%
#   filter(Specimen_Voucher == 'Coen1') %>%
#   pull(Isolate)
