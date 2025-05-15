#' Convert an NCBI BLAST/FASTA alignment to a tidy TSV (and optional FASTA)
#'
#' @param infile  path to the *.aln* / FASTA file exported by NCBI MSA
#' @param outfile optional path to write the tab‑separated table;
#'                if NULL (default) the data.frame is returned
#' @param fasta_seq_names optional character vector `c()` of column names to combine (with '_') for FASTA headers;
#'                        if provided, a FASTA file will be written replacing the .tsv suffix of outfile with .fast
#'                        The column names are "country" and standard genbank field names
#' @param dedup   if TRUE, keep only the first row for each
#'                *database + accession* pair
#'                Takes precedence over \code{unique} and \code{rmdups}.
#' @param unique  if TRUE, drop later records whose complete
#'                *header + sequence* is identical to an earlier one
#'                Takes precedence over \code{rmdups}.
#' @param output_dups    if TRUE, emit only the rows that occur more than once
#'                (definition depends on \code{dedup} or \code{unique}).
#'                if FALSE, emit all rows passing \code{dedup}, \code{unique}, and \code{rmdups}
#' @param rmdups  if TRUE, remove *all* rows whose *database + accession*
#'                appears >1× after the previous filters
#'
#' @return invisibly, the resulting data.frame (also written to
#'         \code{outfile} if that argument is supplied)
#' @export
aln2tsv <- 
  function(
    infile,
    outfile = NULL,
    fasta_seq_names = NULL,
    dedup   = FALSE,
    unique  = FALSE,
    output_dups = FALSE,
    rmdups  = FALSE
  )
  {
    lines <- readLines(infile, warn = FALSE)
    
    # --- split into records ---------------------------------------------------
    recs <- list()
    header <- NULL
    seqbuf <- character()
    flush_rec <- function() {
      if (is.null(header)) return()
      recs[[length(recs) + 1]] <<- list(header = header,
                                        sequence = paste0(seqbuf, collapse = ""))
    }
    
    for (ln in lines) {
      ln <- sub("\r$", "", ln)               # strip DOS CR
      if (startsWith(ln, ">")) {             # new record
        flush_rec()
        header <- sub("^>", "", ln)
        seqbuf <- character()
      } else {
        seqbuf <- c(seqbuf, ln)
      }
    }
    flush_rec()
    
    # --- pass 1: parse each record -------------------------------------------
    all_keys <- character()
    parsed   <- lapply(recs, function(rec) {
      hdr  <- rec$header
      seq  <- rec$sequence
      
      # tokenise header
      tok  <- strsplit(hdr, " ")[[1]]
      idtk <- tok[1]
      db_acc <- strsplit(idtk, "\\|", fixed = FALSE)[[1]]
      database  <- db_acc[1]
      accession <- db_acc[2]
      
      # description (tokens until first '[')
      first_bracket <- which(startsWith(tok, "["))
      descr <- if (length(first_bracket)) {
        paste(tok[2:(first_bracket[1] - 1)], collapse = " ")
      } else {
        paste(tok[-1], collapse = " ")
      }
      
      # key=value brackets
      # tags <- regmatches(hdr, gregexpr("\\[[^][]+\\]", hdr))[[1]]
      # tags <- gsub("^\\[|\\]$", "", tags)

      # allow nested […] by using PCRE recursion (perl = TRUE)
      tags <- regmatches(
        hdr,
        gregexpr("\\[((?>[^\\[\\]]+|(?R))*)\\]", hdr, perl = TRUE)
      )[[1]]
      tags <- gsub("^\\[|\\]$", "", tags)
      
      kv <- if (length(tags)) {
        kvsplit <- strsplit(tags, "=", fixed = TRUE)
        setNames(lapply(kvsplit, `[`, 2), vapply(kvsplit, `[`, "", 1))
      } else list()
      
      all_keys <<- union(all_keys, names(kv))
      
      list(database = database,
           accession = accession,
           description = descr,
           tags = kv,
           sequence = seq,
           header = hdr)               # header kept for --unique comparison
    })
    
    # --- pass 2: build data.frame --------------------------------------------
    base_cols <- c("database", "accession", "description")
    tag_cols  <- sort(setdiff(all_keys, base_cols))
    col_order <- c(base_cols, tag_cols, "sequence")
    
    df <- as.data.frame(
      do.call(rbind, lapply(parsed, function(p) {
        row <- setNames(as.list(rep(NA_character_, length(col_order))), col_order)
        row$database    <- p$database
        row$accession   <- p$accession
        row$description <- p$description
        row$sequence    <- p$sequence
        if (length(p$tags)) row[names(p$tags)] <- p$tags
        row
      })),
      stringsAsFactors = FALSE
    )
    
    # --- pass 3: filtering logic ---------------------------------------------
    keep <- rep(TRUE, nrow(df))
    
    # uniqueness test key
    hdr_seq_key <- paste0(vapply(parsed, `[[`, "", "header"),
                          "\f",                         # unlikely delimiter
                          df$sequence)
    
    if (unique) {
      dup_uni <- duplicated(hdr_seq_key)
      if (!output_dups) keep[dup_uni] <- FALSE  # drop later perfect clones
    }
    
    # de‑dup / duplicate keys
    dbacc_key <- paste(df$database, df$accession, sep = "\f")
    
    if (dedup && !output_dups) {
      dup_dbacc <- duplicated(dbacc_key)
      keep[dup_dbacc] <- FALSE           # keep first only
    }
    
    if (output_dups) {
      # mark singletons for removal according to chosen key
      dup_by <- if (dedup) dbacc_key else hdr_seq_key
      tab <- table(dup_by)
      keep[tab[dup_by] == 1L] <- FALSE
    }
    
    if (rmdups) {
      tab <- table(dbacc_key[keep])      # after all previous filters
      keep[tab[dbacc_key] > 1L] <- FALSE # remove *all* copies of duplicates
    }
    
    # df_out <- df[keep, , drop = FALSE]
    
    df_out <-
      data.frame(
        lapply(df[keep, , drop = FALSE], function(col) {
          if (is.list(col)) {
            # each element should be a length‑1 character, so unlist
            unlist(col, use.names = FALSE)
          } else {
            col
          }
        }),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    
    # add "country" column which is needed by phylogenetic tools rename script
    df_out$country <- sub(":.*$", "", df_out$geo_loc_name)
    
    # --- output TSV--------------------------------------------------------------
    if (!is.null(outfile)) {
      write.table(df_out, outfile, sep = "\t",
                  quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
    
    # --- optional FASTA output ------------------------------------------------
    if (!is.null(fasta_seq_names)) {
      if (is.null(outfile)) stop("outfile must be specified to write FASTA")
      fasta_file <- sub("\\.tsv$", ".fasta", outfile)
      con <- file(fasta_file, "w")
      for (i in seq_len(nrow(df_out))) {
        hdr_vals <- df_out[i, fasta_seq_names, drop = TRUE]
        seq_name <- paste(hdr_vals, collapse = "_")
        writeLines(paste0(">", seq_name), con)
        writeLines(df_out$sequence[i], con)
      }
      close(con)
    }
    
    invisible(df_out)
  }


#' Write a FASTA file from an aln2tsv data.frame
#'
#' @param df               data.frame as returned by aln2tsv()
#' @param fasta_seq_names  character vector of column names to paste (with "_") for each FASTA header
#'                         If there are spaces in the name, they are replaced with "-"
#'                         Column names are "country" and GenBank field names
#' @param fasta_file       path to write the FASTA file (e.g. "output.fa")
#'
#' @return Invisibly returns the path to the written FASTA file.
#' @export

df2fasta <- function(df, fasta_seq_names, fasta_file) {
  # sanity check
  missing_cols <- setdiff(fasta_seq_names, names(df))
  if (length(missing_cols)) {
    stop("Columns not found in data frame: ", paste(missing_cols, collapse = ", "))
  }
  if (!"sequence" %in% names(df)) {
    stop("Data frame must have a 'sequence' column.")
  }
  
  con <- file(fasta_file, "w")
  on.exit(close(con), add = TRUE)
  
  for (i in seq_len(nrow(df))) {
    # extract header fields, replace spaces with dashes
    hdr_vals <- df[i, fasta_seq_names, drop = TRUE]
    hdr_vals <- gsub(" ", "-", as.character(hdr_vals), fixed = TRUE)
    
    seq_name <- paste(hdr_vals, collapse = "_")
    # write FASTA entry
    writeLines(paste0(">", seq_name), con)
    writeLines(df$sequence[i],       con)
  }
  
  invisible(fasta_file)
}

#############################################
##  ischnura_luta_read_fasta_to_dataframe.R
#############################################



## ---- 3.  Helper to convert one FASTA file to a tibble ----
fasta2df <-
  function(file_path) {
    
    dna <- readDNAStringSet(file_path)          # Biostrings object
    gene <- str_extract(basename(file_path),    # “16s”, “28s”, “coi”, “H3”, …
                        "(?<=_)[0-9A-Za-z]+(?=_consensus|$)")
    ## Species name (e.g. "Ischnura luta") from the part
    ## that precedes the first underscore in the filename
    species_raw <- str_extract(basename(file_path), "^[A-Za-z]+_[A-Za-z]+")
    organism    <- species_raw |>
      str_replace("_", " ") |>
      str_to_sentence()      # capitalise genus
    
    tibble(
      locus      = str_to_upper(gene),
      organism  = organism, 
      extraction_id        = names(dna),
      isolate   = str_remove(names(dna), "_E1$"), 
      sequence  = as.character(dna),
      length_bp = width(dna)
    )
  }

#############################################
##  ischnura_luta_read_fasta_to_dataframe_baseR.R
#############################################

## ---- base‑R FASTA → data.frame helper ----
fasta2df <- function(file_path) {
  ## 1. Read the entire FASTA file
  lines       <- readLines(file_path)
  
  ## 2. Identify headers (lines starting with ">") and strip ">"
  hdr_idx     <- grep("^>", lines)
  extraction_id <- sub("^>", "", lines[hdr_idx])
  
  ## 3. Collapse each block of sequence lines into one string
  seqs        <- character(length(hdr_idx))
  for (i in seq_along(hdr_idx)) {
    start_idx <- hdr_idx[i] + 1
    end_idx   <- if (i < length(hdr_idx)) hdr_idx[i + 1] - 1 else length(lines)
    seqs[i]   <- paste(lines[start_idx:end_idx], collapse = "")
  }
  
  ## 4. Sequence lengths
  length_bp   <- nchar(seqs)
  
  ## 5. Extract gene/locus (e.g. "16s","28s","coi","H3") from filename
  # gene        <- sub(".*_([0-9A-Za-z]+)(?:_consensus.*)?\\.fasta$",
  #                    "\\1",
  #                    basename(file_path))
  gene        <- sub("^.*_.*_([0-9A-Za-z]+)_consensus_sequences.*\\.fasta$",
                     "\\1",
                     basename(file_path))
  locus       <- toupper(gene)
  
  ## 6. Extract species ("Ischnura luta") from first two words of filename
  sp_raw      <- sub("^([A-Za-z]+_[A-Za-z]+).*", "\\1", basename(file_path))
  organism    <- sub("_", " ", sp_raw)
  
  ## 7. Derive isolate by dropping trailing "_E1" from the ID
  isolate     <- sub("_E1$", "", extraction_id)
  
  ## 8. Return a plain data.frame
  data.frame(
    locus         = locus,
    organism      = organism,
    extraction_id = extraction_id,
    isolate       = isolate,
    sequence      = seqs,
    length_bp     = length_bp,
    stringsAsFactors = FALSE
  )
}

