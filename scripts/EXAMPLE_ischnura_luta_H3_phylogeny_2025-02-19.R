#### INITIALIZE ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### USER DEFINED VARIABLES ####
library(pegas)
source("../../phylogenetic_tools/functions_sanger.R")
# https://github.com/Ph-IRES/phylogenetic_tools
source("../../phylogenetic_tools/phylogeny_functions.R")

#### Make Consensus Seqs from Fwd-Rev Curated AB1 ####

# this concatonates the curated Fwd and rev seqs for each sample.
# there is an assumption that your sequence file names begin with `SampleID_PrimerName_`.  If they do not, you need to rename them.
# The primer names below should match that in the file name
# the in_files variable is a regular expression for the path that specifies all of your files
# you name the out_file

processCuratedSANGER(
  in_files = "../output/sanger_curated_ab1_ischnura_luta_H3/^.*\\.scf$",
  out_file = "../output/sanger_curated_ab1_ischnura_luta_H3/ischnura_luta_H3_consensus_sequences.fasta",
  fwd_primer_name = "HEXAFIH3",
  rev_primer_name = "HEXARIH3"
)

#### Update GenBank Fasta Seq Names ####

# consult the instructions. You do a blast search on one of your sequences to get all related sequences in genbank.  You download them from the MSA viewer
# the TSV file is created by a bash script that runs on the fasta file you downloaded, `wrangle_blast_output.bash`
# the Fasta file is the one you downloaded from the blast MSA viewer
# you name the output file

renameBlastFastaSeqs(
  inTsvFile = "../output/sanger_curated_ab1_ischnura_luta_H3/blast_rbd_06_E1_all_619.tsv",
  inFastaFile = "../output/sanger_curated_ab1_ischnura_luta_H3/blast_rbd_06_E1_all_619.fasta",
  outFastaFile = "../output/sanger_curated_ab1_ischnura_luta_H3/blast_rbd_06_E1_all_619_renamed.fasta"
)

#### Filter GenBank Fasta Seqs to Keep Only Ingroup & Outgroup Seqs ####

# Use this to keep the ingroup  and outgroup seqs you want and exclude the seqs you don't want from your blast msa download
# if your outgroup is in your blast search, then use this to keep it in analysis, otherwise add the outgroup sequence with concatFastas()
# the Fasta file is the one you downloaded from the blast MSA viewer
# the ingroup and outgroup files have the accession numbers of the seuqences you want to keep, one per line
# you name the output file

filtered <- 
  filterFastaByAccession(
    fasta_file = "../output/sanger_curated_ab1_ischnura_luta_H3/blast_rbd_06_E1_all_619_renamed.fasta",
    ingroup_file = "../output/sanger_curated_ab1_ischnura_luta_H3/accession_list_damsel-dragonflies.txt",
    outgroup_file = "../output/sanger_curated_ab1_ischnura_luta_H3/accession_list_damsel-dragonflies_outgroup.txt",
    output_file = "../output/sanger_curated_ab1_ischnura_luta_H3/blast_rbd_06_E1_filtered_619_renamed.fasta"
  )


#### Concat Fastas, Deduplicate Haplotypes, Align ####

# does what the name implies, concatonates as many fastas as you list.
# you minimally want to add your consenus F&R seqs for each sample from processCuratedSANGER() and the blast MSA seqs you want in your tree
# you could also add the outgroup sequence here.  You should trim down the outgroup sequence to match your consensus seqs from processCuratedSANGER()
# you name the output file

concatFastas(
  inFilePaths = c(
    "../output/sanger_curated_ab1_ischnura_luta_H3/ischnura_luta_H3_consensus_sequences.fasta",
    "../output/sanger_curated_ab1_ischnura_luta_H3/blast_rbd_06_E1_filtered_619_renamed.fasta" #,
    # "../output/sanger_curated_ab1_ischnura_luta_H3/EU055464_dragonfly_Progomphus_borealis.fasta"
  ),
  outFilePath = "../output/sanger_curated_ab1_ischnura_luta_H3/ischnura_luta_H3_filtered_619.fasta"
)


#### REMOVE DUPLICATE HAPLOTYPES ####

# input file should be the output from concatFastas()
# you name the output file

uniqueSeqsFastaFile(
  inFilePath = "../output/sanger_curated_ab1_ischnura_luta_H3/ischnura_luta_H3_filtered_619.fasta",
  outFilePath = "../output/sanger_curated_ab1_ischnura_luta_H3/ischnura_luta_H3_filtered_619_haps.fasta"
)

#### ALIGN UNIQUE HAPLOTYPES ####

# input file should be the output from uniqueSeqsFastaFile()
# you name the output file

fasta <-
  alignFastaFile(
    inFilePath = "../output/sanger_curated_ab1_ischnura_luta_H3/ischnura_luta_H3_filtered_619_haps.fasta",
    outFilePath = "../output/sanger_curated_ab1_ischnura_luta_H3/ischnura_luta_H3_filtered_619_aligned.fasta"
  )

#### MAXIMUM LIKELIHOOD PHYLOGENY ####

# this takes a while to run, it runs modeltest before creating the bootstrapped ML tree
# you are specifying the outgroup with a regular expression that targets the name of the sequence inside of the file.  the accession number should suffice
# you can specify the number of bootstraps and the threshold value for "significance"

tree <-
  fasta %>%
  fasta2tree(
    n_bootstraps = 100,
    threshold_bootstraps = 50,
    my_outgroup = "EU055464"
  )

saveNewickTree(tree, 
               "../output/sanger_curated_ab1_ischnura_luta_H3/tree_ischnura_luta_filtered_619.nwk")

#### PLOT TREE ####
tree %>%
  plotGgTree(
    threshold_bootstraps = 50,
    tip_color = "red", 
    tip_pattern = "^rbd",
    show_node_id = FALSE
  )

#### PLOT HAPLOTYPE NETWORK ####

# here you want to start with a reduced number of sequences, genarally those for a species or closely related species.  
# no outgroups

sequences <- read.dna(file ="../output/sanger_curated_ab1_ischnura_luta_H3/ischnura_luta_H3_rbd_sequences_nooutgroup.fasta", format = "fasta")
haplo_data <- pegas::haplotype(sequences)
haplo_div <- hap.div(haplo_data)

hap_net <- haploNet(haplo_data)
plot(hap_net, size = attr(hap_net, "freq"), label = TRUE)
text(x = -3, y = -2, labels = paste("Hapl. Div.:", round(haplo_div, 2)), pos = 2)


# num_haplotypes <- length(attr(haplo_data, "index"))
# colors <- rainbow(num_haplotypes)
# plot(hap_net, size = attr(hap_net, "freq"), label = TRUE, col = colors)
# 
# greyscales <- gray(seq(0.2, 0.9, length.out = num_haplotypes))
# plot(hap_net, size = attr(hap_net, "freq"), label = TRUE, bg = greyscales)

