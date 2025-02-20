#!/bin/bash

# process fasta file from blast search msa to tidy tsv format

# relies on `parse_blast.awk`

DEFAULT_INFILE="../output/sanger_curated_ab1_ischnura_luta/blast_rbd_06_E1_500.fasta"

# Use $1 if it is not empty and is provided; otherwise, use the default path
INFILE=${1:-$DEFAULT_INFILE}

 cat $INFILE | \
	sed 's/^\(>.*$\)/\1@@@@/' | \
	tr "\n" "\t" | \
	sed 's/@@@@/\n/g' | \
	sed 's/>/\n>/' | \
	sed 's/\t//g' | \
	tail -n+2 | \
	paste - - | \
	grep -v "Query" | \
	sed "s/\([A-Za-z:,'0-9\.]\)\s\s*\([A-Za-z:,'0-9\.]\)/\1_\2/g" | \
	sed "s/\([A-Za-z:,'0-9\.]\)\s\s*\([A-Za-z:,'0-9\.]\)/\1,\2/g" | \
	./parse_blast.awk | \
	sed -e 's/^gb//' \
		-e 's/^ref//' \
		-e 's/^dbj//' \
		-e 's/^emb//' 