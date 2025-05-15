# Create Multilocus Phylogeny



## Brief Instructional Notes

1. Blast target species sequence(s) and return 5000 results on website
	* view with msa viewer in blast results page and download fasta


2. Wrangle the `aln` fasta files from NCBI BLAST MSA with `aln2tsv07.awk` (_NOTE: Step 2 replaced by step 3, but I'm not sure about the `comm` bit)
	* This script replaces `wrangle_blast_output.bash` and `parse_blast.awk`
	* I found that multicopy genes from genomes might have same db and accession but different sequences
	
	```bash
	# isolate records that are not completely identical, but are identical in the db and accession columns
	# run from scripts dir
	comm -3 <(sort -u ../output/multilocus_phylogeny/test07_dedup_dups.tsv) <(sort -u ../output/multilocus_phylogeny/test07_dups.tsv) > ../output/multilocus_phylogeny/test07_difference_dups_dedup-dups.tsv
	```
	
	* I decided that we want to uniq identical rows, then remove all rows that have same db & accession.  In the future, we could try to identify which sequences produce viable protiens to select a seq to keep.
	
	```bash
	# run from scripts dir
	gawk -f aln2tsv08.awk --unique --rmdups ../output/multilocus_phylogeny/blast_rbd_05_E1_5000_h3.aln > ../output/multilocus_phylogeny/test08_uniq_rmdup.tsv
	```
	
3. Wrangle the `aln` fasta files from NCBI BLAST MSA with `aln2tsv*.R` in `prj_rotablue_barcoding\scripts/wrangle_blast_results_multilocus.R`

	* I converted the awk script to an R script `aln2tsv*.R`.  
		* now after downloading `*.aln` from NCBI MSA, they are read directly into R without any bash or awk scripts
		* the `aln2tsv()` function is called in `wrangle_blast_results_multilocus.R`
		* I copied it to `phylogenetic_tools/multilocus_phylogeny_functions.R`
	* In `prj_rotablue_barcoding\scripts` I call this function in `wrangle_blast_results_multilocus.R`
	
4. Wrangle the tsv files in R with `wrangle_blast_results_multilocus.R`
	* the goal here is to group sequences across loci by individual id, such that we know each accession for an individual
	* I converted the awk script to an R script `aln2tsv*.R` that has a function, `aln2tsv()`. It handles steps 2 and 3.  
		* this negates step 2 above, now after downloading `*.aln` from NCBI MSA, they are read directly into this R function
		
5. Wrangle the consensus fasta files from our sanger sequencing in R with `wrangle_blast_results_multilocus.R`
	* 
