# scripts and data to perform some analysis in virome data from the study conducted at the Lemonade Creek at YNP (USA)

Files:

- Taxonomic analysis
Taxonomic annotated example file: Creek_1_all_vir_scaffold_final_classification_edited_example_file.tsv

script to get the most frequent taxonomy consensus: get_most_freq_tax_consensus.py

Diamond BLASTp output file with edition to include the scaffold name: CreekBiofilm_1_assembly_1K_virus_proteins.faa.diamond_blastp_RefSeqComplete.outfmt6.gz_TOPHITS_edit_example.tsv

script to classify scaffold based on protein matches where a TRUE viral scaffold will contain more than 60% best hits with viruses
Get_viral_scaffold_approx_class_V3_publication.py

- Network and clustering analysis
All against all Diamond BLASTp example file as an input to a protein similarity network analysis
All_proteins_METAG_NR_IMGV_CLUSTERS_DB_diamondblastp_out_edit_example.tsv

script to remove best blast reciprocal hits (BBRH)
Remove_best_blast_reciprocal_hit_BBRH_publication.py

- Horizontal gene transfer (HGT) analysis
script to run inside python environment using an diamond BLASTp result, organized according to viral class
script_to_run_HGT_index_analysis_publication.sh

Input file example to run HGT index analysis "script_to_run_HGT_index_analysis_publication.sh" (other files are too large to upload but please contact me for other viral classes files)
List_Arfiviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv

- Compositional analysis
Genomic_protein_composition_analysis_publication.sh
scripts to run composition analysis using as input protein and nucleotide fasta sequences in the R environment using seqinr package

Input file example to run composition analysis using genomic (nucleotides *nucl.fna) and protein data (*prot.faa)
bin.750_prodigal_nucl.fna
bin.750_prodigal_prot.faa




