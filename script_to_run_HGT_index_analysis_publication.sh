# script to run inside python environment using an diamond BLASTp result, organized according to viral class

# target = non_alien
# nontarget = alien

# in python environment
python

# import libraries
import pandas as pd
import numpy
import numpy as np
import sys
import csv

# input csv blast output table (needs to be custom to your blast output data)
data = pd.read_csv("List_Arfiviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Caudoviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Faserviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Herviviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
# nothing for Malgrandaviricetes
data = pd.read_csv("List_Maveriviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Megaviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Pokkesviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Polintoviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
# there is no non-alien sequence for Resentoviricetes
# data = pd.read_csv("List_Resentoviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Revtraviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Tectiliviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Tokiviricetes_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)
data = pd.read_csv("List_Unclassified_vir_vir_scaff_prot.faa.diamond_blastp_RefSeqComplete.outfmt6.tsv", names=['qseqid_scaff', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle', 'staxids', 'sscinames', 'sphylums', 'skingdoms', 'sskingdoms'], delimiter='\t', index_col = 0)

# insert new column with category (non_alien and alien) according to data in sskingdoms (or other higher lineage) column and according with desired target group (Eukaryota or Bacteria etc.) and nontarget (Viruses or Bacteria etc.)
# here I am analyzing virus contigs from metagenomes, therefore they are non-alien and everything outside viruses is alien
data.loc[data.sskingdoms.str.contains("Viruses", na=False), 'category'] = 'non_alien'

# here define who should be the alien (non-target) group to investigate if there was HGT with your genome (non_alien) with this group (alien)
data.loc[data.sskingdoms.str.contains("Eukaryota|Bacteria|Archaea",na=False), 'category'] = 'alien'

## here is important to remove a selected taxa from the same genome you are analyzing, if not the best hit would be from the same group
# remove taxa or group which is same from genome of interest according to column (for instance here sphylums column)
# to remove multiple type for ie. ("Unknown taxid|Cyanidiococcus|Cyanidium|Cyanidioschyzon|Galdieria")
# for metag viruses I am removing all from the same Phylum

data = data[data.sphylums.str.contains("Cressdnaviricota") == False] # Arfiviricetes_vir_vir_scaff_prot
data = data[data.sphylums.str.contains("Uroviricota") == False] # Caudoviricetes_vir_vir_scaff_prot
data = data[data.sphylums.str.contains("Hofneiviricota") == False] # List_Faserviricetes
data = data[data.sphylums.str.contains("Peploviricota") == False] # Herviviricetes_vir_vir_scaff_prot
data = data[data.sphylums.str.contains("Preplasmiviricota") == False] # Maveriviricetes_vir_vir_scaff_prot
data = data[data.sphylums.str.contains("Nucleocytoviricota") == False] # Megaviricetes_vir_vir_scaff_prot
data = data[data.sphylums.str.contains("Nucleocytoviricota") == False] # Pokkesviricetes_vir_vir_scaff_prot
data = data[data.sphylums.str.contains("Preplasmiviricota") == False] # Polintoviricetes_vir_vir_scaff_prot
#data = data[data.sphylums.str.contains("Duplornaviricota") == False] # Resentoviricetes_vir_vir_scaff_prot # there is no non-alien sequence
data = data[data.sphylums.str.contains("Artverviricota") == False] # Revtraviricetes_vir_vir_scaff_prot
data = data[data.sphylums.str.contains("Preplasmiviricota") == False] # Tectiliviricetes_vir_vir_scaff_prot
data = data[data.sphylums.str.contains("Taleaviricota") == False] # Tokiviricetes_vir_vir_scaff_prot
# for unclassified I can't remove a single or all viruses so maybe remove the ones that are from major annotated groups (decided to remove all the major phylae)
data = data[data.sphylums.str.contains("Peploviricota|Nucleocytoviricota|Uroviricota|Cressdnaviricota|Taleaviricota|Preplasmiviricota") == False] # For unclassified viruses

# get best Evalues
BEST_EV= pd.pivot_table(data, index= 'qseqid', columns='category', values='evalue', aggfunc='min')

# replace NaN by 1 (in case there is no target or nontarget
BEST_EV_REPLACE= BEST_EV.replace(np.nan, 1, regex=True)

# get best bitscores
BEST_BIT= pd.pivot_table(data, index= 'qseqid', columns='category', values='bitscore', aggfunc='max')

# replace NaN by 1 (in case there is no target or nontarget
BEST_BIT_REPLACE= BEST_BIT.replace(np.nan, 0, regex=True)

# Original formula to calculate Alien Index (AI) - AI = log((Best E-value for Metazoa) + e-200) - log((Best E-value for Non-Metazoa) + e-200).
"""Based on AI, genes are classified as foreign (AI>30), indeterminate (0<AI<30), or target (AI<0) ... Entries with incomplete taxonomy,
    such as X-ray structures, or hits belonging to the same phylum as the query sequence (
    i.e Rotifera, filter_out_tax_id or Nematoda), were excluded from the analysis - here in my script now the Metazoa or target will be non-alien while the Non-Metazoa or foreign will be alien"""

BEST_EV_REPLACE["Alien index (AI)"] = (numpy.log((BEST_EV_REPLACE["non_alien"])+1e-200) - (numpy.log((BEST_EV_REPLACE["alien"]+ 1e-200))))

# Calculate HGT (hU) index = (best_bitscore_nontarget - best_bitscore_target) / paper formula HGT index (hU) = (best non-metazoan - best metazoan)
"""hU is a measure of how well a given sequence matches to one set of taxa (eg. Metazoa) relative to another, mutually exclusive set of taxa (eg. non-Metazoa). It uses best-hit bitscores and is defined as: (best-hit bitscore for OUTGROUP) - (best-hit bitscore for INGROUP). See Boschetti et al. 2012 for more details. hU >= 30"""
"""HGT index, hU, defined as the difference between the “bitscore” (i.e. score in bits) of the best non-metazoan match and the bitscore of the best metazoan match in the database. from https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003035"""

BEST_BIT_REPLACE["HGT index (hU)"] = BEST_BIT_REPLACE["alien"] - BEST_BIT_REPLACE["non_alien"]

# added here to keep alien and non_alien AI and alien and non_alien hU
BEST_EV_REPLACE.rename(columns={'alien':'alien_AI','non_alien':'non_alien_AI'}, inplace=True)
BEST_BIT_REPLACE.rename(columns={'alien':'alien_hU','non_alien':'non_alien_hU'}, inplace=True)

# concat two objects from AI and Hu (BEST_EV and BEST_BIT)
concat_BEST_EV_BIT = pd.concat([BEST_EV_REPLACE,BEST_BIT_REPLACE], axis=1)


# get best Evalue for indexed sequences (Query=query name)
Best_Evalues1= pd.pivot_table(data, index= 'qseqid', values='evalue', aggfunc='min')

# get best Bitscore for indexed sequences (Query=query name)
Best_Bitscore1= pd.pivot_table(data, index= 'qseqid', values='bitscore', aggfunc='max')

# UPDATED OCTOBER 23 get best (first) full taxa and protein name (stitle) for indexed sequences (Query=query name)
Best_Evalues3= pd.pivot_table(data, index= 'qseqid', values=['sscinames','staxids', 'sphylums','skingdoms', 'sskingdoms', 'stitle'], aggfunc='first')

# group by species name get species name (first)
Best_Evalues3= Best_Evalues3.groupby('qseqid').first()

# concat tables of alien index, HGT index, e-value and species name # here updated droping the drop_met_nonmet
cols_final = pd.concat([concat_BEST_EV_BIT, Best_Evalues1, Best_Bitscore1, Best_Evalues3], axis=1)

# change col names 'species name'
cols_final.rename(columns={'sscinames':'Top sscinames evalue','alien_AI':"alien_AI_best_hit",'non_alien_AI':'non_alien_AI_best_hit','alien_hU':'alien_hU_best_hit','non_alien_hU':'non_alien_hU_best_hit'}, inplace=True)

# print final table final - AI = TRUE if foreign (AI>30) - hU = TRUE if >=30
# in some cases update values for AI based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9678320/ "We found that selecting genes with AI > 10 represented an optimal balance between sensitivity and precision"
cols_final.loc[cols_final['Alien index (AI)'] < 30, 'AI > 30'] = 'False'
cols_final.loc[cols_final['Alien index (AI)'] > 30, 'AI > 30'] = 'True'

cols_final.loc[cols_final['HGT index (hU)'] < 30, 'hU >= 30'] = 'False'
cols_final.loc[cols_final['HGT index (hU)'] >= 30, 'hU >= 30'] = 'True'

# UPDATED OCTOBER 23 reorder cols
cols_final = cols_final[['alien_AI_best_hit','non_alien_AI_best_hit','Alien index (AI)','AI > 30','alien_hU_best_hit','non_alien_hU_best_hit','HGT index (hU)','hU >= 30', 'evalue', 'bitscore', 'Top sscinames evalue', 'staxids', 'sphylums','skingdoms', 'sskingdoms', 'stitle']]

# Replace ".0" to empty "" in the "staxids" column 
cols_final["staxids"] = cols_final["staxids"].astype(str).str.rstrip(".0")

# Replace ".0" to empty "" in the "bitscore" column 
cols_final["bitscore"] = cols_final["bitscore"].astype(str).str.rstrip(".0")

#### print final tables
cols_final.to_csv('Arfiviricetes_less_Cressdnaviricota_HGT_indexes_october.csv')
cols_final.to_csv('Caudoviricetes_less_Uroviricota_HGT_indexes_october.csv')
cols_final.to_csv('Faserviricetes_less_Hofneiviricota_HGT_indexes.csv')
cols_final.to_csv('Herviviricetes_less_Peploviricota_HGT_indexes.csv')
cols_final.to_csv('Maveriviricetes_less_Preplasmiviricota_HGT_indexes_october.csv')
cols_final.to_csv('Megaviricetes_less_Nucleocytoviricota_HGT_indexes_october.csv')
cols_final.to_csv('Pokkesviricetes_less_Nucleocytoviricota_HGT_indexes_october.csv')
cols_final.to_csv('Polintoviricetes_less_Preplasmiviricota_HGT_indexes_october.csv')
cols_final.to_csv('Revtraviricetes_less_Artverviricota_HGT_indexes_october.csv')
cols_final.to_csv('Tectiliviricetes_less_Preplasmiviricota_HGT_indexes_october.csv')
cols_final.to_csv('Tokiviricetes_less_Taleaviricota_HGT_indexes_october.csv')
cols_final.to_csv('Unclassified_all_HGT_indexes.csv')
cols_final.to_csv('Unclassified_less_major_phylum_HGT_indexes_october.csv')