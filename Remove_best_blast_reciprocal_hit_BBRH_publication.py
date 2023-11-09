#!/usr/bin/python

# Usage: python Remove_best_blast_reciprocal_hit_BBRH_publication.py All_proteins_METAG_NR_IMGV_CLUSTERS_DB_diamondblastp_out_edit_example.tsv

# import lib
import pandas as pd
import sys
import csv

# open input file
f = open(sys.argv[1])

# upload blast output edited with headers
df = pd.read_csv((f), delimiter='\t')

# create a column "BBRH" (best blast reciprocal hit) joinin qseqid qnd sseqid columns from blast
df['BBRH'] = df.apply(lambda row: ''.join(sorted([row['qseqid'], row['sseqid']])), axis=1)

# remove duplicates
no_BBRH = df.drop_duplicates('BBRH')


# print result of final dataframe removing the index (integer numeric column)
no_BBRH.to_csv(sys.argv[1] + '_no_BBRH.csv', index=False)