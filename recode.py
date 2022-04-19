# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hršak and Strahil Ristov

import argparse
import pandas as pd
import numpy as np
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='input pedigree file')
parser.add_argument('--output', dest='output', default='',help='filename for recoded pedigree file')

args = parser.parse_args()

if args.output == '':
    output = os.path.splitext(args.filename)[0] + '_recoded.csv'
else:
    output = args.output

data = pd.read_csv(args.filename, dtype=str).fillna('')

columns_to_recode = ['ID','father','mother']
unique_entries, unique_indices = np.unique(data[columns_to_recode],return_index=True)
unique_entries = unique_entries[np.argsort(unique_indices)]

recoded_ids = {}
idx = 1
for entry in unique_entries:
    if entry == '0' or entry == '':
        recoded_ids[entry] = '0'
        continue
    else:
        recoded_ids[entry] = str(idx)
    idx += 1

data.replace({'ID': recoded_ids,'father': recoded_ids, 'mother': recoded_ids},inplace=True)
unnamed_cols = [column for column in data.columns if 'Unnamed' in column]
data.drop(columns=unnamed_cols,inplace=True)
data.to_csv(output,index=False)
