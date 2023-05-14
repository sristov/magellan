# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

import numpy as np
import json

def recode(data):

    columns_to_recode = ['ID','father','mother']
    unique_entries, unique_indices = np.unique(data[columns_to_recode],return_index=True)
    unique_entries = unique_entries[np.argsort(unique_indices)]

    num_digits = int(np.ceil(np.log10(len(unique_entries))))
    addition = 10**(num_digits-1)

    if int(np.ceil(np.log10(len(unique_entries)+addition))) > num_digits:
        addition += 1

    recoded_ids = {}
    idx = 1
    for entry in unique_entries:
        if entry == '0' or entry == '':
            recoded_ids[entry] = '0'
            continue
        else:
            recoded_ids[entry] = f"{idx+addition}"
        idx += 1

    old_ids = dict(zip(recoded_ids.values(), recoded_ids.keys()))

    new_data = data.replace({'ID': recoded_ids,'father': recoded_ids, 'mother': recoded_ids})
    unnamed_cols = [column for column in data.columns if 'Unnamed' in column]
    new_data.drop(columns=unnamed_cols,inplace=True)
    new_data["oldID"] = new_data['ID'].map(old_ids)

    with open('recoded_ids.json', 'w', encoding='utf8') as json_file:
        json.dump(recoded_ids, json_file, indent=4)#, sort_keys=True)

    return new_data
