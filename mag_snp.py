# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

import numpy as np
import pandas as pd
import dask.dataframe as dd
import sys
import json
from os import path
import xgboost as xgb
import tensorflow as tf
import pickle as pkl
from tkinter import messagebox, simpledialog


def analyze_SNP(final_table, snp_map, SNPs_310, samples_population,gc_score_thresh,gt_score_thresh,snp_position,mode):

    gc_score_thresh = float(gc_score_thresh)
    gt_score_thresh = float(gt_score_thresh)

    if snp_position == "AUTO":
        snp_map = snp_map.loc[snp_map['Chromosome'].str.isdigit()]
    else:
        snp_map = snp_map.loc[snp_map['Chromosome'] == snp_position]

    snp_map = snp_map.compute()
    final_table = final_table.loc[final_table["SNP Name"].isin(snp_map["Name"])]
    final_table = final_table.compute()
    final_table = final_table.merge(snp_map, left_on="SNP Name", right_on="Name", how="right",validate='many_to_one')

    final_table.drop(columns=['Name'],inplace=True)
    renaming_dict = {'SNP Name': 'SNP_Name', 
                     'Sample ID': 'Sample_ID',
                     'GC Score': 'GCScore', 
                     'Allele1 - Forward': 'Allele1_Forward',
                     'Allele2 - Forward': 'Allele2_Forward',
                     'Allele1 - Top':'Allele1_Top',
                     'Allele2 - Top':'Allele2_Top',
                     'Chromosome':'Chr',
                     'GenTrain Score':'GTScore'}
    final_table.rename(columns=renaming_dict,inplace=True)

    final_table = pd.merge(final_table, SNPs_310, on="SNP_Name", how="inner")

    final_table = pd.merge(final_table, samples_population, on="Sample_ID", how="inner")

    if gc_score_thresh > 0.0:
        final_table = final_table.loc[final_table['GCScore'] >= gc_score_thresh]
    if gt_score_thresh > 0.0:
        final_table = final_table.loc[final_table['GTScore'] >= gt_score_thresh]

    taurus_subset = final_table[["Population","Sample_ID","SNP_Name", "Chr", "Position", "Allele1_Forward","Allele2_Forward","GCScore","GTScore"]]

    values_filter = (taurus_subset['Position'] == '#REF') | (taurus_subset['Position'] == '#REF!') | (taurus_subset['Position'] == '#Ref') | (taurus_subset['Position'] == '#Ref!')
    taurus_subset = taurus_subset.loc[~values_filter]
    values_filter = (taurus_subset['Position'] == '') | (taurus_subset['Position'] == '0')
    taurus_subset = taurus_subset.loc[~values_filter]

    values_filter = (taurus_subset['Allele1_Forward'] == "-") | (taurus_subset['Allele1_Forward'] == "?") | (taurus_subset['Allele2_Forward'] == "-") | (taurus_subset['Allele2_Forward'] == "?")

    taurus_subset = taurus_subset.loc[~values_filter]

    taurus_LGEN = taurus_subset[["Population","Sample_ID","SNP_Name","Allele1_Forward","Allele2_Forward"]]
    taurus_LGEN.sort_values(by=['SNP_Name','Population','Sample_ID'])

    
    if mode == "gui":
        basename = simpledialog.askstring(
            "Save files", "Give name for output files.", initialvalue="taurus")
    elif mode == "cl":
        basename = "taurus"
    taurus_LGEN.to_csv(f"{basename}.lgen",header=False,index=False,sep=' ')

    taurus_fam = taurus_subset.drop_duplicates(subset="Sample_ID", ignore_index=True)

    taurus_fam = taurus_fam[["Population", "Sample_ID"]]

    taurus_fam["sireID"] = -9
    taurus_fam["damID"] = -9
    taurus_fam["sex"] = 1
    taurus_fam["phenotype"] = -9

    taurus_fam.to_csv(f"{basename}.fam",header=False,index=False,sep=' ')

    #########################################################
    # creating MAP file - columns: chromsome snpname cmunit bp_position
    ########################################################

    tMAP = taurus_subset.drop_duplicates(subset="SNP_Name", ignore_index=True)

    tMAP = tMAP[['Chr','SNP_Name','Position']]

    tMAP['genDist'] = tMAP['Position'].astype(int).div(1000000)
    tMAP['genDist'] = tMAP['genDist'].fillna(0)
    tMAP.sort_values(by=['Position'],inplace=True)

    tMAP[['Chr', 'SNP_Name', 'genDist', 'Position']].to_csv(f"{basename}.map", header=False, index=False, sep=' ')

    proba = final_table[["Sample_ID" ,"Allele1_Forward","Position"]]

    return taurus_subset, proba
    

def save_fasta(filename,proba):

    proba.Position = proba.Position.astype(int)
    _, idx = np.unique(proba['Sample_ID'], return_index=True)
    sample_IDs = list(proba['Sample_ID'].iloc[list(np.sort(idx))])

    with open(filename,'w') as fasta_file:
        for count, sample_ID in enumerate(sample_IDs):
            fasta_file.write(f">{sample_ID}\n")
            subset = proba.loc[proba['Sample_ID'] == sample_ID][['Allele1_Forward','Position']]
            subset.sort_values(by=['Position'],inplace=True)
            fasta_file.write("".join(subset.Allele1_Forward))
            if count < (len(sample_IDs) - 1):
                fasta_file.write("\n")

def check_mutations(taurus_subset,mode):
    here = path.abspath(path.dirname(__file__))
    diseasePath = path.join(here,"model/snp2disease.json")
    humanPath = path.join(here,"model/snp2human.json")
    bovinePath = path.join(here,"model/snp2bovine.json")
    try:
        with open(diseasePath, "r") as snp2disease_file:
            snp2disease = json.load(snp2disease_file)
    except FileNotFoundError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"File {diseasePath} not found!")
            return
        elif mode == "cl":
            sys.stderr.write(f"ERROR: File {diseasePath} not found!\n")
            return
    except PermissionError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"Not allowed to read {diseasePath}.")
            return
        elif mode == "cl":
            sys.stderr.write(f"ERROR: Not allowed to read {diseasePath}.\n")
            return

    try:
        with open(humanPath, "r") as snp2human_file:
            snp2human = json.load(snp2human_file)
    except FileNotFoundError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"File {humanPath} not found!")
            return
        elif mode == "cl":
            sys.stderr.write(f"ERROR: File {humanPath} not found!\n")
            return
    except PermissionError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"Not allowed to read {humanPath}.")
            return
        elif mode == "cl":
            sys.stderr.write(f"ERROR: Not allowed to read {humanPath}.\n")
            return

    try:
        with open(bovinePath, "r") as snp2bovine_file:
            snp2bovine = json.load(snp2bovine_file)
    except FileNotFoundError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"File {bovinePath} not found!")
            return
        elif mode == "cl":
            sys.stderr.write(f"ERROR: File {bovinePath} not found!\n")
            return
    except PermissionError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"Not allowed to read {bovinePath}.")
            return
        elif mode == "cl":
            sys.stderr.write(f"ERROR: Not allowed to read {bovinePath}.\n")
            return

    # detect deleterious mutation

    #indicated_snps = taurus_subset.loc[~taurus_subset["SNP_Name"].str.startswith("SNP")]
    indicated_snps = taurus_subset.loc[taurus_subset["SNP_Name"].isin(snp2disease.keys())]
    indicated_snps = indicated_snps.loc[indicated_snps["SNP_Name"].str[-1] == indicated_snps["Allele1_Forward"]]

    with open ("MUTATIONS_ALERT.TXT",'w') as mutations_file:
        if len(indicated_snps) == 0:
            mutations_file.write("NO DELETERIOUS MUTATIONS FOUND IN THE POPULATION.")
        else:
            mutations_file.write("DELETERIOUS MUTATIONS FOUND IN THE POPULATION:\n\n")
            mutations_file.write("ID            ")
            mutations_file.write("|SNP             ")
            mutations_file.write("|POSITION IN CATTLE   ")
            mutations_file.write("|POSITION IN HUMAN    |DISEASE\n")
            mutations_file.write(103*"-"+"\n")
            for _, row in indicated_snps.iterrows():
                mutations_file.write(f"{row['Sample_ID']:<14}")
                mutations_file.write(f"|{row['SNP_Name']:<16}")
                try:
                    mutations_file.write(f"|{snp2bovine[row['SNP_Name']]:<21}")
                except KeyError:
                    mutations_file.write("|NaN                  ")
                try:
                    mutations_file.write(f"|{snp2human[row['SNP_Name']]:<17}")
                except KeyError:
                    mutations_file.write("|NaN              ")
                try:
                    mutations_file.write(f"    |{snp2disease[row['SNP_Name']]}\n")
                except KeyError:
                    mutations_file.write("    |unknown\n")

def read_fasta(fasta_name,mode):
    sample2seq = {}
    try:
        with open(fasta_name, 'r') as fasta_file:
            for line in fasta_file:
                if line[0] == ">":
                    sample_id = line[1:].strip()

                else:
                    snp_sequence = line.strip()
                    sample2seq[sample_id] = [char.capitalize() for i, char in enumerate(snp_sequence)]
    except FileNotFoundError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"ERROR: Fasta file {fasta_name} found.")
        elif mode == "cl":
            sys.stderr.write(f"ERROR: Fasta file {fasta_name} found.\n")
    except PermissionError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"ERROR: You are not allowed to read {fasta_name}.")
        elif mode == "cl":
            sys.stderr.write(f"ERROR: You are not allowed to read {fasta_name}.\n")
    return sample2seq

letter2ohe = {"A": [1,0,0,0],
              "C": [0,1,0,0],
              "G": [0,0,1,0],
              "T": [0,0,0,1],
              "a": [1,0,0,0],
              "c": [0,1,0,0],
              "g": [0,0,1,0],
              "t": [0,0,0,1],
              "-": [0,0,0,0],
              "N": [0,0,0,0],
              "n": [0,0,0,0],
              "?": [0,0,0,0],}

class SNPClassifier():

    def __init__(self,xgbModelPath,id2seq):
        self.here = path.abspath(path.dirname(__file__))
        self.xgbModelPath = path.join(self.here,"model",xgbModelPath)
        self.id2seq = id2seq
        pass

    def load_xgb(self):
        self.xgbModel = xgb.XGBClassifier()
        self.xgbModel.load_model(self.xgbModelPath)
        pass

    def predict_main_haplogroup(self,mode):
        int2maingroup = {0: "I",
                         1: "P",
                         2: "Q",
                         3: "R",
                         4: "T"}

        self.id2pred = {}
        for id in self.id2seq:
            encoded = np.array([letter2ohe[char] for char in self.id2seq[id]])
            (i, j) = encoded.shape
            encoded = pd.DataFrame(np.reshape(encoded,(1,i*j)))
            try:
                self.id2pred[id] = int2maingroup[self.xgbModel.predict(encoded)[0]]
            except ValueError:
                num_features = self.xgbModel.get_booster().num_features() / j
                message = f"{id} has sequence length {i}, where is must be {int(num_features)}."
                if mode == "gui":
                    messagebox.showerror("ERROR", message)
                elif mode == "cl":
                    sys.stderr.write(f"ERROR: {message}\n")
                self.id2pred[id] = ""
                pass
        
        return

    def predict_subgroup(self,mode):

        for id in self.id2pred:
            try:
                with open(path.join(self.here,f"model/reference_{self.id2pred[id]}_group.pkl"), "rb") as file:
                    reference_dict = pkl.load(file)
            except FileNotFoundError:
                continue

            try:
                model = tf.keras.models.load_model(path.join(self.here,f"model/siamese_{self.id2pred[id]}.mdl"),compile=False)
            except OSError:
                message = f" Directory model/siamese_{self.id2pred[id]}.mdl not found."
                if mode == "gui":
                    messagebox.showerror("ERROR", message)
                elif mode == "cl":
                    sys.stderr.write(f"ERROR: {message}\n")
                continue
            encoded = np.array([letter2ohe[char] for char in self.id2seq[id]])
            score_dict = {}
            for key in reference_dict:
                score_dict[key] = []
                probe = np.broadcast_to(encoded,np.shape(reference_dict[key]))
                reference = np.array(reference_dict[key])
                score = model.predict([probe, reference])
                score_dict[key] = np.mean(score)
            self.id2pred[id] = min(score_dict, key=score_dict.get)

        return

    @staticmethod
    def create_pairs(probe, reference):
        pairs = []
        for i in range(reference.shape[0]):
            pairs.append([reference[i],probe])
        return np.array(pairs)

    def __call__(self,mode):
        self.load_xgb()
        self.predict_main_haplogroup(mode)
        self.predict_subgroup(mode)
        with open("OutputSNP_Classification.txt","w") as outputfile:
            outputfile.write("Haplogroup classification results\n")
            for id in self.id2pred:
                if self.id2pred[id] == "":
                    outputfile.write(f"{id}: not classified\n")
                else:
                    outputfile.write(f"{id}: {self.id2pred[id]}\n")
        pass
