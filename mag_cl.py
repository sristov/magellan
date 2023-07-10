# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

import pandas as pd
import numpy as np
import dask.dataframe as dd
import json
import time
import sys
from os import path
from mag_verif import check_errors, check_haplotype_conflicts, impute_haplotype
from mag_stat import mag_stat
from mag_calc import mag_calc
from mag_sampl import MagSampl
from mag_recode import recode
from mag_viz import line_print
from mag_snp import analyze_SNP, save_fasta, check_mutations, SNPClassifier

class MagellanCL():

    def __init__(self):
        return

    def read_config(self):
        
        try:
            with open("config.json", "r") as parameters_file:
                self.input_parameters = json.load(parameters_file)
        except FileNotFoundError:
            sys.stderr.write("ERROR: File config.json not found! You must run the mag_config.py first!\n")
            sys.exit()
        #for element in self.input_parameters:
        #    print(self.input_parameters[element])
        return

    def open_pedigree_file(self):
    
        try:
            
            self.data = dd.read_csv(self.input_parameters["pedigree"], dtype=np.str_).fillna('')
        except FileNotFoundError as e:
            sys.stderr.write(f'ERROR: File {self.input_parameters["pedigree"]} not found!\n')
            sys.exit()
        except PermissionError:
            sys.stderr.write(f'ERROR: You are not allowed to read {self.input_parameters["pedigree"]}.\n')
            sys.exit()
        except pd.errors.EmptyDataError:
            sys.stderr.write(f'ERROR: File {self.input_parameters["pedigree"]} is empty.\n')
            sys.exit()
        except ValueError:
            sys.stderr.write("ValueError!\n")
            sys.exit()

        if 'x' in self.data.columns.values:
            self.data = self.data.set_index('x')

        unnamed = list(filter(lambda a: 'Unnamed' in a, self.data.columns.values))

        if len(unnamed):
            self.data = self.data.drop(columns=unnamed)

        self.data=self.data.compute()
        self.data.replace({'father':{'':'0'},
                           'mother':{'':'0'},
                           'YOB':{'':'MISSING_YEAR'}},inplace=True)
        if 'available' in self.data.columns.values:
            self.data.replace({'available':{'':'0'}},inplace=True)
        self.make_maps()

        self.data = pd.DataFrame()
        self.IDlist, \
        self.FatherMap, \
        self.MotherMap, \
        self.YobMap, \
        self.GenderMap,\
        FatalError = check_errors(self.IDlist, self.FatherMap, self.MotherMap, self.YobMap, self.GenderMap, "cl")

        if FatalError:
            sys.stderr.write("ERROR: There are fatal errors in pedigree.\n")
            sys.exit()

        return

    def make_maps(self):
    
        try:
            self.IDlist = list(self.data["ID"])
            self.FatherMap = dict(zip(self.data.ID, self.data.father))
            self.MotherMap = dict(zip(self.data.ID, self.data.mother))
            self.YobMap = dict(zip(self.data.ID, self.data.YOB))
            self.GenderMap = dict(zip(self.data.ID, self.data.gender))
        except KeyError:
            sys.stderr.write('ERROR: mandatory column names are ID, father, mother, YOB and gender!\n')
            sys.stderr.write('Correct your pedigree accordingly!\n')
            sys.exit()
        except AttributeError:
            sys.stderr.write('ERROR: mandatory column names are ID, father, mother, YOB and gender!\n')
            sys.stderr.write('Correct your pedigree accordingly!\n')
            sys.exit()
        if 'haplotype' in self.data.columns.values:
            haplotyped_individuals = self.data[self.data['haplotype'] != '']
            self.HaplotypeMap = dict(zip(haplotyped_individuals.ID, haplotyped_individuals.haplotype))
            self.HaplotypedList = list(haplotyped_individuals['ID'])
            self.HaplotypeNamesList = list(set(haplotyped_individuals['haplotype']))
        else:
            self.HaplotypeMap = {}
            self.HaplotypedList = []
            self.HaplotypeNamesList = []
        if 'available' in self.data.columns.values:
            self.AvailableMap = dict(zip(self.data.ID, self.data.available))
        else:
            self.AvailableMap = {}

        return

    def make_df(self):
        self.data = pd.DataFrame(self.FatherMap.items(),columns=['ID','father'])
        self.data['mother'] = self.data['ID'].map(self.MotherMap)
        self.data['YOB'] = self.data['ID'].map(self.YobMap)
        self.data['gender'] = self.data['ID'].map(self.GenderMap)
        if len(self.HaplotypedList):
            self.data['haplotype'] = self.data['ID'].map(self.HaplotypeMap)
        if len(self.AvailableMap):
            self.data['available'] = self.data['ID'].map(self.AvailableMap)
        return

    def run_mag_verif(self):

        sys.stdout.write("Entering Magellan verification module...\n")

        if self.input_parameters["verif_options"]["lineage"] == 'paternal':
            check_haplotype_conflicts(self.IDlist,
                                      self.FatherMap,
                                      self.HaplotypeMap,
                                      self.HaplotypedList,
                                      self.HaplotypeNamesList,
                                      self.input_parameters["verif_options"]["verif_object"],
                                      "cl")
        elif self.input_parameters["verif_options"]["lineage"] == 'maternal':
            check_haplotype_conflicts(self.IDlist,
                                      self.MotherMap,
                                      self.HaplotypeMap,
                                      self.HaplotypedList,
                                      self.HaplotypeNamesList,
                                      self.input_parameters["verif_options"]["verif_object"],
                                      "cl")

        sys.stdout.write("Exiting Magellan verification module...\n")
        return
          
    def run_impute(self):

        sys.stdout.write("Entering Magellan imputation module...\n")
        if not len(self.HaplotypedList):
            sys.stderr.write("ERROR: No haplotypes present in the pedigree!\n")
            return


        ImputedHaplotypeMap,\
        VerifiedList = impute_haplotype(self.IDlist,
                                          self.input_parameters["imput_options"]["lineage"],
                                          self.GenderMap,
                                          self.MotherMap,
                                          self.HaplotypeMap,
                                          self.HaplotypedList,
                                          self.input_parameters["imput_options"]["reliability"],
                                          "cl")

        self.make_df()
        self.data['haplotype'] = self.data['ID'].map(ImputedHaplotypeMap)
        self.data['sequenced'] = self.data['ID'].isin(self.HaplotypedList).astype(int)
        if self.input_parameters["imput_options"]["reliability"] == 'low':
            self.data['verified'] = self.data['ID'].isin(VerifiedList).astype(int)
        try:
            self.data.to_csv(self.input_parameters["imput_options"]["imput_filename"], index=False)
        except AttributeError:
            pass
        except NameError:
            sys.stderr.write("ERROR: no imputed data generated!\n")
            sys.exit()

        sys.stdout.write("Exiting Magellan imputation module...\n")
        return

    def run_mag_stat(self):

        sys.stdout.write("Entering Magellan statistics module...\n")
        if self.input_parameters["stat_options"]["lineage"] == 'paternal':
            mag_stat(self.IDlist,
                     self.FatherMap,
                     self.YobMap,
                     self.GenderMap,
                     self.HaplotypeMap,
                     self.HaplotypedList,
                     self.input_parameters["stat_options"]["first_ref_year"],
                     self.input_parameters["stat_options"]["last_ref_year"],
                     self.input_parameters["stat_options"]["lineage"],
                     "cl")
        elif self.input_parameters["stat_options"]["lineage"] == 'maternal':
            mag_stat(self.IDlist,
                     self.MotherMap,
                     self.YobMap,
                     self.GenderMap,
                     self.HaplotypeMap,
                     self.HaplotypedList,
                     self.input_parameters["stat_options"]["first_ref_year"],
                     self.input_parameters["stat_options"]["last_ref_year"],
                     self.input_parameters["stat_options"]["lineage"],
                     "cl")
        sys.stdout.write("Exiting Magellan statistics module...\n")
        return

    def run_mag_calc(self):

        sys.stdout.write("Entering Magellan calculation module...\n")
        mag_calc(self.IDlist,
                 self.FatherMap,
                 self.MotherMap,
                 self.YobMap,
                 self.GenderMap,
                 self.HaplotypeMap,
                 self.HaplotypedList,
                 self.HaplotypeNamesList,
                 self.input_parameters["calc_options"]["first_ref_year"],
                 self.input_parameters["calc_options"]["last_ref_year"],
                 self.input_parameters["calc_options"]["lineage"],
                 "cl")
        sys.stdout.write("Exiting Magellan calculation module...\n")
        return

    def run_mag_sampl(self):

        sys.stdout.write("Entering Magellan sampling module...\n")
        if self.input_parameters["sampl_options"]["lineage"] == 'paternal':
            mag_sampler = MagSampl(self.IDlist,
                                   self.FatherMap,
                                   self.YobMap,
                                   self.GenderMap,
                                   self.HaplotypeMap,
                                   self.HaplotypedList,
                                   self.HaplotypeNamesList,
                                   self.AvailableMap,
                                   self.input_parameters["sampl_options"]["lineage"])

        elif self.input_parameters["sampl_options"]["lineage"] == 'maternal':
            mag_sampler = MagSampl(self.IDlist,
                                   self.MotherMap,
                                   self.YobMap,
                                   self.GenderMap,
                                   self.HaplotypeMap,
                                   self.HaplotypedList,
                                   self.HaplotypeNamesList,
                                   self.AvailableMap,
                                   self.input_parameters["sampl_options"]["lineage"])

        mag_sampler(self.input_parameters["sampl_options"]["first_ref_year"],
                    self.input_parameters["sampl_options"]["last_ref_year"],
                    self.input_parameters["sampl_options"]["sampling_method"],
                    self.input_parameters["sampl_options"]["K"],
                    "cl")
        sys.stdout.write("Exiting Magellan sampling module...\n")
        return


    def run_recode(self):

        sys.stdout.write("Entering Magellan recoding module...\n")
        self.make_df()
        self.recoded_data = recode(self.data)
        self.data = pd.DataFrame()
        try:
            self.recoded_data.to_csv(self.input_parameters["recode_options"]["recoded_filename"], index=False)
        except AttributeError:
            pass
        except NameError:
            sys.stderr.write("ERROR: no recoded data generated!\n")
            sys.exit()
        self.recoded_data = pd.DataFrame()
        sys.stdout.write("Exiting Magellan recoding module...\n")
        return


    def run_line_print(self):

        sys.stdout.write("Entering Magellan visualization module...\n")
        try:
            data = dd.read_csv(self.input_parameters["visualize_options"]["recoded_filename"], dtype=np.str_).fillna('')
        except FileNotFoundError:
            sys.stderr.write(f'ERROR: File {self.input_parameters["visualize_options"]["recoded_filename"]} not found!\n')
            sys.exit()
        except PermissionError:
            sys.stderr.write(f'ERROR: You are not allowed to read {self.input_parameters["visualize_options"]["recoded_filename"]}.\n')
            sys.exit()
        except pd.errors.EmptyDataError:
            sys.stderr.write(f'ERROR: File {self.input_parameters["visualize_options"]["recoded_filename"]} is empty.\n')
            sys.exit()
        except ValueError:
            sys.stderr.write("ValueError!\n")
            sys.exit()
        except IsADirectoryError:
            sys.stderr.write(f'ERROR: {self.input_parameters["visualize_options"]["recoded_filename"]} is a directory.\n')
            sys.exit()


        data = data.compute()
        lineage = line_print(data,
                             self.input_parameters["visualize_options"]["unit_id"],
                             self.input_parameters["visualize_options"]["start_yob"],
                             self.input_parameters["visualize_options"]["end_yob"],
                             self.input_parameters["visualize_options"]["gen_before"],
                             self.input_parameters["visualize_options"]["gen_after"],
                             self.input_parameters["visualize_options"]["lineage"],
                             "cl",
                             self.input_parameters["visualize_options"]["imagename"])

        data = pd.DataFrame()
        if len(self.input_parameters["visualize_options"]["lineage_filename"]):
            try:
                if path.splitext(self.input_parameters["visualize_options"]["lineage_filename"])[1] == '.csv':
                    lineage.to_csv(self.input_parameters["visualize_options"]["lineage_filename"], index=False)
                elif path.splitext(self.input_parameters["visualize_options"]["lineage_filename"])[1] == '.tsv':
                    lineage.to_csv(self.input_parameters["visualize_options"]["lineage_filename"], index=False, sep="\t")
                else:
                    lineage.to_string(self.input_parameters["visualize_options"]["lineage_filename"], index=False)
            except AttributeError:
                sys.stderr.write("ERROR: no lineage data present, file not created!\n")
                sys.exit()
            except NameError:
                sys.stderr.write("ERROR: no lineage generated!\n")
                sys.exit()
        sys.stdout.write("Exiting Magellan visualization module...\n")
        return
    
    def run_snp(self):

        sys.stdout.write("Entering Magellan SNP module...\n")
        try:
            cols_to_use = ['SNP Name', 'Sample ID', 'Allele1 - Forward', 'Allele2 - Forward', 'GC Score']
            final_table = dd.read_table(self.input_parameters["snp_options"]["report_filename"], 
                                        sep='\t', skiprows=9, header=0, usecols=cols_to_use,dtype=str)

        except FileNotFoundError:
            sys.stderr.write(f'ERROR: File {self.input_parameters["snp_options"]["report_filename"]} not found.\n')
            sys.exit()
        except PermissionError:
            sys.stderr.write(f'ERROR: You are not allowed to read {self.input_parameters["snp_options"]["report_filename"]}.\n')
            sys.exit()
        except pd.errors.EmptyDataError:
            sys.stderr.write(f'ERROR: File {self.input_parameters["snp_options"]["report_filename"]} is empty.\n')
            sys.exit()
        except ValueError:
            sys.exit()
        
        try:
            cols_to_use = ['Name', 'Chromosome', 'Position', 'GenTrain Score']
            snp_map = dd.read_table(self.input_parameters["snp_options"]["map_filename"], sep='\t', header=0, usecols=cols_to_use,dtype=str)

        except FileNotFoundError:
            sys.stderr.write(f'Error in opening file {self.input_parameters["snp_options"]["map_filename"]}.\n')
            sys.exit()
        except PermissionError:
            sys.stderr.write(f'ERROR: You are not allowed to read {self.input_parameters["snp_options"]["map_filename"]}.\n')
            sys.exit()
        except pd.errors.EmptyDataError:
            sys.stderr.write(f'ERROR: File {self.input_parameters["snp_options"]["map_filename"]} is empty.\n')
            sys.exit()
        except ValueError:
            sys.exit()

        try:
            SNPs_310 = pd.read_csv(self.input_parameters["snp_options"]["snp_filename"],usecols=["SNP_Name"],dtype=str)

        except FileNotFoundError:
            sys.stderr.write(f'Error in opening file {self.input_parameters["snp_options"]["snp_filename"]}.\n')
            sys.exit()
        except PermissionError:
            sys.stderr.write(f'ERROR: You are not allowed to read {self.input_parameters["snp_options"]["snp_filename"]}.\n')
            sys.exit()
        except pd.errors.EmptyDataError:
            sys.stderr.write(f'ERROR: File {self.input_parameters["snp_options"]["snp_filename"]} is empty.\n')
            sys.exit()
        except ValueError:
            sys.stderr.write(f'"ERROR: File {self.input_parameters["snp_options"]["snp_filename"]} must contain header SNP_Name.\n')
            sys.exit()

        try:
            cols_to_use = ['Sample_ID', 'Population']
            samples_population = pd.read_table(self.input_parameters["snp_options"]["population_filename"], sep='\t', header=0, usecols=cols_to_use,dtype=str)

        except FileNotFoundError:
            sys.stderr.write(f'Error in opening file {self.input_parameters["snp_options"]["population_filename"]}.\n')
            sys.exit()
        except PermissionError:
            sys.stderr.write(f'ERROR: You are not allowed to read {self.input_parameters["snp_options"]["population_filename"]}.\n')
            sys.exit()
        except pd.errors.EmptyDataError:
            sys.stderr.write(f'ERROR: File {self.input_parameters["snp_options"]["population_filename"]} is empty.\n')
            sys.exit()
        except ValueError:
            sys.stderr.write(f'ERROR: ValueError in reading {self.input_parameters["snp_options"]["population_filename"]}.\n')
            sys.exit()
        except AssertionError:
            sys.stderr.write(f'ERROR: AssertionError in reading {self.input_parameters["snp_options"]["population_filename"]}.\n')
            sys.exit()
    
        taurus_subset, proba = analyze_SNP(final_table,
                                           snp_map,
                                           SNPs_310,
                                           samples_population,
                                           self.input_parameters["snp_options"]["gcscore_thresold"],
                                           self.input_parameters["snp_options"]["gtscore_thresold"],
                                           self.input_parameters["snp_options"]["snp_position"],
                                           "cl")
        
        if len(self.input_parameters["snp_options"]["fasta_filename"]):
            save_fasta(self.input_parameters["snp_options"]["fasta_filename"], proba)
        if self.input_parameters["snp_options"]["mut_check"] == True:
            check_mutations(taurus_subset,"cl")
        sys.stdout.write("Exiting Magellan SNP module...\n")
        return

    def read_fasta(self):
        id2seq = {}
        try:
            with open(self.input_parameters["classify_options"]["fasta_filename"], 'r') as fasta_file:
                for line in fasta_file:
                    if line[0] == ">":
                        sample_id = line[1:].strip()
                    else:
                        snp_sequence = line.strip()
                        id2seq[sample_id] = snp_sequence

        except FileNotFoundError:
            if self.input_parameters["snp_options"]["fasta_filename"] == "":
                return id2seq
            sys.stderr.write(f"Error in opening file {self.input_parameters['classify_options']['fasta_filename']}.\n")
            return id2seq
        except PermissionError:
            sys.stderr.write(f"You are not allowed to read {self.input_parameters['classify_options']['fasta_filename']}.\n")
            return id2seq
        
        return id2seq

    def run_classification(self):
        sys.stdout.write("Entering Magellan classification module...\n")
        if self.input_parameters["classify_options"]["sequence_source"] == "ped":
            if len(self.HaplotypeMap) == 0:
                sys.stderr.write("ERROR: No haplotypes found in the pedigree!\n")
                return
            classifier = SNPClassifier("xgb_main.json",self.HaplotypeMap)
            classifier("cl")
        elif self.input_parameters["classify_options"]["sequence_source"] == "fasta":
            id2seq = self.read_fasta()
            if len(id2seq) == 0:
                return
            classifier = SNPClassifier("xgb_main.json",id2seq)
            classifier("cl")

        if self.input_parameters["classify_options"]["impute_haplogroup"]:
            try:
                self.make_df()
            except AttributeError:
                sys.stderr.write("ERROR: pedigree for inserting predictions is not loaded!\n")
                return
            try:
                self.data["haplotype"] = self.data["ID"].map(classifier.id2pred)
            except AttributeError:
                sys.stderr.write("ERROR: classification not performed.\n")
                return

        if len(self.input_parameters["classify_options"]["csv_name"]):
            try:
                self.data.to_csv(self.input_parameters["classify_options"]["csv_name"],index=False)
                self.data = pd.DataFrame()
            except AttributeError:
                sys.stderr.write(f"ERROR: AttributeError in saving {self.input_parameters['classify_options']['csv_name']}.\n")
                pass
            except NameError:
                sys.stderr.write("ERROR: no data generated, load them and impute them!")
            except PermissionError:
                sys.stderr.write(f"You are not allowed to write {self.input_parameters['classify_options']['csv_name']}.")
        sys.stdout.write("Exiting Magellan classification module...\n")
        return

    def __call__(self):
        sys.stdout.write("Welcome to Maternal Genealogy Lineage Analyser v2.0\n")
        sys.stdout.write("Authors: V. Brajkovic, I. Curik, D. Hrsak, S. Ristov\n")
        sys.stdout.write("MaGelLAn v2.0 is licensed under GNU General Public License v3.0\n")
        self.read_config()
        start_out = time.time()
        if len(self.input_parameters["pedigree"]):
            start_in = time.time()
            self.open_pedigree_file()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent reading pedigree file: {elapsed:12.6f}\n")

        if self.input_parameters["mag_verif"] == True:
            start_in = time.time()
            self.run_mag_verif()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the verification module: {elapsed:12.6f}\n")
        if self.input_parameters["mag_imput"] == True:
            start_in = time.time()
            self.run_impute()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the imputation submodule: {elapsed:12.6f}\n")
        if self.input_parameters["mag_stat"] == True:
            start_in = time.time()
            self.run_mag_stat()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the statistics module: {elapsed:12.6f}\n")
        if self.input_parameters["mag_calc"] == True:
            start_in = time.time()
            self.run_mag_calc()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the calculation module: {elapsed:12.6f}\n")
        if self.input_parameters["mag_sampl"] == True:
            start_in = time.time()
            self.run_mag_sampl()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the sampling module: {elapsed:12.6f}\n")
        if self.input_parameters["mag_recode"] == True:
            start_in = time.time()
            self.run_recode()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the recoding module: {elapsed:12.6f}\n")
        if self.input_parameters["mag_visualize"] == True:
            start_in = time.time()
            self.run_line_print()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the visualization module: {elapsed:12.6f}\n")
        if self.input_parameters["mag_snp"] == True:
            start_in = time.time()
            self.run_snp()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the SNP module: {elapsed:12.6f}\n")
        if self.input_parameters["mag_classify"] == True:
            start_in = time.time()
            self.run_classification()
            end_in = time.time()
            elapsed = end_in - start_in
            sys.stdout.write(f"Time spent in the classification submodule: {elapsed:12.6f}\n")
 
        end_out = time.time()
        elapsed = end_out - start_out
        sys.stdout.write(f"Total time spent in Magellan: {elapsed:12.6f}\n")
        sys.stdout.write("Exiting Maternal Genealogy Lineage Analyser v2.0\n")

        return
