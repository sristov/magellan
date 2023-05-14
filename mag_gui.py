# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

from os import path
import pandas as pd
import numpy as np
import dask.dataframe as dd
import tkinter as tk
#from tk import scrolledtext
from tkinter import filedialog as fd
from tkinter import ttk
from tkinter.messagebox import showerror, showwarning
from pandastable import Table
#from pedigree import Pedigree
from mag_verif import check_errors, check_haplotype_conflicts, impute_haplotype
from mag_stat import mag_stat
from mag_calc import mag_calc
from mag_sampl import MagSampl
from mag_recode import recode
from mag_viz import line_print
from mag_snp import analyze_SNP, save_fasta, check_mutations, SNPClassifier
#from PIL import ImageTk, Image
#from tkPDFViewer import tkPDFViewer

class MagellanGUI(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("Magellan-v2.0demo")
        self.geometry("750x750")

        self.data = pd.DataFrame()
        self.founder = ''
        self.start_year = 0
        self.end_year = 0
        self.Lineage = tk.StringVar(None, 'maternal')
        self.VerifObject = tk.StringVar(None, 'haplogroup')
        self.ImputReliability = tk.StringVar(None, 'low')
        self.FirstRefYear = 0
        self.LastRefYear = 0
        self.K = 100
        self.generations_before = 0
        self.generations_after = 0
        self.gc_score_thresh = 0.0
        self.gt_score_thresh = 0.0
        self.snp_position = tk.StringVar(None,'MT')
        self.sequence_source = tk.StringVar(None,'PED')

        label1 = tk.Label(self,
            text="Welcome to Maternal Genealogy Lineage Analyser v2.0", 
            font=('Aerial 18 bold'),
            foreground='green')
        label1.grid(row=0, column=0, sticky=tk.S, pady=5)
        label2 = tk.Label(self, 
            text="Authors: V. Brajkovic, I. Curik, D. Hrsak, S. Ristov", 
            font=('Aerial 14 bold'))
        label2.grid(row=1, column=0, sticky=tk.S, pady=5)
        label3 = tk.Label(self, 
            text="MaGelLAn v2.0 is licenced under GNU General Public License v3.0", 
            font=('Aerial 12'))
        label3.grid(row=2, column=0, sticky=tk.S, pady=5)
        Buttonbar_main = ttk.Frame(master=self)
        Buttonbar_main.grid(row=3, column=0, sticky=tk.S, pady=5)
    
        btn_open_pdg = tk.ttk.Button(master=Buttonbar_main, 
                                     text="Open pedigree file", 
                                     command=self.open_csv_file)
        btn_open_pdg.grid(row=0, column=0, sticky=tk.E)

        btn_display_pdg = tk.ttk.Button(master=Buttonbar_main,
                                             text="Display pedigree file", 
                                             command=self.display_dataframe)
        btn_display_pdg.grid(row=0, column=1, sticky=tk.E)

        btn_close = tk.ttk.Button(master=Buttonbar_main, 
                                  text='Close program', 
                                  command=self.destroy)

        btn_close.grid(row=0, column=2, sticky=tk.E)

        tabControl = ttk.Notebook(self)
  
        tab1 = ttk.Frame(tabControl)
        tab2 = ttk.Frame(tabControl)
        tab3 = ttk.Frame(tabControl)
        tab4 = ttk.Frame(tabControl)
        tab5 = ttk.Frame(tabControl)
        tab6 = ttk.Frame(tabControl)
        tab7 = ttk.Frame(tabControl)
  
        tabControl.add(tab1, text ='mag_verif')
        tabControl.add(tab2, text ='mag_stat')
        tabControl.add(tab3, text ='mag_calc')
        tabControl.add(tab4, text ='mag_sampl')
        tabControl.add(tab5, text ='recode')
        tabControl.add(tab6, text ='visualize')
        tabControl.add(tab7, text ='SNP analysis')
        tabControl.grid(row=4, column=0, sticky=tk.E, pady=5)
  
        ttk.Label(tab1, 
            text ="Magellan verification module",
            font=('Aerial 10 bold')).grid(column = 0, row = 0, padx = 5, pady = 5)

        Verifbar = ttk.Frame(master=tab1)
        Verifbar.grid(column = 0,row = 1,padx = 5, pady = 5)
        tk.Label(master=Verifbar, text="Lineage for conflict check:").grid(row=0, column=0, sticky=tk.E)
        Radiobar_stat = ttk.Frame(master=Verifbar)
        Radiobar_stat.grid(row=0, column=1, sticky=tk.E)
        tk.Radiobutton(master=Radiobar_stat,
                       text="maternal",
                       variable=self.Lineage,
                       value="maternal").grid(row=0, column=0, sticky=tk.E)

        tk.Radiobutton(master=Radiobar_stat,
                       text="paternal",
                       variable=self.Lineage,
                       value="paternal").grid(row=0, column=1, sticky=tk.E)
        Verifbar1 = ttk.Frame(master=tab1)
        Verifbar1.grid(column = 0,row = 2,padx = 5, pady = 5)
        tk.Label(master=Verifbar1, text="Object of verification:").grid(row=0, column=0, sticky=tk.E)
        Radiobar_verif1 = ttk.Frame(master=Verifbar1)
        Radiobar_verif1.grid(row=0, column=1, sticky=tk.E)
        tk.Radiobutton(master=Radiobar_verif1,
                       text="Haplogroup",
                       variable=self.VerifObject,
                       value="haplogroup").grid(row=0, column=0, sticky=tk.E)

        tk.Radiobutton(master=Radiobar_verif1,
                       text="SNP Sequence",
                       variable=self.VerifObject,
                       value="snpseq").grid(row=0, column=1, sticky=tk.E)
        Buttonbar_verif = ttk.Frame(master=tab1)
        Buttonbar_verif.grid(column=0, row=3,padx=0, pady=0)

        btn_run_verif = tk.ttk.Button(Buttonbar_verif,text="Run mag_verif",
                                      command=self.run_mag_verif)
        btn_run_verif.grid(column = 0, row = 0, padx = 0, pady = 0)

        btn_verif_result = tk.ttk.Button(Buttonbar_verif, text="View report",
                                        command=self.open_verif_report)
        btn_verif_result.grid(column=1, row=0, padx=0, pady=0)

        btn_verif_del = tk.ttk.Button(Buttonbar_verif, text="Clear report",
                                        command=self.delete_verif_report)
        btn_verif_del.grid(column=2, row=0, padx=0, pady=0)

        Buttonbar_verif1 = ttk.Frame(master=tab1)
        Buttonbar_verif1.grid(column=0, row=4, padx=0, pady=0)
        btn_verif_imp = tk.ttk.Button(Buttonbar_verif1, text="Impute haplotype",
                                      command=self.run_impute)
        btn_verif_imp.grid(column=0, row=0, padx=0, pady=0)
        tk.Label(master=Buttonbar_verif1, text="Reliability:").grid(row=0, column=1, sticky=tk.E)
        tk.Radiobutton(master=Buttonbar_verif1, 
                       text="high", 
                       variable=self.ImputReliability,
                       value="high").grid(row=0,column=2, sticky=tk.E)

        tk.Radiobutton(master=Buttonbar_verif1, 
                       text="low", 
                       variable=self.ImputReliability,
                       value="low").grid(row=0,column=3, sticky=tk.E)

        btn_save_imp = tk.ttk.Button(Buttonbar_verif1,
                                     text="Save imputed CSV as",
                                     command=self.save_imputed_csv)
        btn_save_imp.grid(column = 4, row = 0, padx = 0, pady = 0)

        self.txtarea_verif = tk.scrolledtext.ScrolledText(tab1, width=80, height=40)
        self.txtarea_verif.grid(column=0,row=5,padx=0,pady=0)
        self.txtarea_verif.configure(font=(24))

        ttk.Label(master=tab2,
                  text ="Magellan statistics module",
                  font=('Aerial 10 bold')).grid(column = 0, row = 0, padx = 5, pady = 5)
        Statbar = ttk.Frame(master=tab2)
        Statbar.grid(column = 0,row = 1,padx = 0, pady = 0)
        lbl_first_refyear_stat = tk.Label(master=Statbar, text="First reference year:")
        lbl_first_refyear_stat.grid(row=0, column=0, sticky=tk.E)
        ent_first_refyear_stat = tk.Entry(master=Statbar, width=20)
        ent_first_refyear_stat.insert(tk.END, '0')
        ent_first_refyear_stat.grid(row=0, column=1)

        def set_fry_stat():
            try:
                self.FirstRefYear = int(ent_first_refyear_stat.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid integer!")
        tk.Button(master=Statbar, text="Set", command=set_fry_stat).grid(column = 2, row = 0, padx = 0, pady = 0)

        Statbar1 = ttk.Frame(master=tab2)
        Statbar1.grid(column = 0, row = 2, padx = 0, pady = 0)
        lbl_last_refyear_stat = tk.Label(master=Statbar1, text="Last reference year:")
        lbl_last_refyear_stat.grid(row=0, column=0, sticky=tk.E)
        ent_last_refyear_stat = tk.Entry(master=Statbar1, width=20)
        ent_last_refyear_stat.insert(tk.END, '0')
        ent_last_refyear_stat.grid(row=0, column=1)
        
        def set_lry_stat():
            try:
                self.LastRefYear = int(ent_last_refyear_stat.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid integer!")
        tk.Button(master=Statbar1, text="Set", command=set_lry_stat).grid(column = 2, row = 0, padx = 0, pady = 0)
        
        Radiobar_stat = ttk.Frame(master=tab2)
        Radiobar_stat.grid(row=3,column=0, sticky=tk.S)
        tk.Radiobutton(master=Radiobar_stat, 
                       text="maternal", 
                       variable=self.Lineage,
                       value="maternal").grid(row=0,column=0, sticky=tk.E)

        tk.Radiobutton(master=Radiobar_stat, 
                       text="paternal", 
                       variable=self.Lineage,
                       value="paternal").grid(row=0,column=1, sticky=tk.E)

        Buttonbar_stat = ttk.Frame(master=tab2)
        Buttonbar_stat.grid(column=0, row=4, padx=0, pady=0)
        btn_run_mag_stat = tk.ttk.Button(Buttonbar_stat, text="Run mag_stat",
                                         command=self.run_mag_stat)
        btn_run_mag_stat.grid(column = 0, row = 0, padx = 0, pady = 0)
        btn_stat_result = tk.ttk.Button(Buttonbar_stat, text="View report",
                                        command=self.open_stat_report)
        btn_stat_result.grid(column=1,row=0,padx=0,pady=0)

        btn_stat_clear = tk.ttk.Button(Buttonbar_stat, text="Clear report",
                                       command=self.delete_stat_report)
        btn_stat_clear.grid(column=2,row=0,padx=0,pady=0)

        self.txtarea_stat = tk.scrolledtext.ScrolledText(master=tab2, width=80, height=40)
        self.txtarea_stat.grid(column=0, row=5, padx=0, pady=0)
        self.txtarea_stat.configure(font=(24))

        ttk.Label(tab3,
            text ="Magellan calculation module",
            font=('Aerial 10 bold')).grid(column = 0, row = 0, padx = 0, pady = 5, sticky=tk.E)
        Calcbar = ttk.Frame(master=tab3)
        Calcbar.grid(column = 0,row = 1,padx = 0, pady = 0, sticky=tk.E)
        lbl_first_refyear_calc = tk.Label(master=Calcbar, text="First reference year:")
        ent_first_refyear_calc = tk.Entry(master=Calcbar, width=20)
        ent_first_refyear_calc.insert(tk.END, '0')
        lbl_first_refyear_calc.grid(row=0, column=0, sticky=tk.E)
        ent_first_refyear_calc.grid(row=0, column=1)
        def set_fry_calc():
            try:
                self.FirstRefYear = int(ent_first_refyear_calc.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid integer!")
        tk.Button(Calcbar, text="Set", command=set_fry_calc).grid(column = 2, row = 0, padx = 0, pady = 0)
        
        Calcbar1 = ttk.Frame(master=tab3)
        Calcbar1.grid(column = 0,row = 2,padx = 0, pady = 0, sticky=tk.E)

        lbl_last_refyear_calc = tk.Label(master=Calcbar1, text="Last reference year:")
        lbl_last_refyear_calc.grid(row=0, column=0, sticky=tk.E)
        ent_last_refyear_calc = tk.Entry(master=Calcbar1, width=20)
        ent_last_refyear_calc.insert(tk.END, '0')
        ent_last_refyear_calc.grid(row=0, column=1)
        def set_lry_calc():
            try:
                self.LastRefYear = int(ent_last_refyear_calc.get())
            except ValueError:
                showerror("ERROR", "ERROR: the value you have entered is not a valid integer!")
        tk.Button(Calcbar1, text="Set", command=set_lry_calc).grid(column = 2, row = 0, padx = 0, pady = 0)


        Radiobar_calc = ttk.Frame(master=tab3)
        Radiobar_calc.grid(row=3, column=0, sticky=tk.E)
        tk.Radiobutton(master=Radiobar_calc,
                       text="maternal",
                       variable=self.Lineage,
                       value="maternal").grid(row=0, column=0, sticky=tk.E)

        tk.Radiobutton(master=Radiobar_calc,
                       text="paternal",
                       variable=self.Lineage,
                       value="paternal").grid(row=0, column=1, sticky=tk.E)

        Buttonbar_calc = ttk.Frame(master=tab3)
        Buttonbar_calc.grid(column=0, row=4, padx=0, pady=0, sticky=tk.E)

        btn_run_mag_calc = tk.ttk.Button(Buttonbar_calc,
                                         text="Run mag_calc",
                                         command=self.run_mag_calc)
        btn_run_mag_calc.grid(column = 0, row = 0, padx = 0, pady = 5, sticky = tk.E)

        btn_calc_result = tk.ttk.Button(Buttonbar_calc, text="View report",
                                        command=self.open_calc_report)
        btn_calc_result.grid(column=1, row=0, padx=0, pady=0)

        btn_calc_clear = tk.ttk.Button(Buttonbar_calc, text="Clear report",
                                        command=self.delete_calc_report)
        btn_calc_clear.grid(column=2, row=0, padx=0, pady=0)

        self.txtarea_calc = tk.scrolledtext.ScrolledText(tab3, width=80, height=40)
        self.txtarea_calc.grid(column=0, row=5, columnspan=3, padx=0, pady=0)
        self.txtarea_calc.configure(font=(24))

        ttk.Label(tab4,
            text ="Magellan sampling module",
            font=('Aerial 10 bold')).grid(column = 0, row = 0, padx = 5, pady = 5, sticky=tk.E)

        Samplbar = ttk.Frame(master=tab4)
        Samplbar.grid(column = 0,row = 1,padx = 0, pady = 0, sticky=tk.E)
        lbl_K = tk.Label(master=Samplbar, text="How many to sequence:")
        lbl_K.grid(row=0, column=0, sticky=tk.E)
        ent_K = tk.Entry(master=Samplbar, width=20)
        ent_K.insert(tk.END, '100')
        ent_K.grid(row=0, column=1)
        def set_hmts_sampl():
            try:
                self.K = int(ent_K.get())
            except ValueError:
                showerror("ERROR", "ERROR: the value you have entered is not a valid integer!")
        tk.Button(master=Samplbar, text="Set", command=set_hmts_sampl).grid(column = 2, row = 0, padx = 0, pady = 0)
        
        
        Samplbar1 = ttk.Frame(master=tab4)
        Samplbar1.grid(column = 0,row = 2,padx = 0, pady = 0, sticky=tk.E)
        lbl_first_refyear_sampl = tk.Label(master=Samplbar1, text="First reference year:")
        lbl_first_refyear_sampl.grid(row=0, column=0, sticky=tk.E)
        ent_first_refyear_sampl = tk.Entry(master=Samplbar1, width=20)
        ent_first_refyear_sampl.insert(tk.END, '0')
        ent_first_refyear_sampl.grid(row=0, column=1)
        def set_fry_sampl():
            try:
                self.FirstRefYear = int(ent_first_refyear_sampl.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid integer!")
        tk.Button(master=Samplbar1, text="Set", command=set_fry_sampl).grid(column = 2, row = 0, padx = 0, pady = 0)

        Samplbar2 = ttk.Frame(master=tab4)
        Samplbar2.grid(column = 0,row = 3,padx = 0, pady = 0, sticky=tk.E)
        lbl_last_refyear_sampl = tk.Label(master=Samplbar2, text="Last reference year:")
        lbl_last_refyear_sampl.grid(row=0, column=0, sticky=tk.E)
        ent_last_refyear_sampl = tk.Entry(master=Samplbar2, width=20)
        ent_last_refyear_sampl.insert(tk.END, '0')
        ent_last_refyear_sampl.grid(row=0, column=1)
        def set_lry_sampl():
            try:
                self.LastRefYear = int(ent_last_refyear_sampl.get())
            except ValueError:
                showerror("ERROR", "ERROR: the value you have entered is not a valid integer!")
        tk.Button(master=Samplbar2, text="Set", command=set_lry_sampl).grid(column = 2, row = 0, padx = 0, pady = 0)

        Samplbar3 = ttk.Frame(master=tab4)
        Samplbar3.grid(column = 0,row = 4,padx = 0, pady = 5, sticky=tk.E)
        tk.Label(master=Samplbar3, text="Lineage to sample:").grid(row=0, column=0, sticky=tk.E)
        Radiobar1 = ttk.Frame(master=Samplbar3)
        Radiobar1.grid(row=0,column=1, sticky=tk.E)
        tk.Radiobutton(master=Radiobar1, 
                       text="maternal", 
                       variable=self.Lineage,
                       value="maternal").grid(row=0,column=0, sticky=tk.E)

        tk.Radiobutton(master=Radiobar1, 
                       text="paternal", 
                       variable=self.Lineage,
                       value="paternal").grid(row=0,column=1, sticky=tk.E)
        
        Samplbar4 = ttk.Frame(master=tab4)
        Samplbar4.grid(column = 0,row = 5,padx = 0, pady = 0, sticky=tk.E)
        tk.Label(master=Samplbar4, text="Sampling method:").grid(row=0, column=0, sticky=tk.E)
        self.SamplingMethod = tk.StringVar(None,'greedy')
        Radiobar2 = ttk.Frame(master=Samplbar4)
        Radiobar2.grid(row=0,column=1, sticky=tk.E)
        tk.Radiobutton(master=Radiobar2, 
                       text="greedy", 
                       variable=self.SamplingMethod,
                       value="greedy").grid(row=0,column=0, sticky=tk.E)

        tk.Radiobutton(master=Radiobar2, 
                       text="optimal", 
                       variable=self.SamplingMethod,
                       value="optimal").grid(row=0,column=1, sticky=tk.E)
        
        Buttonbar_sampl = ttk.Frame(master=tab4)
        Buttonbar_sampl.grid(column=0, row=6, padx=0, pady=0, sticky=tk.E)

        btn_run_sampl = tk.ttk.Button(Buttonbar_sampl,
                                      text="Run mag_sampl",
                                      command=self.run_mag_sampl)
        btn_run_sampl.grid(column = 0, row = 0, padx = 0, pady = 0)

        btn_sampl_result = tk.ttk.Button(Buttonbar_sampl, text="View report",
                                        command=self.open_sampl_report)
        btn_sampl_result.grid(column=1, row=0, padx=0, pady=0)

        btn_sampl_clear = tk.ttk.Button(Buttonbar_sampl, text="Clear report",
                                        command=self.delete_sampl_report)
        btn_sampl_clear.grid(column=2, row=0, padx=0, pady=0)

        self.txtarea_sampl = tk.scrolledtext.ScrolledText(tab4, width=80, height=40)
        self.txtarea_sampl.grid(column=0, row=7, columnspan=3, padx=0, pady=0)
        self.txtarea_sampl.configure(font=(24))

        ttk.Label(tab5,
            text ="Magellan auxiliary module recode",
            font=('Aerial 10 bold')).grid(column = 0, row = 0, padx = 5, pady = 5)
        Buttonbar_recode = ttk.Frame(master=tab5)
        Buttonbar_recode.grid(column=0, row=1, columnspan=3, padx=0, pady=0)

        btn_run_recode = tk.ttk.Button(Buttonbar_recode,
                                       text="Run recode",
                                       command=self.run_recode)
        btn_run_recode.grid(column = 0, row = 0, padx = 0, pady = 0)
        btn_save_recode = tk.ttk.Button(Buttonbar_recode,
                                        text="Save recoded CSV as",
                                        command=self.save_recoded_csv)
        btn_save_recode.grid(column = 1, row = 0, padx = 0, pady = 0)
        ttk.Label(tab6, text ="Magellan auxiliary module visualize line",
                  font=('Aerial 10 bold')).grid(column = 0, row = 0, padx = 5, pady = 5)
        lbl_founder = tk.Label(master=tab6, text="Unit ID:")
        ent_founder = tk.Entry(master=tab6, width=10)
        lbl_start_year = tk.Label(master=tab6, text="Start year of birth:")
        ent_start_year = tk.Entry(master=tab6, width=10)
        ent_start_year.insert(tk.END, '0')
        lbl_end_year = tk.Label(master=tab6, text="End year of birth:")
        ent_end_year = tk.Entry(master=tab6, width=10)
        ent_end_year.insert(tk.END, '0')
        lbl_generations_before = tk.Label(master=tab6, text="Generations before:")
        ent_generations_before = tk.Entry(master=tab6, width=10)
        ent_generations_before.insert(tk.END, '0')
        lbl_generations_after = tk.Label(master=tab6, text="Generations after:")
        ent_generations_after = tk.Entry(master=tab6, width=10)
        ent_generations_after.insert(tk.END, '0')
        lbl_founder.grid(row=1, column=0, sticky=tk.E)
        ent_founder.grid(row=1, column=1)
        lbl_start_year.grid(row=2, column=0, sticky=tk.E)
        ent_start_year.grid(row=2, column=1)
        lbl_end_year.grid(row=3, column=0, sticky=tk.E)
        ent_end_year.grid(row=3, column=1)
        lbl_generations_before.grid(row=4, column=0, sticky=tk.E)
        ent_generations_before.grid(row=4, column=1)
        lbl_generations_after.grid(row=5, column=0, sticky=tk.E)
        ent_generations_after.grid(row=5, column=1)
        tk.Label(master=tab6, text="Lineage to visualize:").grid(row=6, column=0, sticky=tk.E)
        Radiobar_viz = ttk.Frame(master=tab6)
        Radiobar_viz.grid(row=6,column=1, sticky=tk.E)
        tk.Radiobutton(master=Radiobar_viz, 
                       text="maternal", 
                       variable=self.Lineage,
                       value="maternal").grid(row=0,column=0, sticky=tk.E)

        tk.Radiobutton(master=Radiobar_viz, 
                       text="paternal", 
                       variable=self.Lineage,
                       value="paternal").grid(row=0,column=1, sticky=tk.E)
        Buttonbar_viz = ttk.Frame(master=tab6)
        Buttonbar_viz.grid(column=0, row=7, columnspan=3, padx=0, pady=0)

        btn_run_line_print = tk.ttk.Button(Buttonbar_viz,
                                       text="Run lineage print",
                                       command=self.run_line_print)
        btn_run_line_print.grid(column = 0, row = 0, padx = 0, pady = 5, sticky = tk.W+tk.E)
        #self.btn_display_pdg = tk.ttk.Button(Buttonbar_viz,
        #                                     text="Display founder line", 
        #                                     command=self.show_image)
        #self.btn_display_pdg.grid(column = 1, row = 0, padx = 0, pady = 0, sticky = tk.W+tk.E)
        btn_save_founderline = tk.ttk.Button(Buttonbar_viz,
                                             text="Save lineage as",
                                             command=self.save_lineage_csv)
        btn_save_founderline.grid(column = 1, row = 0, padx = 0, pady = 0, sticky = tk.W+tk.E)

        #self.canvas_ped = tk.Canvas(tab6,width=480,height=480,bg="white")
        #self.canvas_ped.grid(column=0, row=8, columnspan=3, padx=0, pady=0, sticky=tk.W+tk.E)
      
        ttk.Label(tab7, 
            text ="Magellan SNP analysis module",
            font=('Aerial 10 bold')).grid(column = 0, row = 0, padx = 5, pady = 5)
        SNPbar = ttk.Frame(master=tab7)
        SNPbar.grid(column = 0,row = 1,padx = 0, pady = 5, sticky=tk.S)
        lbl_set_gcscore = tk.Label(master=SNPbar, text="Set GC Score Threshold:")
        lbl_set_gcscore.grid(row=0, column=0, sticky=tk.W)
        ent_set_gcscore = tk.Entry(master=SNPbar, width=10)
        ent_set_gcscore.insert(tk.END, '0.0')
        ent_set_gcscore.grid(row=0, column=1)
        def set_gc_score():
            try:
                if float(ent_set_gcscore.get()) < 0.0 or float(ent_set_gcscore.get()) > 1.0:
                    showerror("ERROR", f"ERROR: Allowed values in range [0.0,1.0]!")
                    return
                self.gc_score_thresh = float(ent_set_gcscore.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid floating point number!")
        tk.Button(master=SNPbar, text="Set", command=set_gc_score).grid(column = 2, row = 0, padx = 0, pady = 0)
        

        SNPbar1 = ttk.Frame(master=tab7)
        SNPbar1.grid(column = 0,row = 2,padx = 0, pady = 5, sticky=tk.S)
        lbl_set_gtscore = tk.Label(master=SNPbar1, text="Set GT Score Threshold:")
        lbl_set_gtscore.grid(row=0, column=0, sticky=tk.E)
        ent_set_gtscore = tk.Entry(master=SNPbar1, width=10)
        ent_set_gtscore.insert(tk.END, '0.0')
        ent_set_gtscore.grid(row=0, column=1)
        def set_gt_score():
            try:
                if float(ent_set_gtscore.get()) < 0.0 or float(ent_set_gtscore.get()) > 1.0:
                    showerror("ERROR", f"ERROR: Allowed values in range [0.0,1.0]!")
                    return
                self.gt_score_thresh = float(ent_set_gtscore.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid floating point number!")
        tk.Button(master=SNPbar1, text="Set", command=set_gt_score).grid(column = 2, row = 0, padx = 0, pady = 0)

        SNPbar2 = ttk.Frame(master=tab7)
        SNPbar2.grid(column = 0,row = 3,padx = 0, pady = 5, sticky=tk.S)
        lbl_snp_position = tk.Label(master=SNPbar2, text="Part of genome for analysis:")
        lbl_snp_position.grid(row=0, column=0, sticky=tk.E)
        Radiobar_snp = ttk.Frame(master=SNPbar2)
        Radiobar_snp.grid(row=0,column=1, sticky=tk.E)
        
        tk.Radiobutton(master=Radiobar_snp, 
                       text="mitochondrial", 
                       variable=self.snp_position,
                       value="MT").grid(row=0,column=0, sticky=tk.E)

        tk.Radiobutton(master=Radiobar_snp, 
                       text="Y chromosome", 
                       variable=self.snp_position,
                       value="Y").grid(row=0,column=1, sticky=tk.E)
        tk.Radiobutton(master=Radiobar_snp, 
                       text="X chromosome", 
                       variable=self.snp_position,
                       value="X").grid(row=0,column=2, sticky=tk.E)

        tk.Radiobutton(master=Radiobar_snp, 
                       text="autosome", 
                       variable=self.snp_position,
                       value="AUTO").grid(row=0,column=3, sticky=tk.E)

        Buttonbar_snp = ttk.Frame(master=tab7)
        Buttonbar_snp.grid(column=0, row=4, padx=5, pady=5)

        self.btn_open_rpt = tk.ttk.Button(master=Buttonbar_snp,
                                          text="Load Report file",
                                          command=self.load_final_report)
        self.btn_open_rpt.grid(row=0, column=0, sticky=tk.E)

        self.btn_open_map = tk.ttk.Button(master=Buttonbar_snp,
                                             text="Load MAP file",
                                             command=self.load_snp_map)
        self.btn_open_map.grid(row=0, column=1, sticky=tk.E)

        self.btn_open_snp = tk.ttk.Button(master=Buttonbar_snp,
                                       text='Load SNP list',
                                       command=self.load_snp_list)
        self.btn_open_snp.grid(row=0, column=2, sticky=tk.E)
        self.btn_open_pop = tk.ttk.Button(master=Buttonbar_snp,
                                       text='Load population list',
                                       command=self.load_population)
        self.btn_open_pop.grid(row=0, column=3, sticky=tk.E)
        Buttonbar_snp1 = ttk.Frame(master=tab7)
        Buttonbar_snp1.grid(column=0, row=5, padx=0, pady=5)

        self.btn_run_snp = tk.ttk.Button(master=Buttonbar_snp1,
                                          text="Run SNP analysis",
                                          command=self.run_snp)
        self.btn_run_snp.grid(row=0, column=0, sticky=tk.E)

        self.btn_save_fasta = tk.ttk.Button(master=Buttonbar_snp1,
                                          text="Save FASTA file",
                                          command=self.write_fasta)
        self.btn_save_fasta.grid(row=0, column=1, sticky=tk.E)

        self.btn_mutations = tk.ttk.Button(master=Buttonbar_snp1,
                                          text='Check deleterious mutations',
                                          command=self.mutations_check)
        self.btn_mutations.grid(row=0, column=2, sticky=tk.E)

        Predbar_snp = ttk.Frame(master=tab7)
        Predbar_snp.grid(column=0, row=6, padx=0, pady=5)

        self.btn_classify_snp = tk.ttk.Button(master=Predbar_snp,
                                          text="Run classification",
                                          command=self.run_classify)
        self.btn_classify_snp.grid(row=0, column=0, sticky=tk.E)

        self.btn_impute_snp = tk.ttk.Button(master=Predbar_snp,
                                          text="Insert predictions",
                                          command=self.run_to_haplogroup)
        self.btn_impute_snp.grid(row=0, column=1, sticky=tk.E)
        self.btn_save_csv = tk.ttk.Button(master=Predbar_snp,
                                          text="Save pedigree",
                                          command=self.save_imputed_csv)
        self.btn_save_csv.grid(row=0, column=2, sticky=tk.E)

        tk.Radiobutton(master=Predbar_snp, 
                       text="pedigree", 
                       variable=self.sequence_source,
                       value="PED").grid(row=1,column=0, sticky=tk.E)

        tk.Radiobutton(master=Predbar_snp, 
                       text="FASTA", 
                       variable=self.sequence_source,
                       value="FASTA").grid(row=1,column=1, sticky=tk.E)

        self.txtarea_snp = tk.scrolledtext.ScrolledText(tab7, width=80, height=40)
        self.txtarea_snp.grid(column=0, row=7, columnspan=3, padx=5, pady=5)
        self.txtarea_snp.configure(font=(24))
        
        def set_founder():
            if ent_founder.get() == "":
                showerror("ERROR", "ERROR: The ID for visualization not defined!")
                return
            elif ent_founder.get().isdigit() == False:
                showerror("ERROR", "ERROR: The ID for visualization must contain digits only!")
                return
            else:
                self.founder = ent_founder.get()
        tk.Button(tab6, text="Set", command=set_founder).grid(column = 2, row = 1, padx = 5, pady = 5)

        def set_startyear():
            try:
                self.start_year = int(ent_start_year.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid integer!")
        tk.Button(tab6, text="Set", command=set_startyear).grid(column = 2, row = 2, padx = 5, pady = 5)
        def set_endyear():
            try:
                self.end_year = int(ent_end_year.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid integer!")
        tk.Button(tab6, text="Set", command=set_endyear).grid(column = 2, row = 3, padx = 5, pady = 5)
        def set_generations_before():
            try:
                self.generations_before = int(ent_generations_before.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid integer!")
        tk.Button(tab6, text="Set", command=set_generations_before).grid(column = 2, row = 4, padx = 5, pady = 5)
        def set_generations_after():
            try:
                self.generations_after = int(ent_generations_after.get())
            except ValueError:
                showerror("ERROR", f"ERROR: the value you have entered is not a valid integer!")
        tk.Button(tab6, text="Set", command=set_generations_after).grid(column = 2, row = 5, padx = 5, pady = 5)


    def open_csv_file(self):

        try:
            filetypes = (('CSV files', '*.csv'),('All files', '*.*'))
            self.filename = fd.askopenfilename(title='Open a file',filetypes=filetypes)
            self.data = pd.read_csv(self.filename,dtype=np.str_).fillna('')
        except FileNotFoundError:
            if self.filename == "":
                return
            showerror("ERROR", f"Error in opening file {self.filename}.")
            return
        except PermissionError:
            showerror("ERROR", f"You are not allowed to read {self.filename}.")
            return
        except pd.errors.EmptyDataError:
            showerror("ERROR", f"File {self.filename} is empty.")
            return
        except ValueError:
            showerror("ERROR", f"ValueError in reading file {self.filename}.")
            return

        if 'x' in self.data.columns.values:
            self.data.set_index('x',inplace=True)

        unnamed = list(filter(lambda a: 'Unnamed' in a, self.data.columns.values))

        if len(unnamed):
            self.data.drop(columns=unnamed,inplace=True)

        self.data.replace({'father':{'':'0'},
                           'mother':{'':'0'},
                           'YOB':{'':'MISSING_YEAR'}},inplace=True)
        if 'available' in self.data.columns.values:
            self.data.replace({'available':{'':'0'}},inplace=True)
        #self.MyPedigree = Pedigree()
        #self.MyPedigree(self.data)
        if self.make_maps():
            return

        self.IDlist, \
        self.FatherMap, \
        self.MotherMap, \
        self.YobMap, \
        self.GenderMap,\
        FatalError = check_errors(self.IDlist, self.FatherMap, self.MotherMap, self.YobMap, self.GenderMap, "gui")


        if FatalError:
            showerror("ERROR", "ERROR: There are fatal errors in pedigree.")

        return

    def open_verif_report(self):
        try: 
            txtfile = "OutputVerif_Summary.txt"
            with open(txtfile, 'r') as txt:
                content = txt.read()
        except FileNotFoundError:
            showerror("ERROR", f"File {txtfile} not found!")
            return
        except PermissionError:
            showerror(
                "ERROR", f"You are not allowed to read {txtfile}.")
            return
        self.txtarea_verif.insert(tk.END, content)

        return

    def delete_verif_report(self):
        self.txtarea_verif.delete("1.0", tk.END)

    def open_stat_report(self):
        try:
            if self.Lineage.get() == 'maternal':
                txtfile = "OutputStat_DamLineMembership_1.txt"
            elif self.Lineage.get() == 'paternal':
                txtfile = "OutputStat_SireLineMembership_1.txt"
            with open(txtfile, 'r') as txt:
                content = ""
                for line in txt.readlines():
                    if 'founder' in line:
                        content += line

        except FileNotFoundError:
            showerror("ERROR", f"File {txtfile} not found!")
            return
        except PermissionError:
            showerror("ERROR", f"You are not allowed to read {txtfile}.")
            return
        try:
            self.txtarea_stat.delete(0, tk.END)
        except tk.TclError:
            pass
        self.txtarea_stat.insert(tk.END, content)

        return

    def delete_stat_report(self):
        self.txtarea_stat.delete("1.0", tk.END)

    def open_calc_report(self):
        try:
            txtfile = "OutputCalc_InputAndResults.txt"

            with open(txtfile, 'r') as txt:
                content = txt.read()

        except FileNotFoundError:
            showerror("ERROR", f"File {txtfile} not found!")
            return
        except PermissionError:
            showerror("ERROR", f"You are not allowed to read {txtfile}.")
            return
        self.txtarea_calc.insert(tk.END, content)

        return

    def delete_calc_report(self):
        self.txtarea_calc.delete("1.0", tk.END)

    def open_sampl_report(self):
        try:
            txtfile = "OutputSampl_DetailedInfo.txt"

            with open(txtfile, 'r') as txt:
                content = txt.read()

        except FileNotFoundError:
            showerror("ERROR", f"File {txtfile} not found!")
            return
        except PermissionError:
            showerror("ERROR", f"You are not allowed to read {txtfile}.")
            return
        self.txtarea_sampl.insert(tk.END, content)

        return

    #def show_image(self):
    #    try:
    #        #global img
    #        filetypes = [('PNG Image', '.png')]
    #        img_filename = fd.askopenfilename(
    #            title='Select image', filetypes=filetypes)
    #        print(img_filename)
    #        image_file = Image.open(img_filename)
    #        self.img = ImageTk.PhotoImage(image_file)
    #        self.canvas_ped.create_image(0, 0, anchor=tk.NW, image=self.img)
    #        self.canvas_ped.image = self.img
    #    except FileNotFoundError:
    #        if img_filename == "":
    #            return
    #        showerror("ERROR", f"Error in opening {img_filename}!")
    #        return
    #    except PermissionError:
    #        showerror("ERROR", f"You are not allowed to read {img_filename}.")
    #        return
    #    except AttributeError:
    #        showerror("ERROR", f"Image file not selected.")
    #        return
    #    return

    def delete_sampl_report(self):
        self.txtarea_sampl.delete("1.0", tk.END)

    def load_final_report(self):
        try:
            filetypes = (('TXT files', '*.txt'), ('All files', '*.*'))
            cols_to_use = ['SNP Name','Sample ID','Allele1 - Forward','Allele2 - Forward','GC Score']
            report_filename = fd.askopenfilename(title='Open a file', filetypes=filetypes)
            #self.final_table = pd.read_table(report_filename, sep='\t', skiprows=9, header=0, usecols=cols_to_use)
            self.final_table = dd.read_table(report_filename, sep='\t', skiprows=9, header=0, usecols=cols_to_use,dtype=str)

        except FileNotFoundError:
            if report_filename == "":
                return
            showerror("ERROR", f"Error in opening {report_filename}.")
            return
        except PermissionError:
            showerror("ERROR", f"You are not allowed to read {report_filename}.")
            return
        except pd.errors.EmptyDataError:
            showerror("ERROR", f"File {report_filename} is empty.")
            return
        except ValueError:
            showerror("ERROR", f"ValueError in reading {report_filename}.")
            return
        except IsADirectoryError:
            showwarning("", f"No file selected.")
            return

    def load_snp_map(self):
        try:
            filetypes = (('TXT files', '*.txt'), ('All files', '*.*'))
            cols_to_use = ['Name', 'Chromosome', 'Position', 'GenTrain Score']
            snp_map_file = fd.askopenfilename(
                title='Open a file', filetypes=filetypes)
            self.snp_map = dd.read_table(snp_map_file, sep='\t', header=0, usecols=cols_to_use,dtype=str)

        except FileNotFoundError:
            if snp_map_file == "":
                return
            showerror("ERROR", f"Error in opening {snp_map_file}.")
            return
        except PermissionError:
            showerror(
                "ERROR", f"You are not allowed to read {snp_map_file}.")
            return
        except pd.errors.EmptyDataError:
            showerror("ERROR", f"File {snp_map_file} is empty.")
            return
        except ValueError:
            showerror("ERROR", f"ValueError in reading {snp_map_file}.")
            return
        except IsADirectoryError:
            showwarning("", f"No file selected.")
            return

    def load_snp_list(self):
        try:
            filetypes = (('TXT files', '*.txt'), ('All files', '*.*'))
            snp_310_file = fd.askopenfilename(
                title='Open a file', filetypes=filetypes)
            self.SNPs_310 = pd.read_csv(snp_310_file,usecols=["SNP_Name"],dtype=str)

        except FileNotFoundError:
            if snp_310_file == "":
                return
            showerror("ERROR", f"Error in opening file {snp_310_file}.")
            return
        except PermissionError:
            showerror(
                "ERROR", f"You are not allowed to read {snp_310_file}.")
            return
        except pd.errors.EmptyDataError:
            showerror("ERROR", f"File {snp_310_file} is empty.")
            return
        except ValueError:
            showerror("ERROR", f"File {snp_310_file} must contain header SNP_Name.")
            return
        except IsADirectoryError:
            showwarning("", f"No file selected.")
            return

    def load_population(self):
        try:
            filetypes = (('TXT files', '*.txt'), ('All files', '*.*'))
            cols_to_use = ['Sample_ID', 'Population']
            population_file = fd.askopenfilename(title='Open a file', filetypes=filetypes)
            self.samples_population = pd.read_table(population_file, sep='\t', header=0, usecols=cols_to_use,dtype=str)

        except FileNotFoundError:
            if population_file == "":
                return
            showerror("ERROR", f"Error in opening file {population_file}.")
            return
        except PermissionError:
            showerror(
                "ERROR", f"You are not allowed to read {population_file}.")
            return
        except pd.errors.EmptyDataError:
            showerror("ERROR", f"File {population_file} is empty.")
            return
        except ValueError:
            showerror("ERROR", f"ValueError in reading {population_file}.")
            return
        except AssertionError:
            showerror("ERROR", f"AssertionError in reading {population_file}.")
            return
        except IsADirectoryError:
            showwarning("", f"No file selected.")
            return

    def run_snp(self):

        try:
            self.taurus_subset, self.proba = analyze_SNP(self.final_table,
                                                         self.snp_map, 
                                                         self.SNPs_310,
                                                         self.samples_population,
                                                         self.gc_score_thresh,
                                                         self.gt_score_thresh,
                                                         self.snp_position.get(),
                                                         "gui")
        except AttributeError:
            showerror("ERROR", f"You have not loaded all 4 necessary files!")

    def write_fasta(self):
        filepath = fd.asksaveasfilename(
            defaultextension=".fas",
            filetypes=[("FASTA Files", "*.fas"), ("All Files", "*.*")],
        )
        if not filepath:
            return
        try:
            save_fasta(filepath, self.proba)
        except AttributeError:
            showerror(
                "ERROR", f"You have not generated data to save to {filepath}!")
        return

    def read_fasta(self):
        id2seq = {}
        try:
            filetypes = (('FASTA files', '*.fasta'), ('All files', '*.*'))
            fasta_name = fd.askopenfilename(
                title='Open a file', filetypes=filetypes)
            with open(fasta_name, 'r') as fasta_file:
                for line in fasta_file:
                    if line[0] == ">":
                        sample_id = line[1:].strip()
                    else:
                        snp_sequence = line.strip()
                        id2seq[sample_id] = snp_sequence

        except FileNotFoundError:
            if fasta_name == "":
                return id2seq
            showerror("ERROR", f"Error in opening file {fasta_name}.")
            return id2seq
        except PermissionError:
            showerror(
                "ERROR", f"You are not allowed to read {fasta_name}.")
            return id2seq
        
        return id2seq

    def run_classify(self):

        if self.sequence_source.get() == "PED":
            if len(self.HaplotypeMap) == 0:
                showerror("ERROR", f"No haplotypes found in the pedigree!.")
                return
            self.classifier = SNPClassifier("xgb_main.json",self.HaplotypeMap)
            self.classifier("gui")
        elif self.sequence_source.get() == "FASTA":
            id2seq = self.read_fasta()
            if len(id2seq) == 0:
                return
            self.classifier = SNPClassifier("xgb_main.json",id2seq)
            self.classifier("gui")
        pass

    def run_to_haplogroup(self):
        if len(self.data) == 0:
            showerror("ERROR", f"Pedigree for inserting predictions not loaded!")
            return
        try:
            self.data["haplotype"] = self.data["ID"].map(self.classifier.id2pred)
        except AttributeError:
            showerror("ERROR", f"ERROR: classification not performed.")
            return

        pass

    def mutations_check(self):
            
        try:
            check_mutations(self.taurus_subset,"gui")
        except AttributeError:
            showerror(
                "ERROR", f"You have to run SNP analysis first!")
            return
        try:
            txtfile = "MUTATIONS_ALERT.TXT"

            with open(txtfile, 'r') as txt:
                content = txt.read()

        except FileNotFoundError:
            showerror("ERROR", f"File {txtfile} not found!")
            return
        except PermissionError:
            showerror("ERROR", f"You are not allowed to read {txtfile}.")
            return
        self.txtarea_snp.insert(tk.END, content)


    def make_maps(self):

        try:
            self.IDlist = list(self.data["ID"])
            self.FatherMap = dict(zip(self.data.ID, self.data.father))
            self.MotherMap = dict(zip(self.data.ID, self.data.mother))
            self.YobMap = dict(zip(self.data.ID, self.data.YOB))
            self.GenderMap = dict(zip(self.data.ID, self.data.gender))
        except KeyError:
            msg = "Mandatory column names are ID, father, mother, YOB and gender!\nCorrect your pedigree accordingly!"
            showerror("ERROR", msg)
            return 1
        except AttributeError:
            msg = "Mandatory column names are ID, father, mother, YOB and gender!\nCorrect your pedigree accordingly!"
            showerror("ERROR", msg)
            return 1
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

        return 0

    def display_dataframe(self):
        
        self.df_display = tk.Toplevel(self)
        self.df_display.geometry("750x250")
        self.df_display.title("Pedigree")
        #self.f2.pack(fill="both",expand=1)
        try:
            self.table = Table(self.df_display, dataframe=self.data,read_only=True)
            self.table.show()
        except pd.errors.EmptyDataError:
            showerror("ERROR", "ERROR: No pedigree file is loaded!")

    def run_mag_verif(self):

        try:
            if self.Lineage.get() == 'paternal':
                check_haplotype_conflicts(self.IDlist,
                                          self.FatherMap,
                                          self.HaplotypeMap,
                                          self.HaplotypedList,
                                          self.HaplotypeNamesList,
                                          self.VerifObject.get(),
                                          "gui")
            elif self.Lineage.get() == 'maternal':
                check_haplotype_conflicts(self.IDlist,
                                          self.MotherMap,
                                          self.HaplotypeMap,
                                          self.HaplotypedList,
                                          self.HaplotypeNamesList,
                                          self.VerifObject.get(),
                                          "gui")
        except AttributeError:
            showerror("ERROR", "ERROR: No pedigree file is loaded!")
 
        return
 
    def run_impute(self):

        if not len(self.HaplotypedList):
            showerror("ERROR", "ERROR: No haplotypes present in the pedigree!")
            return


        ImputedHaplotypeMap,\
        VerifiedList = impute_haplotype(self.IDlist,
                                          self.Lineage.get(),
                                          self.GenderMap,
                                          self.MotherMap,
                                          self.HaplotypeMap,
                                          self.HaplotypedList,
                                          self.HaplotypeNamesList,
                                          self.VerifObject.get(),
                                          self.ImputReliability.get(),
                                          "gui")

        self.data['haplotype'] = self.data['ID'].map(ImputedHaplotypeMap)
        self.data['sequenced'] = self.data['ID'].isin(self.HaplotypedList).astype(int)
        if self.ImputReliability.get() == 'low':
            self.data['verified'] = self.data['ID'].isin(VerifiedList).astype(int)

        return

    def run_mag_stat(self):
        try:
            if self.Lineage.get() == 'paternal':
                mag_stat(self.IDlist,
                         self.FatherMap,
                         self.YobMap,
                         self.GenderMap,
                         self.HaplotypeMap,
                         self.HaplotypedList,
                         self.FirstRefYear,
                         self.LastRefYear,
                         self.Lineage.get(),
                         "gui")
            elif self.Lineage.get() == 'maternal':
                mag_stat(self.IDlist,
                         self.MotherMap,
                         self.YobMap,
                         self.GenderMap,
                         self.HaplotypeMap,
                         self.HaplotypedList,
                         self.FirstRefYear,
                         self.LastRefYear,
                         self.Lineage.get(),
                         "gui")
        except AttributeError:
            showerror("ERROR", "ERROR: No pedigree file is loaded!")
 
        return

    def run_mag_calc(self):
        try:
            mag_calc(self.IDlist,
                     self.FatherMap,
                     self.MotherMap,
                     self.YobMap,
                     self.GenderMap,
                     self.HaplotypeMap,
                     self.HaplotypedList,
                     self.HaplotypeNamesList,
                     self.FirstRefYear,
                     self.LastRefYear,
                     self.Lineage.get(),
                     "gui")
        except AttributeError:
            showerror("ERROR", "ERROR: No pedigree file is loaded!")
            return

    def run_mag_sampl(self):
        try:
            if self.Lineage.get() == 'paternal':
                mag_sampler = MagSampl(self.IDlist,
                                       self.FatherMap,
                                       self.YobMap,
                                       self.GenderMap,
                                       self.HaplotypeMap,
                                       self.HaplotypedList,
                                       self.HaplotypeNamesList,
                                       self.AvailableMap,
                                       self.Lineage.get())
            elif self.Lineage.get() == 'maternal':
                mag_sampler = MagSampl(self.IDlist,
                                       self.MotherMap,
                                       self.YobMap,
                                       self.GenderMap,
                                       self.HaplotypeMap,
                                       self.HaplotypedList,
                                       self.HaplotypeNamesList,
                                       self.AvailableMap,
                                       self.Lineage.get())
        except AttributeError:
            showerror("ERROR", "ERROR: No pedigree file is loaded!")
            return


        mag_sampler(self.FirstRefYear,self.LastRefYear,self.SamplingMethod.get(),self.K,"gui")
 

    def run_recode(self):
        try:
            self.recoded_data = recode(self.data)
        except UnicodeDecodeError:
            showerror("ERROR", "ERROR: No pedigree file is loaded!")
            return

    def run_line_print(self):

        self.lineage = line_print(self.data,
                                  self.founder,
                                  self.start_year,
                                  self.end_year,
                                  self.generations_before,
                                  self.generations_after,
                                  self.Lineage.get(),
                                  "gui",
                                  "R_plots.png")

    def save_recoded_csv(self):
        filepath = fd.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")],
        )
        if not filepath:
            return
        try:
            self.recoded_data.to_csv(filepath,index=False)
        except AttributeError:
            showerror("ERROR", f"AttributeError in saving {filepath}.")
            pass
        except NameError:
            showerror(title="No data error", 
                      message="ERROR: no recoded data generated, run recode!")

    def save_imputed_csv(self):
        showwarning("WARNING", f"Be careful not to overwrite your original pedigree!.")
        filepath = fd.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")],
        )
        if not filepath:
            return
        try:
            self.data.to_csv(filepath,index=False)
        except AttributeError:
            showerror("ERROR", f"AttributeError in saving {filepath}.")
            pass
        except NameError:
            showerror(title="No data error", 
                      message="ERROR: no data generated, load them and impute them!")

    def save_lineage_csv(self):
        filepath = fd.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv"), 
                       ("TSV Files", "*.tsv"), 
                       ("All Files", "*.*")],
        )
        if not filepath:
            return

        try:
            if path.splitext(filepath)[1] == '.csv':
                self.lineage.to_csv(filepath,index=False)
            elif path.splitext(filepath)[1] == '.tsv':
                self.lineage.to_csv(filepath,index=False,sep="\t")
            else:
                self.lineage.to_string(filepath,index=False)
        except AttributeError:
            showerror(title="No data error",
                      message="ERROR: no lineage data present, file not created!")
            return
        except NameError:
            showerror(title="No data error", 
                      message="ERROR: no lineage generated, run print_line!")

def runGUI():
 
    mag_gui = MagellanGUI()
    mag_gui.mainloop()

if __name__ == "__main__":
    runGUI()
