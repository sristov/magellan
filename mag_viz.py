# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov



import platform

if "Windows" not in platform.system():
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import StrVector, IntVector
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects import r
    import rpy2.rinterface_lib.callbacks

import pandas as pd
import numpy as np
import sys
from tkinter import messagebox, simpledialog


def build_tree(data,founder,ancestor_line):

    graph_tree = data.loc[data['ID']==founder].values.tolist()
    ancestors = [founder]
    ancestor_offspring_pairs = data.loc[data[ancestor_line].isin(ancestors)]#[['ID',ancestor_line]]

    generation = 0
    generation_list = [0]
    while len(ancestor_offspring_pairs):
        generation += 1
        for row in ancestor_offspring_pairs.values.tolist():
            graph_tree.append(row)
        generation_list.extend(len(ancestor_offspring_pairs)*[generation])
        ancestors = list(ancestor_offspring_pairs['ID'])
        ancestor_offspring_pairs = data.loc[data[ancestor_line].isin(ancestors)]#[['ID',ancestor_line]]

    affected = np.ones(np.shape(np.array(graph_tree)[:,0]),dtype=int)#1*(np.array(graph_tree)[:,0] == ID)
    graph_tree = np.hstack((graph_tree,np.expand_dims(generation_list,axis=1),np.expand_dims(affected,axis=1)))

    return np.array(graph_tree)

def fix_data(data,lineage,ID):
    
    #kinship2 requires all members to be defined and at least 2 founders

    founder_fathers = pd.DataFrame(data.loc[~data['father'].isin(data['ID'])]['father'])
    founder_mothers = pd.DataFrame(data.loc[~data['mother'].isin(data['ID'])]['mother'])

    founder_fathers.rename(columns={'father':'ID'},inplace=True)
    founder_mothers.rename(columns={'mother':'ID'},inplace=True)

    # remove duplicate items
    founder_fathers.drop_duplicates(subset=['ID'],inplace=True)
    founder_mothers.drop_duplicates(subset=['ID'],inplace=True)
    
    founder_fathers.insert(1,'father','0')
    founder_fathers.insert(2,'mother','0')
    founder_fathers.insert(3,'gender','1')
    founder_fathers.insert(4,'oldID','0')
    founder_mothers.insert(1, 'father', '0')
    founder_mothers.insert(2, 'mother', '0')
    founder_mothers.insert(3, 'gender', '2')
    founder_mothers.insert(4,'oldID','0')
    if lineage == 'paternal':
        founder_fathers.insert(5,'affected','1')
        founder_mothers.insert(5,'affected','0')
    elif lineage == 'maternal':
        founder_fathers.insert(5,'affected','0')
        founder_mothers.insert(5,'affected','1')


    founder_fathers.drop(founder_fathers.loc[founder_fathers['ID'] == '0'].index, inplace=True)
    founder_mothers.drop(founder_mothers.loc[founder_mothers['ID'] == '0'].index, inplace=True)

    data = pd.concat([data, founder_fathers], ignore_index=True)
    data = pd.concat([data, founder_mothers], ignore_index=True)

    dummyID = "9"*(len(ID)+1)

    filter1 = (data["father"] != "0") & (data["mother"] == "0")
    filter2 = (data["father"] == "0") & (data["mother"] != "0")
    # kinship2 doesn't allow only one parent being defined
    # insert a dummy missing parent
    if lineage == 'paternal':
        data['mother'] = data['mother'].where(~filter1,other=dummyID)
        data['mother'] = data['mother'].where(~filter2,other='0')
        new_row = {'ID':[dummyID],'father':['0'],'mother':['0'],'gender':['2'],'oldID':['0'],'affected':['0']}
        new_row = pd.DataFrame(new_row)
        data = pd.concat([data,new_row],ignore_index = True)
    elif lineage == 'maternal':
        data['father'] = data['father'].where(~filter2, other=dummyID)
        data['father'] = data['father'].where(~filter1, other='0')
        new_row = {'ID':[dummyID],'father':['0'],'mother':['0'],'gender':['1'],'oldID':['0'],'affected':['0']}
        new_row = pd.DataFrame(new_row)
        data = pd.concat([data,new_row],ignore_index = True)

    #ID to id, gender to sex
    data.rename(columns={'ID':'id','gender':'sex'},inplace=True)

    return data


def visualize_line(data,mode,image_filename):
    utils = importr('utils')
    utils.chooseCRANmirror(ind=1)

    if not rpackages.isinstalled('kinship2'):
        utils.install_packages(StrVector(['kinship2']))

    kinship2 = importr('kinship2')

    r_id = IntVector(data.id)
    r_father = IntVector(data.father)
    r_mother = IntVector(data.mother)
    r_sex = IntVector(data.sex)
    r_affected = IntVector(data.affected)

    pedAll = kinship2.pedigree(id=r_id,dadid=r_father,momid=r_mother,sex=r_sex,affected=r_affected)

    #plot pedigree
    if mode == "gui":
        image_file = simpledialog.askstring(
            "Save image", "Give filename for PNG image.", initialvalue=image_filename)
    elif mode == "cl":
        image_file = image_filename
    r.png(image_file,width=2100,height=1500)

    buf1 = []
    def f1(x):
        buf1.append(x)

    buf2 = []
    def f2(x):
        buf2.append(x)

    consolewrite_print_backup = rpy2.rinterface_lib.callbacks.consolewrite_print
    rpy2.rinterface_lib.callbacks.consolewrite_print = f1

    consolewrite_warn_backup = rpy2.rinterface_lib.callbacks.consolewrite_warnerror
    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = f2
    
    grdevices = importr('grDevices')
    kinship2.plot_pedigree(pedAll,col=4)  
    grdevices.dev_off()

    output_message = "".join(buf1)
    if len(output_message):
        if mode == "gui":
            messagebox.showinfo("INFO", output_message)
        elif mode == "cl":
            sys.stdout.write(output_message+'\n')

    output_warning = "".join(buf2)
    if len(output_warning):
        if mode == "gui":
            messagebox.showwarning("WARNING", output_warning)
        elif mode == "cl":
            sys.stderr.write(f"WARNING: {output_warning}\n")

    # restore default function
    rpy2.rinterface_lib.callbacks.consolewrite_print = consolewrite_print_backup
    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = consolewrite_warn_backup

    return


def line_print(data, ID, startyear, endyear, gen_before, gen_after, lineage, mode, image_filename):

    if len(data) == 0 and mode == "gui":
        messagebox.showerror("ERROR", "No pedigree is loaded!")
        return pd.DataFrame()

    if ID == '' and mode == "gui":
        messagebox.showerror("ERROR", "Individual for visualization not defined!")
        return pd.DataFrame()

    if (endyear != 0) and (startyear > endyear):
        if mode == "gui":
            messagebox.showerror("ERROR", f"Defined start year {startyear} is later than defined end year {endyear}!")
        elif mode == "cl":
            sys.stderr.write(f"ERROR: Defined start year {startyear} is later than defined end year {endyear}!\n")
        return pd.DataFrame()

    try:
        if list(data.loc[data['ID'] == ID]['gender'])[0] == '2':
            ancestor_line = 'mother'
            gender = '2'
            gender1 = '2'
            str_gender = 'female'
        elif list(data.loc[data['ID'] == ID]['gender'])[0] == '1':
            if lineage == 'maternal':
                ancestor_line = 'mother'
                gender1 = '2'
                str_gender = 'female'
            elif lineage == 'paternal':
                ancestor_line = 'father'
                gender1 = '1'
                str_gender = 'male'
            gender = '1'
    except IndexError:
        if mode == "gui":
            messagebox.showerror("ERROR", f"Individual {ID} not found in the pedigree!")
        elif mode == "cl":
            sys.stderr.write(f"ERROR: Individual {ID} not found in the pedigree!\n")
        return pd.DataFrame()

    allowed_column_names = ['x','ID','father','mother','YOB','gender','haplotype','available','oldID']
    cols_to_use = [name for name in data.columns.values if name in allowed_column_names]

    data_types = {}
    for col in cols_to_use:
        if col == 'YOB' or col == 'generation' or col == 'affected':
            data_types[col] = int
        else:
            data_types[col] = str

    if gender == '1' and ancestor_line == 'father':
        data = data.loc[data['gender'] == gender]

    # check if ID is defined in ID column
    if len(list(data[data.ID == ID][ancestor_line])):
        if (list(data[data.ID == ID][ancestor_line])[0] != '0'):
            current = ID
            ancestor = list(data[data.ID == current][ancestor_line])[0]
            while ancestor != '0':
                current = ancestor
                # check if current is defined in ID column, if not add it
                if len(list(data[data.ID == current][ancestor_line])):
                    ancestor = list(data[data.ID == current][ancestor_line])[0]
                else:
                    ancestor = '0'
                    new_row = pd.DataFrame({'ID':[current],'father':['0'],'mother':['0'],'gender':[gender1],'YOB':['0']})
                    data = pd.concat([data,new_row],ignore_index=True)
            founder = current
            message = f"Individual {ID} is not a founder in the {str_gender} line.\nThe founder is {founder}."
            if mode == "gui":
                messagebox.showinfo("INFO", message)
            elif mode == "cl":
                sys.stdout.write(f"INFO: {message}\n")
        else:
            founder = ID
            message = f"Individual {ID} is the founder in the {str_gender} line."
            if mode == "gui":
                messagebox.showinfo("INFO", message)
            elif mode == "cl":
                sys.stdout.write(f"INFO: {message}\n")
    else:
        founder = ID
        # add missing founder entry
        new_row = pd.DataFrame({'ID':[ID],'father':['0'],'mother':['0'],'gender':[gender1],'YOB':['0']})
        data = pd.concat([data,new_row],ignore_index=True)
        message = f"Individual {ID} is the founder in the {str_gender} line."
        if mode == "gui":
            messagebox.showinfo("INFO", message)
        elif mode == "cl":
            sys.stdout.write(f"INFO: {message}\n")

    try:
        graph_tree = pd.DataFrame(build_tree(data,founder,ancestor_line),columns=cols_to_use+['generation','affected'])
    except ValueError:
        message = f"Founder {founder} has no {str_gender} offspring, no output is produced."
        if mode == "gui":
            messagebox.showerror("ERROR", message)
        elif mode == "cl":
            sys.stderr.write(f"ERROR: {message}\n")
        return pd.DataFrame()
 
    graph_tree = graph_tree.astype({'generation':int})
    ID_generation = list(graph_tree[graph_tree.ID == ID]['generation'])[0]
    if gen_before != 0:
        first_gen = int(ID_generation) - int(gen_before)
        graph_tree = graph_tree[graph_tree['generation'] >= first_gen]
    if gen_after != 0:
        last_gen = int(ID_generation) + int(gen_after)
        graph_tree = graph_tree[graph_tree['generation'] <= last_gen]

    try:
        graph_tree.drop(columns='generation',inplace=True)
    except KeyError:
        if mode == "gui":
            messagebox.showerror("ERROR", "KeyError when erasing generation column.")
            pass
        elif mode == "cl":
            sys.stderr.write("ERROR: KeyError when erasing generation column.\n")
            pass

    try:
        data = pd.DataFrame(graph_tree, columns=cols_to_use+["affected"], dtype=str)
    except ValueError:
        message = f"Founder {founder} has no {str_gender} offspring, no output is produced."
        if mode == "gui":
            messagebox.showerror("ERROR", message)
        elif mode == "cl":
            sys.stderr.write(f"ERROR: {message}\n")
        return pd.DataFrame()

    data = data.astype(data_types)

    if startyear != 0:
        data = data.loc[data['YOB'] >= int(startyear)]
    if endyear != 0:
        data = data.loc[data['YOB'] <= int(endyear)]

    if len(data) <= 1:
        if ancestor_line == 'mother':
            message = f"Founder {founder} has no offspring born from {startyear} to {endyear}."
        elif ancestor_line == 'father':
            message = f"Founder {founder} has no {str_gender} offspring born from {startyear} to {endyear}."
        if mode == "gui":
            messagebox.showerror("ERROR", message)
        elif mode == "cl":
            sys.stderr.write(f"ERROR: {message}\n")
        return pd.DataFrame()

    #try:
    #    data.drop(columns=['YOB'], inplace=True)
    #except KeyError:
    #    pass

    #try:
    #    data.drop(columns=['haplotype'], inplace=True)
    #except KeyError:
    #    pass

    data = fix_data(data,lineage,ID)
    print(data)
    if "Windows" not in platform.system():
        visualize_line(data, mode, image_filename)
    else:
        message = f"You are running Magellan on {platform.system()}.\nVisualization is not available but you can save the lineage to file."
        if mode == "gui":
            messagebox.showinfo("INFO", message)
        elif mode == "cl":
            sys.stdout.write(f"INFO: {message}\n")

    return data
