# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

import sys
import json

def inp_int(msg):
    while True:
        try:
            answer = input(msg)
        except KeyboardInterrupt:
            sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
            sys.exit()

        if answer == "":
           value = 0
           break
        else:
            try:
                value = int(answer)
                break
            except ValueError:
                sys.stderr.write("ERROR: the value you have entered is not a valid integer!\n")
                continue
    return value
    
def inp_float(msg):
    while True:
        try:
            answer = input(msg)
        except KeyboardInterrupt:
            sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
            sys.exit()

        if answer == "":
           value = 0.0
           break
        else:
            try:
                value = float(answer)
                break
            except ValueError:
                sys.stderr.write("ERROR: the value you have entered is not a valid floating point number!\n")
                continue
    return value

def inp_str(msg,allowed_values):
    while True:
        try:
            answer = input(msg)
        except KeyboardInterrupt:
            sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
            sys.exit()

        if answer == "":
            value = allowed_values[0]
            break
        if answer.lower() not in allowed_values:
            sys.stderr.write(f"ERROR: Option {answer} not recognized!\n")
            continue
        else:
            value = answer.lower()
            break
    return value

def inp_yes_no(msg="",warnmsg=""):
    while True:
        try:
            answer = input(msg)
        except KeyboardInterrupt:
            sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
            sys.exit()

        if answer.capitalize() == "Yes":
            value = True
            break
        elif answer.capitalize() == "No":
            value = False
            sys.stdout.write(warnmsg)
            break
        else:
            sys.stderr.write(f"ERROR: Option {answer} not recognized!\nAllowed options are Yes and No!")
            continue
    return value

def inp_lineage():
    while True:
        try:
            lineage = input("What lineage do you want to analyze? [maternal, paternal]\n")
        except KeyboardInterrupt:
            sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
            sys.exit()
        if lineage.lower() != "maternal" and lineage.lower() != "paternal":
            sys.stderr.write(f"ERROR: Option {lineage} not recognized!\n")
            continue
        else:
            value = lineage.lower()
            break
    return value


sys.stdout.write("+==========================================+\n")
sys.stdout.write("|    MaGelLan2.0 job preparation script    |\n")
sys.stdout.write("+==========================================+\n\n\n")

input_parameters = {}
try:
    pedigree_filename = input("Provide the input pedigree file. Press Enter for no file\n")
except KeyboardInterrupt:
    sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
    sys.exit()

if len(pedigree_filename):

    input_parameters["pedigree"] = pedigree_filename

    input_parameters["mag_verif"] = inp_yes_no(msg="Do you want to run Magellan verification module? [Yes, No]\n",
                                               warnmsg="WARNING: mag_verif will not be executed, be sure your pedigree has no internal conflicts!\n")

    input_parameters["mag_imput"] = inp_yes_no(msg="Do you want to impute haplotypes? [Yes, No]\n")

    input_parameters["mag_stat"] = inp_yes_no(msg="Do you want to run Magellan statistics module? [Yes, No]\n")

    input_parameters["mag_calc"] = inp_yes_no(msg="Do you want to run Magellan calculation module? [Yes, No]\n")

    input_parameters["mag_sampl"] = inp_yes_no(msg="Do you want to run Magellan sampling module? [Yes, No]\n")

    input_parameters["mag_recode"] = inp_yes_no(msg="Do you want to recode your pedigree for visualization? [Yes/No]\n")

    input_parameters["mag_visualize"] = inp_yes_no(msg="Do you want to visualize a lineage in your pedigree? [Yes/No]\n")

else:
    sys.stdout.write("No pedigree file chosen.\n")
    sys.stdout.write("Verification, statistics, calculation, sampling, recoding and visualization modules will not be executed\n")
    input_parameters["pedigree"] = ""
    input_parameters["mag_verif"] = False
    input_parameters["mag_imput"] = False
    input_parameters["mag_stat"] = False
    input_parameters["mag_calc"] = False
    input_parameters["mag_sampl"] = False
    input_parameters["mag_recode"] = False
    input_parameters["mag_visualize"] = False


input_parameters["mag_snp"] = inp_yes_no(msg="Do you want to run SNP analysis? [Yes/No]\n")

input_parameters["mag_classify"] = inp_yes_no(msg="Do you want to run haplotype classification? [Yes/No]\n")

if input_parameters["mag_verif"] == True:
    sys.stdout.write("\n\n----------------------------------------------------------\n")
    sys.stdout.write("Defining input options for Magellan verification module...\n\n")
    input_parameters["verif_options"] = {}
    
    input_parameters["verif_options"]["lineage"] = inp_str("What lineage do you want to analyze? [maternal, paternal]\n", ["maternal","paternal"])

    input_parameters["verif_options"]["verif_object"] = inp_str("Do you want to check conflicts for haplogroups or SNP sequences? [haplogroup, snpseq]\n",
                                                           ["haplogroup","snpseq"])

if input_parameters["mag_imput"] == True:
    sys.stdout.write("\n\n----------------------------------------------------------\n")
    sys.stdout.write("Defining input options for haplotype imputation...\n\n")
    input_parameters["imput_options"] = {}
    
    input_parameters["imput_options"]["lineage"] = inp_str("What lineage do you want to impute? [maternal, paternal]\n", ["maternal","paternal"])
    input_parameters["imput_options"]["reliability"] = inp_str("What level of reliability in imputation? [low, high]\n", ["low","high"])

    try:
        imput_filename = input("Provide a name for the imputed pedigree file: (press Enter if you don't want to save it)\n")
    except KeyboardInterrupt:
        sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
        sys.exit()
    input_parameters["imput_options"]["imput_filename"] = imput_filename


if input_parameters["mag_stat"] == True:
    sys.stdout.write("\n\n--------------------------------------------------------\n")
    sys.stdout.write("Defining input options for Magellan statistics module...\n\n")
    input_parameters["stat_options"] = {}
    
    input_parameters["stat_options"]["first_ref_year"] = inp_int("Provide the start year of birth for reference population (press Enter for no lower limit):\n")
    
    input_parameters["stat_options"]["last_ref_year"] = inp_int("Provide the end year of birth for reference population (press Enter for no upper limit):\n")

    input_parameters["stat_options"]["lineage"] = inp_str("What lineage do you want to analyze? [maternal, paternal]\n", ["maternal","paternal"])



if input_parameters["mag_calc"] == True:
    sys.stdout.write("\n\n---------------------------------------------------------\n")
    sys.stdout.write("Defining input options for Magellan calculation module...\n\n")
    input_parameters["calc_options"] = {}

    if input_parameters["mag_stat"] == True:

        same_settings = inp_yes_no("Do you want to put the same settings as for mag_stat? [Yes/No]\n")

        if same_settings:
            input_parameters["calc_options"] = input_parameters["stat_options"]

        else:
            
            input_parameters["calc_options"]["first_ref_year"] = inp_int("Provide the start year of birth for reference population (press Enter no lower limit):\n")
            
            input_parameters["calc_options"]["last_ref_year"] = inp_int("Provide the end year of birth for reference population (press Enter no upper limit):\n")
            
            input_parameters["calc_options"]["lineage"] = inp_str("What lineage do you want to analyze? [maternal, paternal]\n", ["maternal","paternal"])

    else:
        
        input_parameters["calc_options"]["first_ref_year"] = inp_int("Provide the start year of birth for reference population (press Enter no lower limit):\n")
            
        input_parameters["calc_options"]["last_ref_year"] = inp_int("Provide the end year of birth for reference population (press Enter no upper limit):\n")
            
        input_parameters["calc_options"]["lineage"] = inp_str("What lineage do you want to analyze? [maternal, paternal]\n", ["maternal","paternal"])


if input_parameters["mag_sampl"] == True:
    sys.stdout.write("\n\n------------------------------------------------------\n")
    sys.stdout.write("Defining input options for Magellan sampling module...\n\n")
    input_parameters["sampl_options"] = {}

    input_parameters["sampl_options"]["K"] = inp_int("How many individuals do you want to select for sampling?\n")
    
    input_parameters["sampl_options"]["first_ref_year"] = inp_int("Provide the start year of birth for reference population (press Enter no lower limit):\n")
            
    input_parameters["sampl_options"]["last_ref_year"] = inp_int("Provide the end year of birth for reference population (press Enter no upper limit):\n")
            
    input_parameters["sampl_options"]["lineage"] = inp_str("What lineage do you want to sample? [maternal, paternal]\n", ["maternal","paternal"])

    input_parameters["sampl_options"]["sampling_method"] = inp_str("What sampling algorithm do you want to use? [Default= greedy]\n",['greedy','optimal','cubic'])
    

if input_parameters["mag_recode"] == True:
    sys.stdout.write("\n\n------------------------------------------------------\n")
    sys.stdout.write("Defining input options for Magellan recoding module...\n\n")
    input_parameters["recode_options"] = {}
    recoded_filename = input("Provide the filename to save the recoded pedigree:\n")
    input_parameters["recode_options"]["recoded_filename"] = recoded_filename

if input_parameters["mag_visualize"] == True:
    sys.stdout.write("\n\n-----------------------------------------------------------\n")
    sys.stdout.write("Defining input options for Magellan visualization module...\n\n")
    input_parameters["visualize_options"] = {}
    if input_parameters["mag_recode"] == True:
        same_recoded = inp_yes_no("Is the name of the recoded file the same as for mag_recode module? [Yes/No]\n")
        if same_recoded:
            input_parameters["visualize_options"]["recoded_filename"] = input_parameters["recode_options"]["recoded_filename"]
        else:
            recoded_filename = input("Provide the filename of the recoded pedigree:\n")
            input_parameters["visualize_options"]["recoded_filename"] = recoded_filename
    else:
        recoded_filename = input("Provide the filename of the recoded pedigree:\n")
        input_parameters["visualize_options"]["recoded_filename"] = recoded_filename
    
    while True:
        unit_id = input("Provide the ID of the individual for lineage visualization (numeric ID only).\n")
        if unit_id == "":
            sys.stderr.write("ERROR: The ID for visualization not defined!\n")
            continue
        elif unit_id.isdigit() == False:
            sys.stderr.write("ERROR: The ID for visualization must contain digits only!\n")
            continue
        else:
            input_parameters["visualize_options"]["unit_id"] = unit_id
            break

    input_parameters["visualize_options"]["start_yob"] = inp_int("Provide start year of birth for visualization. (press Enter for all years)\n")
    
    input_parameters["visualize_options"]["end_yob"] = inp_int("Provide end year of birth for visualization. (press Enter for all years)\n")
    
    input_parameters["visualize_options"]["gen_before"] = inp_int(f"How many generations before {unit_id} should be visualized. (press Enter for all generations)\n")
    
    input_parameters["visualize_options"]["gen_after"] = inp_int(f"How many generations after {unit_id} should be visualized. (press Enter for all generations)\n")

    input_parameters["visualize_options"]["lineage"] = inp_str("What lineage do you want to visualize? [maternal, paternal]\n", ["maternal","paternal"])

    try:
        lineage_imagename = input("Provide a name for the PNG image file: [default: R_plots.png]\n")
    except KeyboardInterrupt:
        sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
        sys.exit()

    if lineage_imagename == "":
        input_parameters["visualize_options"]["imagename"] = "R_plots.png"
    else:
        input_parameters["visualize_options"]["imagename"] = lineage_imagename

    try:
        lineage_filename = input("Provide a name for the lineage CSV file: (press Enter if you don't want to save it)\n")
    except KeyboardInterrupt:
        sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
        sys.exit()
    input_parameters["visualize_options"]["lineage_filename"] = lineage_filename

if input_parameters["mag_snp"] == True:
    sys.stdout.write("\n\n----------------------------------------------------------\n")
    sys.stdout.write("Defining input options for Magellan SNP analysis module...\n\n")
    input_parameters["snp_options"] = {}
    try:
        report_filename = input("Provide the name of the Report file:\n")
        map_filename = input("Provide the name of the Map file:\n")
        snp_filename = input("Provide the name of the SNP list file:\n")
        population_filename = input("Provide the name of the Population list file:\n")
    except KeyboardInterrupt:
        sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
        sys.exit()
    input_parameters["snp_options"]["report_filename"] = report_filename
    input_parameters["snp_options"]["map_filename"] = map_filename
    input_parameters["snp_options"]["snp_filename"] = snp_filename
    input_parameters["snp_options"]["population_filename"] = population_filename

    input_parameters["snp_options"]["gcscore_thresold"] = inp_float("Provide the threshold for the GC Score: [default: 0.0]\n")
    
    input_parameters["snp_options"]["gtscore_thresold"] = inp_float("Provide the threshold for the GT Score: [default: 0.0]\n")
    
    input_parameters["snp_options"]["snp_position"] = inp_str("Provide the position of SNPs: [default: MT]\n",["MT", "Y", "X", "AUTO"])
    
    try:
        fasta_filename = input("Provide a name to save FASTA sequence. [Press Enter to not save]\n")
    except KeyboardInterrupt:
        sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
        sys.exit()
    input_parameters["snp_options"]["fasta_filename"] = fasta_filename
    
    input_parameters["snp_options"]["mut_check"] = inp_yes_no("Do you want to check for deleterious mutations? [Yes/No]\n")


if input_parameters["mag_classify"] == True:
    sys.stdout.write("\n\n----------------------------------------------------------\n")
    sys.stdout.write("Defining input options for Magellan SNP classification module...\n\n")
    input_parameters["classify_options"] = {}
    
    input_parameters["classify_options"]["sequence_source"] = inp_str("Provide the source of SNP sequences [default: ped]\n",["ped", "fasta"])
    
    if input_parameters["classify_options"]["sequence_source"] == "fasta":
        same_fasta = inp_yes_no("Is the fasta file same as you saved in SNP analysis? [Yes,No]\n")
        if same_fasta:
            if len(input_parameters["snp_options"]["fasta_filename"]) == 0:
                sys.stderr.write("ERROR: no FASTA file has been saved in the Magellan SNP module.\n")
                sys.exit()
            else:
                input_parameters["classify_options"]["fasta_filename"] = input_parameters["snp_options"]["fasta_filename"]
        else:
            try:
                input_parameters["classify_options"]["fasta_filename"] = input("Provide a name of the FASTA file to read.\n")
            except KeyboardInterrupt:
                sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
                sys.exit()
            
    input_parameters["classify_options"]["impute_haplogroup"] = inp_yes_no("Do you want to insert haplogroups into your pedigree? [Yes/No]\n")
    
    if input_parameters["classify_options"]["impute_haplogroup"]:
        sys.stdout.write("WARNING: be careful not to overwrite your original pedigree!\n")
        try:
            input_parameters["classify_options"]["csv_name"] = input("Provide a name of the CSV file to save the imputed pedigree.\n")
        except KeyboardInterrupt:
            sys.stderr.write("\nScript interrupted.\nConfiguration file not generated\n")
            sys.exit()
    else:
        input_parameters["classify_options"]["csv_name"] = ""

    
sys.stdout.write("+=====================================================+\n")
sys.stdout.write("| Saving the selected options to the config.json file |\n")
sys.stdout.write("+=====================================================+\n\n")

with open("config.json", "w") as jsfile:
    json.dump(input_parameters, jsfile, indent=4)


sys.stdout.write("+================================+\n")
sys.stdout.write("| Sucessfully exiting the script |\n")
sys.stdout.write("+================================+\n")
