# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hršak and Strahil Ristov

import argparse
import pandas as pd
import numpy as np
import sys

def build_tree(data,founder,ancestor_line):

    graph_tree = []
    ancestors = [founder]
    ancestor_offspring_pairs = data.loc[data[ancestor_line].isin(ancestors)]#[['ID',ancestor_line]]

    while len(ancestor_offspring_pairs):
        for row in ancestor_offspring_pairs.values.tolist():
           graph_tree.append(row)
        ancestors = list(ancestor_offspring_pairs['ID'])
        ancestor_offspring_pairs = data.loc[data[ancestor_line].isin(ancestors)]#[['ID',ancestor_line]]

    return np.array(graph_tree)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='input pedigree file')
    parser.add_argument('--output', dest='output', default='', help='define output name (without .csv extension)')
    parser.add_argument('--founder', dest='founder', default='', help='define which founder line to analyze')
    parser.add_argument('--gender', dest='gender', default='female', help='define whether we follow female or male line')
    parser.add_argument('--start-year',default=0,dest='startyear',type=int,
                        help='''Give start year for the reference population.''')
    parser.add_argument('--end-year', default=0, dest='endyear', type=int,
                        help='''Give end year for the reference population.''')

    args = parser.parse_args()

    if args.gender == 'female':
        ancestor_line = 'mother'
        gender = '2'
    elif args.gender == 'male':
        ancestor_line = 'father'
        gender = '1'
    else:
        print("ERROR: Unrecognized option for gender! Allowed options: male, female.")
        sys.exit()

    if args.founder == '':
        print("ERROR: Founder for visualization not defined!")
        sys.exit()

    if (args.endyear != 0) and (args.startyear > args.endyear):
        print('ERROR: Defined start year {0} is later than defined end year {1}!'.format(args.startyear,args.endyear))
        sys.exit()

    if args.output == '':
        output = args.founder
    else:
        output = args.output

    data_types = {'x': str,
                  'ID': str,
                  'father': str,
                  'mother': str,
                  'YOB': int,
                  'gender': str,}

    cols_to_use = ['x','ID','father','mother','YOB','gender']

    try:
        data = pd.read_csv(args.filename, usecols=cols_to_use, dtype=data_types).fillna('0')
    except FileNotFoundError:
        print("File {0} not found.".format(args.filename))
    except pd.errors.EmptyDataError:
        print("File {0} is empty.".format(args.filename))

    data = data.loc[data['gender'] == gender]

    if args.founder not in list(data.ID):
        print("ERROR: Individual {0} is not defined in the pedigree!".format(args.founder))
        print("Choose another founder individual.")
        sys.exit()

    if (list(data[data.ID == args.founder][ancestor_line])[0] != '0'):
        print("ERROR: Individual {0} is not a founder in the {1} line!".format(args.founder,args.gender))
        print("Choose another individual that is a real founder, i.e. has no defined {0}.".format(ancestor_line))
        sys.exit()

    graph_tree = build_tree(data,args.founder,ancestor_line)

    try:
        data = pd.DataFrame(graph_tree, columns=cols_to_use, dtype=str)
    except ValueError:
        print('Founder {0} has no {1} offspring, no output is produced.'.format(
            args.founder, args.gender))
        sys.exit()

    data = data.astype(data_types)

    if args.startyear != 0:
        data = data.loc[data['YOB'] >= args.startyear]
        output = output + '_from{0}'.format(args.startyear)
    if args.endyear != 0:
        data = data.loc[data['YOB'] <= args.endyear]
        output = output + '_to{0}'.format(args.endyear)

    if len(data) == 0:
        print('Founder {0} has no {1} offspring born from {2} to {3}.'.format(args.founder,args.gender,args.startyear,args.endyear))
        print('No output is produced.')
        sys.exit()

    output = output + '.csv'

    data.to_csv(output,index=False)
    


if __name__ == "__main__":
    main()
