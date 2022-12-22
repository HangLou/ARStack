import numpy as np
import pandas as pd
import os
from collections import Counter
import traceback
from tqdm import tqdm
numbers = np.arange(0, 1000, 100)
numbers = dict(zip(numbers, numbers))


def get_data(dataset):
    target_csv = pd.read_csv(
        'data/input/processed/targets_Dataset'+str(dataset)+'.csv')
    del target_csv['Unnamed: 0']
    target_csv.set_index('0', inplace=True)

    keep_track = pd.read_csv(
        'data/input/processed/atom_counts_Dataset'+str(dataset)+'.csv')
    del keep_track['Unnamed: 0']
    keep_track.set_index('0', inplace=True)
    # first_atoms=keep_track['1']

    return(target_csv, keep_track)


def get_MSAs(file, target, first_atoms, length):
    #global gs
    global final_df
    global final_solo

    pdb = file[:4]
    # print(pdb)
    df = pd.read_csv('data/input/raw/MSA_csv/'+file)
    df.set_index('Unnamed: 0', inplace=True)

    non_valid = set(df.columns)-set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                                     'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'])  # ['X', 'B', 'Z', 'J', 'U', 'O']

    if len(non_valid) > 0:
        for aa in non_valid:
            try:
                del df[aa]  # X comes from alignments and can be removed,
            except:
                non_present = 0

    correct_order = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
                     'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

#    input(df.head())
#    print(df.columns)
    gs1 = df.columns
    df = df[correct_order]
    df = df.round(5)
    gs2 = df.columns
    if list(gs1) != list(gs2):
        print(gs1, '\n', gs2)
        # input('no')
#    print(df.columns)
#    input(df.head())

    values = target.loc[pdb]
    values.dropna(how='any', inplace=True)

    # take target values for each atom type
    midvals = values[0::3]  # N
    midvals2 = values[1::3]  # CA
    midvals3 = values[2::3]  # C

    try:
        start_atom = int(first_atoms.loc[pdb])
        start_amino = int((start_atom)+1/3)
        # print(start_atom)
        # input(start_amino)
        # print(df)
        # print(len(df))

        #print(df.head(), df.columns)
        # adding tails (of 0s) at both ends to include amino acids at the start and end.
        origo = np.zeros(21)
        origo = [origo, origo, origo, origo, origo,
                 origo, origo, origo, origo, origo]
        tails = pd.DataFrame(origo)
        tails.columns = df.columns

        df = pd.concat([tails, df], ignore_index=True)
        df = pd.concat([df, tails], ignore_index=True)

        for i in range(length):
            try:
                concat = []
                for m in range(i, i+21):
                    concat += list(df.iloc[m])

                final_df.append(
                    [str(pdb), str(i), midvals[i], 1, 0, 0] + concat)
                final_df.append(
                    [str(pdb), str(i), midvals2[i], 0, 1, 0] + concat)
                final_df.append(
                    [str(pdb), str(i), midvals3[i], 0, 0, 1] + concat)

                final_solo.append([str(pdb), str(i), midvals2[i]] + concat)
                # gs[str(df.columns)]=1
            except Exception as e:
                print(traceback.print_exc())
                print('nah', e)

    except Exception as e:
        print(e, traceback.print_exc())


def MSA_df(dataset):
    target, atom_count = get_data(dataset)
    first_atoms = atom_count['1']

    #global gs
    # gs={}
    k = 0
    csv = pd.DataFrame()
    global final_df
    global final_solo
    final_df = []
    final_solo = []
    # print(len(target))
    numbers = np.arange(0, len(target), 50)
    numbers = dict(zip(numbers, numbers))

    #selection = open('..\\clean scripts combined sets\\one_prot_pfam_working.txt','r').readlines()
    selection = list(target.index)
    #selection = [i[:4].upper() for i in selection]
    #selection = ['6ILD']

    for file in tqdm(os.listdir('data/input/raw/MSA_csv/')):
        if file[:4] in selection:
            try:
                if k in numbers:
                    print(k)

                pdb = file[:4]
                # print(len(atom_count.loc[pdb].dropna()))
                length = int((len(atom_count.loc[pdb].dropna())+1)/3)

                MSAs = get_MSAs(file, target, first_atoms, length)

                k += 1
            except Exception as e:
                print('outer', e)

    df2 = pd.DataFrame(final_df)
#    df3=pd.DataFrame(final_solo)
    # input(df2.head())
    df2.to_csv(
        'data/input/processed/MSA_Dataset'+str(dataset)+'.csv')
    print('printed')


if __name__ == '__main__':  # TRY BOTH 5 AND 10 AA ON THE FLANKS
    dataset = 1
    MSA_df(dataset)
    dataset = 2
    MSA_df(dataset)
    dataset = 3
    MSA_df(dataset)
#    df3.to_csv('..\\clean scripts combined sets\\one_per_fam_comb_20_atom_C_train_MSA_uniprot_descriptors_5flank.csv')
