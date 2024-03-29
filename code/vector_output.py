from Bio.PDB import *
from Bio import pairwise2
from Bio.SVDSuperimposer import SVDSuperimposer
from itertools import groupby
from operator import itemgetter
import os
import numpy as np
import warnings
import pandas as pd
from tqdm import tqdm
warnings.filterwarnings('ignore')

VALID_ATOMS = ('C', 'CA', 'N')


def get_protein_list(input_dir='data/input/raw/Dataset1'):
    ls = []
    for name in os.listdir(input_dir):
        if 'alpha' in name:
            ls.append(name[:4])
    return(ls)


def get_baker_coordinates(protein, dataset):
    atom_coordinates = []
    p = PDBParser()
    structure = p.get_structure(
        '', f'data/input/raw/Dataset{dataset}/{protein}_rosetta.pdb')
    for atom in structure.get_atoms():
        atom_coordinates.append((atom.get_name(), list(atom.get_coord())))

    return atom_coordinates


def get_alpha_ppsequence_and_coordinates(protein, dataset):
    atom_coordinates = []
    p = PDBParser()
    structure = p.get_structure(
        '', f'data/input/raw/Dataset{dataset}/{protein}_alpha.pdb')
    for atom in structure.get_atoms():
        atom_name = atom.get_name()
        if atom_name in VALID_ATOMS:
            atom_coordinates.append((atom_name, list(atom.get_coord())))

    ppb = PPBuilder()
    pp_sequence = ''
    for polypeptide in ppb.build_peptides(structure):
        pp_sequence = str(polypeptide.get_sequence())

    return atom_coordinates, pp_sequence


def get_backbone_indices(sequence, overlapping_indices):
    backbone_indices = []
    oi_set = set(overlapping_indices)
    current_end = 0
    for atom_index in range(len(sequence)):
        current_end += 3
        if atom_index in oi_set:
            for i in range(current_end-3, current_end):
                backbone_indices.append(i)

    return backbone_indices


def align_alpha_reference_seqs(protein, alpha_polypeptide_sequence, dataset):
    p = MMCIFParser()
    structure = p.get_structure(
        '', f'data/input/raw/ref_Dataset{dataset}/{protein.lower()}.cif')

    all_chains = []
    done = []
    for chains in structure:
        for chain1 in chains:
            all_chains.append([chain1, len(chain1)])
            ppb = PPBuilder()
            # print(ppb.build_peptides(chain1)[0].get_sequence())
    # print(all_chains)
    if 1 == 1:
        # sort all available chains according to their length allowing us to pick out the longest possible match
        all_chains.sort(key=lambda x: x[1], reverse=True)
        # input(all_chains)
        all_chains = [c[0] for c in all_chains]
        for chain in all_chains:
            atom_coordinates = []
            for atom in chain.get_atoms():
                atom_name = atom.get_name()
                if atom_name in VALID_ATOMS:
                    atom_coordinates.append(
                        (atom_name, list(atom.get_coord())))

            ppb = PPBuilder()
            for polypeptide in ppb.build_peptides(chain):
                pp_sequence = str(polypeptide.get_sequence())
                # input(pp_sequence)

                # int(len(alpha_polypeptide_sequence)/2): ensuring that we're at minimum 50 (0.8*62)
                if int(len(pp_sequence)) < 63:
                    #print('too small', protein, int(len(pp_sequence)), int(len(alpha_polypeptide_sequence)))
                    continue
                else:
                    # input(pp_sequence)
                    # Align reference and AlphaFold polypeptide sequences.
                    alignment = pairwise2.align.globalxx(
                        alpha_polypeptide_sequence, pp_sequence)
                    # input(alignment)
                    aligned_alpha_sequence = alignment[0][0]
                    aligned_ref_sequence = alignment[0][1]
    # if protein =='1B6U':
    ##                    print(len(aligned_alpha_sequence), len(aligned_ref_sequence))
                    assert(len(aligned_alpha_sequence)
                           == len(aligned_ref_sequence))

                    # @NOTE: The overlapping indices here are the ones where the
                    # amino acid in the AlphaFold sequence matches that in the
                    # reference, because it is the reference that is moved.
                    index = 0
                    overlapping_indices = []
                    for i, j in zip(aligned_alpha_sequence, aligned_ref_sequence):
                        if i == j:
                            overlapping_indices.append(index)
                        index += 1

                    # Group overlapping indices by whether there's a gap in the
                    # alignment.
                    grouped_overlapping_indices = []
                    for h, g in groupby(enumerate(overlapping_indices), lambda x: x[0] - x[1]):
                        grouped_overlapping_indices.append(
                            list(map(itemgetter(1), g)))
                    # Choose longest continuous sequence if there are multiple
                    # sequences after grouping them.
                    longest_overlapping_indices = overlapping_indices

                    if len(grouped_overlapping_indices) > 1:
                        longest_seen = 0
                        longest_index = 0

                        for i, subset_overlapping_indices in enumerate(grouped_overlapping_indices):
                            soi_len = len(subset_overlapping_indices)

                            if soi_len > longest_seen:
                                longest_seen = soi_len
                                longest_index = i
                        longest_overlapping_indices = grouped_overlapping_indices[longest_index]

                    # Get sequence of longest overlapping indices.
                    longest_overlapping_sequence_list = []
                    for i in set(longest_overlapping_indices):
                        longest_overlapping_sequence_list.append(
                            aligned_ref_sequence[i])
                    longest_overlapping_sequence = ''.join(
                        longest_overlapping_sequence_list)
                    # Get longest non-overlapping indices.
                    ref_overlapping_indices = []
                    for starting_index in range(len(pp_sequence)):

                        pp_sequence_subset = pp_sequence[starting_index:len(
                            longest_overlapping_sequence)]
    # if protein =='1B6U':
    # print(pp_sequence_subset)
    # input(longest_overlapping_sequence)
                        if pp_sequence_subset == longest_overlapping_sequence:

                            los_len = len(longest_overlapping_sequence)

                            ref_overlapping_indices = list(
                                range(starting_index, los_len))
    # if protein =='1B6U':
    # print(ref_overlapping_indices)
                    # Get predicted atoms covered by the alignments at the 80%
                    # coverage threshold.
                    alignment2 = int(pairwise2.align.globalxx(
                        pp_sequence, alpha_polypeptide_sequence)[0][2])
                    longest_aligned = alignment2
                    # input(longest_aligned)
                    if longest_aligned > (0.8 * len(pp_sequence)):

                        ref_backbone_indices = get_backbone_indices(pp_sequence,
                                                                    ref_overlapping_indices)
                        alpha_backbone_indices = get_backbone_indices(
                            alpha_polypeptide_sequence,
                            longest_overlapping_indices)

                        # if protein =='1B6U':
                        #print(ref_backbone_indices, alpha_backbone_indices)
                        # print(alpha_backbone_indices)
                        # input(ref_backbone_indices)
                        return (atom_coordinates, ref_backbone_indices, alpha_backbone_indices)

    return None


def get_transformed_coordinates(primary_array, secondary_array):
    sup = SVDSuperimposer()
    sup.set(primary_array, secondary_array)
    sup.run()
    #rms = sup.get_rms()
    #rms = '{:.3f}'.format(rms)
    return sup.get_transformed()


def calculate_distances(alpha_trans, baker_trans, ref_array):
    distances = []
    input_vectors = []

    for i in range(len(ref_array)):
        alpha_p = alpha_trans[i]
        baker_p = baker_trans[i]
        ref_p = ref_array[i]

        dir_alpha_to_baker = baker_p - alpha_p

        # direction1 *direction2 = 0 # due to 90 degree angle at intersection
        #direction1 = dir_alpha_to_baker

        # direction2 = (alpha coordinate + t*direction1)-reference_coordinate

        # as stated above: #direction1 *direction2 = 0 # due to 90 degree angle at intersection

        # therefore: dir_alpha_to_baker*((alpha coordinate + t*direction1)-reference_coordinate))=0

        # multiply in direction1 and solve for t

        t = -1*np.sum((alpha_p-ref_p)*dir_alpha_to_baker) / \
            np.sum(dir_alpha_to_baker**2)

        # wrong: t = -1 * (dir_alpha_to_baker*alpha_p-dir_alpha_to_baker*ref_p)/dir_alpha_to_baker**2#(np.sum(mult_dir))/(np.sum(mult_dir_exp))

        # input(t)
        complete_point = alpha_p + dir_alpha_to_baker * t

        distance = complete_point - alpha_p
        # @NOTE: Not a fan of this re-assignment her, but it is what it is. -Jk
        distance = np.sqrt(distance[0]**2 + distance[1]**2 + distance[2]**2)

        # Same as dir_alpha_from_baker. Kept for continuity, will update later.
        ref1 = baker_p - alpha_p
        # Direction of the complete point to the baker coordinate.
        ref2 = baker_p - complete_point

        ref_dist_1 = np.sqrt(ref1[0]**2 + ref1[1]**2 + ref1[2]**2)
        ref_dist_2 = np.sqrt(ref2[0]**2 + ref2[1]**2 + ref2[2]**2)

        if ref_dist_2 > ref_dist_1:
            if distance < ref_dist_2:
                distance = distance * -1

        input_vector = list(ref_p) + list(alpha_p) + list(ref1)

        distances.append(distance)
        input_vectors.append(input_vector)

    return distances, input_vectors


def superimpose(alpha_atom_coordinates, baker_atom_coordinates, ref_atom_coordinates):
    alpha_array = np.array([i[1] for i in alpha_atom_coordinates])
    baker_array = np.array([i[1] for i in baker_atom_coordinates])
    ref_array = np.array([i[1] for i in ref_atom_coordinates])

    alpha_trans = get_transformed_coordinates(ref_array, alpha_array)
    baker_trans = get_transformed_coordinates(alpha_trans, baker_array)

    return calculate_distances(alpha_trans, baker_trans, ref_array)


def main(protein, dataset):
    #protein = '1A1B'

    base_baker_atom_coordinates = get_baker_coordinates(protein, dataset)

    alpha_prediction = get_alpha_ppsequence_and_coordinates(protein, dataset)
    base_alpha_atom_coordinates = alpha_prediction[0]
    alpha_polypeptide_sequence = alpha_prediction[1]

    alignment_outcome = align_alpha_reference_seqs(
        protein, alpha_polypeptide_sequence, dataset)
   # print(alignment_outcome[2])
    ref_atom_coordinates = alignment_outcome[0]
    ref_backbone_indices = alignment_outcome[1]
    alpha_backbone_indices = alignment_outcome[2]

    ax = alpha_backbone_indices[0]
    ay = alpha_backbone_indices[-1]+1
    alpha_atom_coordinates = base_alpha_atom_coordinates[ax:ay]
    baker_atom_coordinates = base_baker_atom_coordinates[ax:ay]

    rx = ref_backbone_indices[0]
    ry = ref_backbone_indices[-1]+1
    ref_atom_coordinates = ref_atom_coordinates[rx:ry]

    la = len(alpha_atom_coordinates)
    lb = len(baker_atom_coordinates)
    lr = len(ref_atom_coordinates)
    assert(la == lb == lr)

    #print(f'Base Baker Coords:     {len(base_baker_atom_coordinates)}')
    #print(f'Base AlphaFold Coords: {len(base_alpha_atom_coordinates)}')
    #print(f'Base Ref Coords:       {len(ref_atom_coordinates)}')
    # print()
    #print(f'Ref backbone:          {len(ref_backbone_indices)}')
    #print(f'Alpha backbone:        {len(alpha_backbone_indices)}')
    # print()
    #print(f'Base Baker Coords:     {len(baker_atom_coordinates)}')
    #print(f'Base AlphaFold Coords: {len(alpha_atom_coordinates)}')
    # print()
    #print(f'AlphaFold Seq:         {len(alpha_polypeptide_sequence)}')

    distances, vector_info = superimpose(
        alpha_atom_coordinates, baker_atom_coordinates, ref_atom_coordinates)
    assert(len(distances) == len(vector_info))
    return distances, vector_info, alpha_backbone_indices


def create_vector_info(dataset: int = 1):
    input_dir = 'data/input/raw/Dataset' + \
        str(dataset)+'/'
    prot_ls = get_protein_list(input_dir)
    vector_info = []
    atom_counts = []
    total = []
    numbers = np.arange(0, 1500, 100)
    numbers = dict(zip(numbers, numbers))
    tp = 0
    #tp = 1
    for protein in tqdm(prot_ls):
        # if protein == '4B3W':
        # input(protein)

        # print(protein)
        # if tp in numbers:
        #   print(tp)
        #tp += 1
        #protein = '4B3W'
        try:

            ls, vector_data, alpha_backbone_indices = main(protein, dataset)
    # total.append([prot]+ls)
    # atom_counts.append([prot]+rel_atm_alpha)
            vector_info.append([protein]+vector_data)
            atom_counts.append([protein]+alpha_backbone_indices)
            total.append([protein]+ls)
        except Exception as e:
            print(protein, e)
            # print(traceback.print_exc())
    df = pd.DataFrame(vector_info)
    df1 = pd.DataFrame(atom_counts)
    df2 = pd.DataFrame(total)
    df.to_csv('data/input/processed/vector_Dataset'+str(dataset)+'.csv')
    df1.to_csv('data/input/processed/atom_counts_Dataset'+str(dataset)+'.csv')
    df2.to_csv('data/input/processed/targets_Dataset'+str(dataset)+'.csv')


if __name__ == '__main__':
    dataset = 1
    create_vector_info(dataset=dataset)
    dataset = 2
    create_vector_info(dataset=dataset)
    dataset = 3
    create_vector_info(dataset=dataset)
