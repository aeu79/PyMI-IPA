import sys
import os
import scipy.io as spio
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import numpy as np


### glob parameters
replicate=1
Nincrement=1600
LengthA=64
msa_fasta_filename = 'Standard_HKRR_dataset.fasta'
mat=spio.loadmat('SpeciesNumbering_Standard_HKRR_dataset.mat',squeeze_me=True)

##from matlab data to df
col1 = list(range(len(mat['SpeciesNumbering_extr'])))
col2 = [x[1] for x in mat['SpeciesNumbering_extr']]
SpeciesNumbering_extr = pd.DataFrame([*zip(col1,col2)])


def main():
    encoded_focus_alignment, encoded_focus_alignment_headers, L = readAlignment_and_NumberSpecies(msa_fasta_filename, SpeciesNumbering_extr)





















def readAlignment_and_NumberSpecies(msa_fasta_filename,SpeciesNumbering_extr):
    '''
    This assumes that the first sequence is a reference sequence and that the last one is a dummy
    '''

    #read in alignment
    full_alignment = []
    encoded_focus_alignment_headers = []
    with open(msa_fasta_filename) as handle:
        for values in SimpleFastaParser(handle):
            encoded_focus_alignment_headers.append(values[0])
            full_alignment.append(values[1])

    alignment_height = len(full_alignment)
    alignment_width = len(full_alignment[0])

    # initialize
    encoded_focus_alignment = np.zeros((alignment_height, alignment_width+2))

    # encoded the sequence
    letter2number_map = {'-': 1, 'A': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8, 'I': 9, 'K': 10, 'L': 11,
                         'M': 12, 'N': 13, 'P': 14, 'Q': 15, 'R': 16, 'S': 17, 'T': 18, 'V': 19, 'W': 20, 'Y': 21,
                         'B': -1, 'Z': -1, 'J': -1, 'O': -1, 'U': -1, 'X': -1, 'a': -2, 'c': -2, 'd': -2, 'e': -2,
                         'f': -2, 'g': -2, 'h': -2, 'i': -2, 'k': -2, 'l': -2, 'm': -2, 'n': -2, 'p': -2, 'q': -2,
                         'r': -2, 's': -2,
                         't': -2, 'v': -2, 'w': -2, 'y': -2, 'b': -2, 'z': -2, 'j': -2, 'x': -2, 'u': -2, 'o': -2,
                         '.': -3}

    for i in range(alignment_height):
        encoded_focus_alignment[i, :-2] = np.array([letter2number_map.get(ch,0) for ch in full_alignment[i]])

        if i == 0 :
            encoded_focus_alignment[i, -2] = -1

        if 0 < i < alignment_height-1 :
            species_id = encoded_focus_alignment_headers[i].split('|')[1]
            #print(i,species_id)
            encoded_focus_alignment[i, -2] =SpeciesNumbering_extr.loc[SpeciesNumbering_extr[1]==species_id,0].iloc[0]

        if i < alignment_height -1:
            encoded_focus_alignment[i, -1] = i
        else:
            encoded_focus_alignment[i, -2] = float("NaN")
            encoded_focus_alignment[i, -1] = float("NaN")

    encoded_focus_alignment = np.delete(encoded_focus_alignment, np.where(~encoded_focus_alignment.any(axis=1)), axis=0)
    return encoded_focus_alignment, encoded_focus_alignment_headers,alignment_width








