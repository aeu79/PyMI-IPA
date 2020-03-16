import sys
import os
import scipy.io as spio
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import numpy as np
import math


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
    
    # read sequences, adding species number in L+1 and sequence number in L+2
    # L is the full length of concatenated sequences, without supplementary indicators such as species and initial index.

    encoded_focus_alignment, encoded_focus_alignment_headers, L = readAlignment_and_NumberSpecies(msa_fasta_filename, SpeciesNumbering_extr)
    
    # suppress species with one pair
    encoded_focus_alignment, encoded_focus_alignment_headers = SuppressSpeciesWithOnePair(encoded_focus_alignment, encoded_focus_alignment_headers)
    # tabulate species and sequences within species
    table_count_species = count_species(encoded_focus_alignment)

    #number of sequences
    N = encoded_focus_alignment.shape[0]
#    print(N)

    # number of rouns(last one -> all sequences are in the training set)
    Nrounds = math.ceil(encoded_focus_alignment.shape[0]/Nincrement + 1)
#    print(Nrounds)

    #start from random within-speicies pairings: scrable the pairings for this. 
    encoded_training_alignment = ScrambleSeqs(encoded_focus_alignment, LengthA, table_count_species)
    encoded_training_alignment = np.delete(encoded_training_alignment,[ L , L + 1 , L + 2 , L + 3 ],axis = 1)
    print(encoded_training_alignment.shape)
 #   np.save('encoded_focus_alignment.npy',encoded_focus_alignment)
















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

def count_species(encoded_focus_alignment):
    '''
    return a table with species id, 1st index in encoded_focus_alignment, last index in encoded_focus_alignment
    :param
    encoded_focus_alignment: multiple sequence alignment
    :return: a table contains index start & end points for each species.
    '''

    (N, alignment_width) = encoded_focus_alignment.shape
    L = alignment_width - 2

    table_count_species = np.zeros((N, 3))
    #print(table_count_species.shape)
    species_id = encoded_focus_alignment[1, L]
    iini = 1
    count = -1
    for i in range(1, N-1):
        species_id_prev = species_id
        species_id = encoded_focus_alignment[i, L]

        if species_id != species_id_prev:
            count = count + 1
            ifin = i - 1

            table_count_species[count, 0] = species_id_prev
            table_count_species[count, 1] = iini
            table_count_species[count, 2] = ifin
            iini = i

    count = count + 1
    ifin = N - 1
    table_count_species[count, 0] = species_id
    table_count_species[count, 1] = iini
    table_count_species[count, 2] = ifin - 1

    table_count_species = np.delete(table_count_species, np.where(~table_count_species.any(axis=1)), axis=0)

    return table_count_species


def SuppressSpeciesWithOnePair(encoded_focus_alignment, encoded_focus_alignment_headers):
    '''
    produce alignment where the species with only one pair are suppressed.
    :param encoded_focus_alignment: multiple sequence alignment
    :param encoded_focus_alignment_headers: sequence header (species info)
    :return: sequence alignments & header info after removing species with only one protein.
    '''

    # according to readAlignment_and_NumberSpecies function, the second last column is for speicies ID.
    # colval: the value in the species ID
    # count : how many proteins in this species
    # to do : remove species within one protein from the alignment

    colval , count = np.unique(encoded_focus_alignment[:, -2], return_counts=True)
    colval = colval.astype(int)

    #col stores speciesID which has only one protein
    col=[]
    for speciesID, occu in dict(zip(colval, count)).items():
        if occu == 1:
            col.append(speciesID)

    #delete the first sequence  and last dummy...?(she defined it!)
    del col[0]
    del col[-1]

    NUP_alignment = encoded_focus_alignment
    mask = ~np.isin(NUP_alignment[:,-2],col)
    NUP_alignment = NUP_alignment[mask,:]
    NUP_alignment_headers = np.array(encoded_focus_alignment_headers)[mask].tolist()

    return NUP_alignment, NUP_alignment_headers

def ScrambleSeqs(encoded_focus_alignment, LengthA, table_count_species):
    '''
    scramble the pairs.
    :param encoded_focus_alignment: multiple sequence alignment
    :param LengthA: length of first protein of the pair.
    :param table_count_species: tabulate species and sequences within species
    :return: scrambled multiple sequence alignment.
    '''
    N, alignment_width = encoded_focus_alignment.shape
    L = alignment_width - 2   #last 2 cols contain species index and initial sequence index!

    # added columns will contain the same data but for the RR
    scrambled_alignment = np.zeros((N,L+4))

    # loop over species
    for i in range(len(table_count_species)):
        # number of sequences in this species
        Nseqs = int(table_count_species[i,2] -table_count_species[i,1] + 1)
        #random permutation
        thePerm = np.arange(Nseqs)
        np.random.shuffle(thePerm)

        for j in range(Nseqs):

            HKIndex = int(table_count_species[i, 1] + j)
            RRIndex = int(table_count_species[i, 1] + thePerm[j])
            scrambled_alignment[HKIndex, :LengthA] = encoded_focus_alignment[HKIndex, :LengthA]
            scrambled_alignment[HKIndex, LengthA :L] = encoded_focus_alignment[RRIndex, LengthA :L]
            scrambled_alignment[HKIndex, L :L + 2] = encoded_focus_alignment[HKIndex, L :L + 2]
            scrambled_alignment[HKIndex, L + 2:L + 4] = encoded_focus_alignment[RRIndex, L :L + 2]

    scrambled_alignment = np.delete(scrambled_alignment, np.where(~scrambled_alignment.any(axis=1)), axis=0)
    return scrambled_alignment




def Compute_pairing_scores(test_seqs,Nseqs,PMIs, LengthA, L):
    '''
    calculate PMI-based pairing scores between all pairs of HKs and RRs in test_seqs
    Line: HK; Col: RR.
    :param test_seqs: seqs within one species.
    :param NSeqs: number of sequences in this species.
    :param PMIs: previously-generated
    :param LengthA: length of 1st protein in the pair.
    :param L: length of two proteins.
    :return: pairing score.
    '''

    Pairing_scores = np.zeros((Nseqs, Nseqs))

    ## index coombination of two individual proteins: a[0,63] b[64,175]
    ## k,l appear in pairs: [0,64],[0,65],[0,66]...[63,175]
    indlst = []
    for k,l in list(product(list(range(LengthA)),list(range(LengthA,L)))):
        indlst.append([k,l])

    ## length of a_lst/b_lst : 64*112=7168
    a_lst = [x[0] for x in indlst]   # element in a_lst: position in protein A (from 0 to 63)
    b_lst = [x[1] for x in indlst]   # element in b_lst: position in protein B (from 64 to 175)


    ## fill in Pairing_scores matrix
    for index in np.ndindex(Pairing_scores.shape):
        aa1_lst = test_seqs[index[0],a_lst]
        aa2_lst = test_seqs[index[1],b_lst]

        ## -1 : aa1_lst.max is 21.0 and need to be 20 and integer.
        aa1_lst = map(lambda x: int(x)- 1 , aa1_lst)
        aa2_lst = map(lambda x: int(x)- 1 , aa2_lst)

        Pairing_scores[index] = PMIs[a_lst, b_lst, aa1_lst, aa2_lst].sum()

    return Pairing_scores



if __name__ == "__main__":
    main()




