import math
from itertools import product

import numpy as np
import pandas as pd
import scipy.io as spio
from Bio.SeqIO.FastaIO import SimpleFastaParser
from munkres import Munkres

import pmis

# setting
np.set_printoptions(suppress=True)

### glob parameters
replicate = 1
Nincrement = 6
LengthA = 64
msa_fasta_filename = 'Standard_HKRR_dataset.fasta'
mat = spio.loadmat('SpeciesNumbering_Standard_HKRR_dataset.mat', squeeze_me=True)

##from matlab data to df
col1 = list(range(len(mat['SpeciesNumbering_extr'])))
col2 = [x[1] for x in mat['SpeciesNumbering_extr']]
SpeciesNumbering_extr = pd.DataFrame([*zip(col1, col2)])


def main():
    # read sequences, adding species number in L+1 and sequence number in L+2
    # L is the full length of concatenated sequences, without supplementary indicators such as species and initial index.

    encoded_focus_alignment, encoded_focus_alignment_headers, L = readAlignment_and_NumberSpecies(msa_fasta_filename,
                                                                                                  SpeciesNumbering_extr)

    # suppress species with one pair
    encoded_focus_alignment, encoded_focus_alignment_headers = SuppressSpeciesWithOnePair(encoded_focus_alignment,
                                                                                          encoded_focus_alignment_headers)
    # tabulate species and sequences within species
    table_count_species = count_species(encoded_focus_alignment)

    # number of sequences
    N = encoded_focus_alignment.shape[0]
    print("Number of sequences: " + str(N))

    # number of runs(last one -> all sequences are in the training set)
    Nrounds = math.ceil(encoded_focus_alignment.shape[0] / Nincrement + 1)
    print("Number of rounds: " + str(Nrounds))

    # start from random within-speicies pairings: scrable the pairings for this.
    encoded_training_alignment = ScrambleSeqs(encoded_focus_alignment, LengthA, table_count_species)

    ##save the species and initial indices of the sequences in the scrambled alignment we start from
    filename = 'Res_python/IniScrambling_Ninc' + str(Nincrement) + '_rep' + str(replicate) + '.txt'
    np.savetxt(filename, encoded_training_alignment[:, L:], fmt='%d', delimiter='\t')
    encoded_training_alignment = np.delete(encoded_training_alignment, [L, L + 1, L + 2, L + 3], axis=1)
    #   np.save('encoded_focus_alignment.npy',encoded_focus_alignment)

    # initialize
    NSeqs_new = 0
    Output = np.zeros((Nrounds, 6))

    # iterate
    for rounds in range(Nrounds):

        print("Round: " + str(rounds + 1))

        if rounds > 0:
            # Use the gap to rank pairs
            Results = Results[Results[:, -1].argsort()[::-1]]

            # number of sequences that will be added to form the training set for this round
            NSeqs_new = NSeqs_new + Nincrement
            if NSeqs_new >= Results.shape[0]:
                NSeqs_new = Results.shape[0]  # for the last round, all paired sequences will be in the training set

            # save to Output the number of TP or FP in the training set
            Output[rounds, 4] = np.count_nonzero(Results[:NSeqs_new, 1] == Results[:NSeqs_new, 2])
            Output[rounds, 5] = np.count_nonzero(Results[:NSeqs_new, 1] != Results[:NSeqs_new, 2])

            # construct new training set
            newseqs = np.zeros((NSeqs_new, L))

            sorted_lst1 = np.searchsorted(encoded_focus_alignment[:, L + 1], Results[:NSeqs_new, 1])
            sorted_lst2 = np.searchsorted(encoded_focus_alignment[:, L + 1], Results[:NSeqs_new, 2])

            newseqs[:, :LengthA] = encoded_focus_alignment[sorted_lst1, :LengthA]
            newseqs[:, LengthA:L] = encoded_focus_alignment[sorted_lst2, LengthA:L]

            encoded_training_alignment = newseqs

        PMIs, Meff = pmis.compute_pmis(encoded_training_alignment, 0.15, 0.15)

        # compute pairings and gap scores for all pairs
        Results = Predict_pairs(encoded_focus_alignment, -PMIs, LengthA, table_count_species)
        Output[rounds, 0] = NSeqs_new
        Output[rounds, 1] = Meff
        Output[rounds, 2] = np.count_nonzero(Results[:, 1] == Results[:, 2])
        Output[rounds, 3] = np.count_nonzero(Results[:, 1] != Results[:, 2])

    filename = 'Res_python/TP_data_Ninc' + str(Nincrement) + '_rep' + str(replicate) + '.txt'
    np.savetxt(filename, Output, fmt=['%d', '%1.3f', '%d', '%d', '%d', '%d'], delimiter='\t')

    filename = 'Res_python/Resf_Ninc' + str(Nincrement) + '_rep' + str(replicate) + '.txt'
    np.savetxt(filename, Results, fmt=['%d', '%d', '%d', '%1.3f', '%1.3f'], delimiter='\t')



def readAlignment_and_NumberSpecies(msa_fasta_filename, SpeciesNumbering_extr):
    '''
    This assumes that the first sequence is a reference sequence and that the last one is a dummy
    '''

    # read in alignment
    full_alignment = []
    encoded_focus_alignment_headers = []
    with open(msa_fasta_filename) as handle:
        for values in SimpleFastaParser(handle):
            encoded_focus_alignment_headers.append(values[0])
            full_alignment.append(values[1])

    alignment_height = len(full_alignment)
    alignment_width = len(full_alignment[0])

    # initialize
    encoded_focus_alignment = np.zeros((alignment_height, alignment_width + 2))

    # encoded the sequence
    letter2number_map = {'-': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'K': 9, 'L': 10,
                         'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 20,
                         'B': -2, 'Z': -2, 'J': -2, 'O': -2, 'U': -2, 'X': -2, 'a': -3, 'c': -3, 'd': -3, 'e': -3,
                         'f': -3, 'g': -3, 'h': -3, 'i': -3, 'k': -3, 'l': -3, 'm': -3, 'n': -3, 'p': -3, 'q': -3,
                         'r': -3, 's': -3,
                         't': -3, 'v': -3, 'w': -3, 'y': -3, 'b': -3, 'z': -3, 'j': -3, 'x': -3, 'u': -3, 'o': -3,
                         '.': -4}

    for i in range(alignment_height):
        encoded_focus_alignment[i, :-2] = np.array([letter2number_map.get(ch, 0) for ch in full_alignment[i]])

        if i == 0:
            encoded_focus_alignment[i, -2] = -1

        if 0 < i < alignment_height - 1:
            species_id = encoded_focus_alignment_headers[i].split('|')[1]
            # print(i,species_id)
            encoded_focus_alignment[i, -2] = SpeciesNumbering_extr.loc[SpeciesNumbering_extr[1] == species_id, 0].iloc[
                0]

        if i < alignment_height - 1:
            encoded_focus_alignment[i, -1] = i
        else:
            encoded_focus_alignment[i, -2] = float("NaN")
            encoded_focus_alignment[i, -1] = float("NaN")

    encoded_focus_alignment = np.delete(encoded_focus_alignment, np.where(~encoded_focus_alignment.any(axis=1)), axis=0)
    return encoded_focus_alignment, encoded_focus_alignment_headers, alignment_width


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
    # print(table_count_species.shape)
    species_id = encoded_focus_alignment[1, L]
    iini = 1
    count = -1
    for i in range(1, N - 1):
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

    colval, count = np.unique(encoded_focus_alignment[:, -2], return_counts=True)
    colval = colval.astype(int)

    # col stores speciesID which has only one protein
    col = []
    for speciesID, occu in dict(zip(colval, count)).items():
        if occu == 1:
            col.append(speciesID)

    # delete the first sequence  and last dummy...?(she defined it!)
    del col[0]
    del col[-1]

    NUP_alignment = encoded_focus_alignment
    mask = ~np.isin(NUP_alignment[:, -2], col)
    NUP_alignment = NUP_alignment[mask, :]
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
    L = alignment_width - 2  # last 2 cols contain species index and initial sequence index!

    # added columns will contain the same data but for the RR
    scrambled_alignment = np.zeros((N, L + 4))

    # loop over species
    for i in range(len(table_count_species)):
        # number of sequences in this species
        Nseqs = int(table_count_species[i, 2] - table_count_species[i, 1] + 1)
        # random permutation
        thePerm = np.arange(Nseqs)
        np.random.shuffle(thePerm)

        for j in range(Nseqs):
            HKIndex = int(table_count_species[i, 1] + j)
            RRIndex = int(table_count_species[i, 1] + thePerm[j])
            scrambled_alignment[HKIndex, :LengthA] = encoded_focus_alignment[HKIndex, :LengthA]
            scrambled_alignment[HKIndex, LengthA:L] = encoded_focus_alignment[RRIndex, LengthA:L]
            scrambled_alignment[HKIndex, L:L + 2] = encoded_focus_alignment[HKIndex, L:L + 2]
            scrambled_alignment[HKIndex, L + 2:L + 4] = encoded_focus_alignment[RRIndex, L:L + 2]

    scrambled_alignment = np.delete(scrambled_alignment, np.where(~scrambled_alignment.any(axis=1)), axis=0)
    return scrambled_alignment


def Predict_pairs(encoded_focus_alignment, PMIs, LengthA, table_count_species):
    '''
    makes pairing predictions on encoded_focus_alignment using PMI scores.
    :param encoded_focus_alignment:
    :param PMIs:
    :param LengthA: length of the first protein in the pair.
    :param table_count_species:
    :return:
    '''
    encoded_focus_alignment = encoded_focus_alignment.astype(int)
    table_count_species = table_count_species.astype(int)
    N, alignment_width = encoded_focus_alignment.shape
    L = alignment_width - 2

    # initialize the Results array, used for saving data
    # col 1: species
    # col 2: HK index in initial alignment
    # col 3: RR index in initial alignment
    # col 4: score of pairing
    # col 5: gap
    Results = np.zeros((5, N - 2)).astype(object)  # need to transform to get the same format as matlab version.
    # total pair counter
    pre_tol = 0
    cur_tol = 0

    for i in range(table_count_species.shape[0]):
        test_seqs = encoded_focus_alignment[table_count_species[i, 1]:table_count_species[i, 2] + 1, :]
        Nseqs = table_count_species[i, 2] - table_count_species[i, 1] + 1
        species_id = table_count_species[i, 0]

        # now compute the PMI-based pairing score of all the HK-RR pairs within the species corresponding to i
        Pairing_scores = Compute_pairing_scores(test_seqs, Nseqs, PMIs, LengthA, L)

        if Nseqs == 1:
            assignment = 1
            Pairing_scores_b = Pairing_scores - Pairing_scores.min()  # ensure that all elements are >= 0

        elif Pairing_scores.min() == Pairing_scores.max():

            # random permutation
            thePerm = np.arange(Nseqs)
            np.random.shuffle(thePerm)  # avoid spurious positive results
            Pairing_scores_b = Pairing_scores - Pairing_scores.min()

        else:
            # use the Hungarian algorithm
            Pairing_scores_b = Pairing_scores - Pairing_scores.min()
            assignment, score = assignmentoptimal(Pairing_scores_b.tolist())  # input requires list not np.array
            # deal with identical rows, i.e.effectively identical HKs
            uEn = np.unique(Pairing_scores_b, axis=0)
            if uEn.shape[0] != Pairing_scores_b.shape[0]:
                assignment = randomize_equal(assignment, Pairing_scores_b, uEn, 'row')

            uEn = np.unique(Pairing_scores_b, axis=1)
            if uEn.shape[1] != Pairing_scores_b.shape[1]:
                assignment = randomize_equal(assignment, Pairing_scores_b, uEn, 'col')

        bigval = 1e3 * Pairing_scores_b.max()

        cur_tol = pre_tol + Nseqs
        # print(i)
        Results[:-1, pre_tol:cur_tol] = np.array(
            [np.repeat(species_id, Nseqs), test_seqs[:, L + 1], test_seqs[assignment, L + 1],
             Pairing_scores[np.arange(Nseqs), assignment]])

        if Nseqs == 1:
            Results[4, pre_tol:cur_tol] = abs(Pairing_scores)

        elif Pairing_scores.min() == Pairing_scores.max():
            Results[4, pre_tol:cur_tol] = 0

        else:
            for j in range(Nseqs):
                Pairing_scores_mod = np.copy(Pairing_scores_b)
                Pairing_scores_mod[j, assignment[j]] = bigval
                useless, score_mod = assignmentoptimal(Pairing_scores_mod.tolist())
                Results[4, pre_tol + j] = score_mod - score

        pre_tol = cur_tol

    return Results.T


def Compute_pairing_scores(test_seqs, Nseqs, PMIs, LengthA, L):
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
    for k, l in list(product(list(range(LengthA)), list(range(LengthA, L)))):
        indlst.append([k, l])

    ## length of a_lst/b_lst : 64*112=7168
    a_lst = [x[0] for x in indlst]  # element in a_lst: position in protein A (from 0 to 63)
    b_lst = [x[1] for x in indlst]  # element in b_lst: position in protein B (from 64 to 175)

    ## fill in Pairing_scores matrix
    for index in np.ndindex(Pairing_scores.shape):
        aa1_lst = test_seqs[index[0], a_lst]
        aa2_lst = test_seqs[index[1], b_lst]

        # ## -1 : aa1_lst.max is 21.0 and need to be 20 and integer.
        # aa1_lst = [int(x)-1 for x in aa1_lst]
        # aa2_lst = [int(x)-1 for x in aa2_lst]

        Pairing_scores[index] = PMIs[a_lst, b_lst, aa1_lst, aa2_lst].sum()

    return Pairing_scores


def assignmentoptimal(distMatrix):
    '''
    Compute optimal assignment by Munkres algorithm
    ASSIGNMENTOPTIMAL(DISTMATRIX) computes the optimal assignment (minimum
    overall costs) for the given rectangular distance or cost matrix, for
    example the assignment of tracks (in rows) to observations (in
    columns). The result is a column vector containing the assigned column
    number in each row (or 0 if no assignment could be done).

    [ASSIGNMENT, COST] = ASSIGNMENTOPTIMAL(DISTMATRIX) returns the
    assignment vector and the overall cost.

    Here we use a python package for Hungarian Algorithm(https://github.com/bmc/munkres).

    :param distMatrix: a cost matrix but in list format(Note: this package does not support numpy array!))
    :return:assignment vector(np.array), minimum cost.
    '''

    m = Munkres()
    indexes = m.compute(distMatrix)
    assignment = np.array([x[1] for x in indexes])
    cost = 0
    for row, column in indexes:
        value = distMatrix[row][column]
        cost += value

    return assignment, cost


def randomize_equal(assignment, HKRR_energy_b, uEn, vec):
    '''
    deal with identical columns, i.e. effectively identical RRs
    :param assignment: assignment vector from Hungarian algorithm
    :param HKRR_energy_b: Pairing_scores_b where all elements are positive.
    :param uEn: Pairing_scores_b after applying unique rows or columns
    :param vec: argument of rows or columns.
    :return: a random permutation of these repeated vectors in assignment
    '''

    if vec == 'col':
        uEn = uEn.T
        HKRR_energy_b = HKRR_energy_b.T

    for i in range(uEn.shape[0]):
        rows = np.where((uEn[i] == HKRR_energy_b).all(axis=1))
        if len(rows[0]) > 1:
            thePerm = np.arange(len(rows[0]))
            np.random.shuffle(thePerm)
            #  print(thePerm)
            assignment[rows] = assignment[rows[0][thePerm]]

    return assignment


if __name__ == "__main__":
    main()