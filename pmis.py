import math
from itertools import product

import numpy as np
import scipy.spatial as scispa


def compute_pmis(encoded_focus_alignment, pseudocount_weight, theta):
    encoded_focus_alignment = encoded_focus_alignment.astype(int)

    q = 21
    # compute the empirical frequencies from the sequence data
    meff, Pij_true, Pi_true, alignment_width = count_alignment(encoded_focus_alignment, theta, q)

    # include pseudocounts
    Pij, Pi = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q)
    # compute PMIs
    PMIs = Get_PMIs(Pij, Pi, alignment_width, q)
    # np.save('PMIsssfile.npy',PMIs)
    return PMIs, meff


##auxiliary function definitions
def count_alignment(encoded_focus_alignment, theta, q):
    # calculate Meff
    alignment_height = encoded_focus_alignment.shape[0]
    alignment_width = encoded_focus_alignment.shape[1]
    W = np.ones((1, alignment_height))
    # whether you weight or not
    if theta > 0.0:
        encoded_focus_alignment_pdist = scispa.distance.pdist(encoded_focus_alignment, 'hamming')
        booleanTF = theta > encoded_focus_alignment_pdist
        boolean01_encoded_focus = booleanTF.astype(int)
        W = (1. / (1 + sum(scispa.distance.squareform(boolean01_encoded_focus))))
    Meff = sum(W)
    # compute the frequencies
    Pij_true = np.zeros((alignment_width, alignment_width, q, q))
    Pi_true = np.zeros((alignment_width, q))

    # single-site frequencies

    def vectorize(data):
        result = np.zeros((q,))
        for i in range(q):
            result[i] = W[data == i].sum()
        return result

    Pi_true = np.apply_along_axis(vectorize, 0, encoded_focus_alignment)
    Pi_true = Pi_true.T

    # normalization
    Pi_true = Pi_true / Meff

    # two-site frequencies

    Pij_ind = np.triu_indices(alignment_width, 1)

    for k in range(len(Pij_ind[0])):
        e, j = Pij_ind[0][k], Pij_ind[1][k]
        e_col = encoded_focus_alignment[:, e]
        j_col = encoded_focus_alignment[:, j]
        np.add.at(Pij_true[e, j, :, :], [e_col.tolist(), j_col.tolist()], W)
        Pij_true[j, e, :, :] = Pij_true[e, j, :, :].T

    Pij_true = Pij_true / Meff

    # fill the diagonal of Pij_true
    diag_ind = product(list(range(alignment_width)), list(range(q)))
    for k in diag_ind:
        Pij_true[k[0], k[0], k[1], k[1]] = Pi_true[k[0], k[1]]

    return (Meff, Pij_true, Pi_true, alignment_width)


def with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q):
    # add pseudocounts to deal with some finite size effects
    Pij = (1. - pseudocount_weight) * Pij_true + pseudocount_weight / q / q * (
        np.ones((alignment_width, alignment_width, q, q)))
    Pi = (1. - pseudocount_weight) * Pi_true + pseudocount_weight / q * (np.ones((alignment_width, q)))
    # correct things on the diagonal
    # scra = np.identity(q)
    # for i in range(alignment_width):
    #     for alpha in range(1,q):
    #         for beta in range(1,q):
    #             Pij[i, i, alpha, beta] = (1.-pseudocount_weight)*Pij_true[i, i, alpha, beta] + pseudocount_weight/q*scra[alpha, beta]
    #
    psedo_ind = product(list(range(alignment_width)), list(range(1, q)))
    for k in psedo_ind:
        Pij[k[0], k[0], k[1], k[1]] = (1. - pseudocount_weight) * Pij_true[
            k[0], k[0], k[1], k[1]] + pseudocount_weight / q

    ##keep non-diagnal elements zero
    for l in range(alignment_width):
        Pij[l, l, 1:, 1:][np.triu_indices(q - 1, 1)] = 0
        Pij[l, l, 1:, 1:][np.triu_indices(q - 1, 1)[::-1]] = 0

    return (Pij, Pi)


def Get_PMIs(Pij, Pi, alignment_width, q):
    # Compute PMIs from the matrices of frequencies.
    # Here, the 1st aa type(=gap) is ignored (slight improvement)
    C = np.zeros((alignment_width, alignment_width, q, q))

    ## my code
    C_ind = np.triu_indices(alignment_width, 0)
    np.seterr(divide='ignore') # Avoid warning log(0)
    for m in range(len(C_ind[0])):
        i, j = C_ind[0][m], C_ind[1][m]
        Pi_mat = np.dot(Pi[i, :][1:, None], Pi[j, :][None, 1:])
        C[i, j, 1:, 1:] = np.log(Pij[i, j, 1:, 1:] / Pi_mat)
        # Alternative to avoid the devide by zero warning: TODO
        # resulti = np.where(Pij[i, j, 1:, 1:] / Pi_mat > 0.0000000001, Pij[i, j, 1:, 1:] / Pi_mat, -math.inf)
        # C[i, j, 1:, 1:] = np.log(resulti, out=resulti, where=resulti > 0)
        C[j, i, 1:, 1:] = np.log(Pij[j, i, 1:, 1:] / Pi_mat.T)
        #resultj = np.where(Pij[j, i, 1:, 1:] / Pi_mat.T > 0.0000000001, Pij[j, i, 1:, 1:] / Pi_mat.T, -math.inf)
        #C[j, i, 1:, 1:] = np.log(resultj, out=resultj, where=resultj > 0)
    np.seterr(divide='warn') # warnings back on

    # mask for Pij=0
    C[Pij == 0] = -math.inf

    # for i in range(alignment_width):
    #     for j in range(alignment_width):
    #         for alpha in range(1,q):
    #             for beta in range(1,q):
    #                 if Pij[i, j, alpha, beta]!=0:
    #                     C[i,j,alpha,beta]= math.log((Pij[i, j, alpha, beta]) / (Pi[i, alpha]*Pi[j, beta]))
    #                 else:
    #                     C[i,j,alpha,beta]=-math.inf

    return C


if __name__ == "__main__":
    fi = np.load('encoded_training_alignment.npy')
    # fi= os.sys.argv[1]
    pseudocount_weight = 0.15
    theta = 0.15
    compute_pmis(fi, pseudocount_weight, theta)
