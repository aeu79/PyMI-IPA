import numpy as np
import scipy.spatial as scispa
import os,sys, math
def Compute_PMIs(encoded_focus_alignment, pseudocount_weight, theta):
    encoded_focus_alignment=np.load('encoded_focus_alignment.npy')
    print(np.shape(encoded_focus_alignment))
    q= 21
#compute the empirical frequencies from the sequence data
    Meff, Pij_true, Pi_true, alignment_width = count_alignment(encoded_focus_alignment, theta, q)
#include pseudocounts
    Pij,Pi= with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q)
#compute PMIs
    PMIs = Get_PMIs(Pij, Pi, alignment_width, q)
    np.save('PMIs_file.npy',PMIs)
    #print(PMIs,Meff)
    return(PMIs, Meff)


##auxiliary function definitions
def count_alignment(encoded_focus_alignment, theta, q):
#calculate Meff
#    encoded_focus_alignment=np.load('encoded_focus_alignment.npy')
    alignment_height = encoded_focus_alignment.shape[0]
    print(alignment_height)
    alignment_width = encoded_focus_alignment.shape[1]
    W = np.ones((1, alignment_height))
    #print(np.shape(W))
    if theta > 0.0: # whether you weight or not
        encoded_focus_alignment_pdist=scispa.distance.pdist(encoded_focus_alignment,'hamming')
        booleanTF= theta > encoded_focus_alignment_pdist
        #print(booleanTF)
        boolean01_encoded_focus= booleanTF.astype(int)
        #print(boolean01_encoded_focus)
        W = (1./(1+sum(scispa.distance.squareform(boolean01_encoded_focus))))
    Meff=sum(W)
    print(Meff)
#compute the frequencies
    Pij_true = np.zeros((alignment_width, alignment_width, q, q))
    #print(np.shape(Pij_true))
    Pi_true = np.zeros((alignment_width, q))
    #print(np.shape(Pi_true))
    #print(np.shape(Pi_true))
#single-site frequencies
    for j in range(alignment_height):
        for i in range(alignment_width):
            #print(encoded_focus_alignment[j, i])
            Pi_true[i,(int(encoded_focus_alignment[j, i]-1))] += W[j] #increment the proba to have this residue at position i by the weight W(j) of the sequence j considered
#normalization
    Pi_true = Pi_true/Meff
    #print(Pi_true)
#two-site frequencies
    for l in range(int(alignment_height)):
        for e in range(int(alignment_width)):
            for j in range(int(e+1),int(alignment_width)):
                #print(encoded_focus_alignment[l, j])
                Pij_true[e, j, (int(encoded_focus_alignment[l, e])-1),(int(encoded_focus_alignment[l, j])-1)] += W[l]
                Pij_true[j, e, (int(encoded_focus_alignment[l, j])-1),(int(encoded_focus_alignment[l, e])-1)] = Pij_true[e, j, (int(encoded_focus_alignment[l, e])-1),(int(encoded_focus_alignment[l, j])-1)]
    Pij_true = Pij_true/Meff
    scra = np.identity(q)
    for i in range(int(alignment_width)):
        for alpha in range(q):
            for beta in range(q):
                Pij_true[i, i, alpha, beta] = Pi_true[i, alpha] * scra[alpha, beta]

    return(Meff, Pij_true, Pi_true, alignment_width)
def with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q):
    Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*(np.ones((alignment_width, alignment_width, q, q)))
    Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*(np.ones((alignment_width, q)))
    scra = np.identity(q)
    for i in range(alignment_width):
        for alpha in range(1,q):
            for beta in range(1,q):
                Pij[i, i, alpha, beta] = (1.-pseudocount_weight)*Pij_true[i, i, alpha, beta] + pseudocount_weight/q*scra[alpha, beta]
    print(np.amin(Pij))
    print(np.amin(Pi))
    return(Pij, Pi)

def Get_PMIs(Pij, Pi, alignment_width, q):
    C=np.zeros((alignment_width,alignment_width,q,q))
    for i in range(alignment_width):
        for j in range(alignment_width):
            for alpha in range(1,q):
                for beta in range(1,q):
                    if Pij[i, j, alpha, beta]!=0:
                        C[i,j,alpha,beta]= math.log((Pij[i, j, alpha, beta]) / (Pi[i, alpha]*Pi[j, beta]))
                    else:
                        C[i,j,alpha,beta]=-math.inf

    return(C)



#fill the diagonal of Pij_true

if __name__ == "__main__":
    fi= os.sys.argv[1]
    pseudocount_weight=0.15
    theta=0.15
    Compute_PMIs(fi, pseudocount_weight, theta)
