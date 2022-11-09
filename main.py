import math

import numpy as np
import argparse

#Argomenti
parser = argparse.ArgumentParser()
parser.add_argument("--panel", help="Raw panel", default=".")

args = parser.parse_args()

H1Blocks = []


def print_hi():
    print('Hi')
    blocks = []
    block = []
    block.append([1, 2, 3, 8, 11, 12, 13]) #indici blocco
    block.extend([4, 5]) #inizio e fine blocco
    block.append('10') #sequenza
    blocks.append(block);
    print(blocks)

    X = read_matrix(args.panel)

    order_mat = np.empty((len(X), len(X[0]))).astype('int32')
    divergence_mat = np.empty((len(X), len(X[0]))).astype('int32')

    order_r_mat = np.empty((len(X), len(X[0]))).astype('int32')
    divergence_r_mat = np.empty((len(X), len(X[0]))).astype('int32')

    #print("Matrix")
    #matprint(X)
    #print()
    #print()

    index_mat = []
    diver_mat = []
    index_r_mat = []
    diver_r_mat = []

    for element in range(len(X[0])+1):
        ppa, div = prefix_divergence_arrays(X, element)
        index_mat.append(ppa)
        diver_mat.append(div)
    del index_mat[0]
    del diver_mat[0]

    index_numpy = np.array([np.array(indexi) for indexi in index_mat])
    for element in range(len(index_numpy)):
        order_mat[:, element] = index_numpy[element]

    diver_numpy = np.array([np.array(indexi) for indexi in diver_mat])
    for element in range(len(diver_numpy)):
        divergence_mat[:, element] = diver_numpy[element]


    print("Index")
    matprint(order_mat)
    print()
    print()
    print("Divergence")
    matprint(divergence_mat)

    X_r = np.fliplr(X)
    for element in range(len(X_r[0])+1):
        ppa, div = prefix_divergence_arrays(X_r, element)
        index_r_mat.append(ppa)
        diver_r_mat.append(div)
    del index_r_mat[0]
    del diver_r_mat[0]

    index_r_numpy = np.array([np.array(indexi) for indexi in index_r_mat])
    for element in range(len(index_r_numpy)):
        order_r_mat[:, element] = index_r_numpy[element]

    diver_r_numpy = np.array([np.array(indexi) for indexi in diver_r_mat])
    for element in range(len(diver_r_numpy)):
        divergence_r_mat[:, element] = diver_r_numpy[element]

    block_width(blocks, X, order_mat, order_r_mat, divergence_mat, divergence_r_mat)

def read_matrix(path):
    X = np.loadtxt(path, dtype='i', delimiter=',')
    return X

def prefix_divergence_arrays(X, pos):
    M = len(X)
    N = pos

    # initialise positional prefix array
    ppa = list(range(M))
    # initialise divergence array
    div = [0] * M

    for k in range(N - 1):

        # setup intermediates
        a = list()
        b = list()
        d = list()
        e = list()
        p = q = k + 1

        # iterate over haplotypes in reverse prefix sorted order
        for index, match_start in zip(ppa, div):
            haplotype = X[index]
            allele = haplotype[k]

            # update intermediates
            if match_start > p:
                p = match_start
            if match_start > q:
                q = match_start

            # update intermediates
            if allele == 0:
                a.append(index)
                d.append(p)
                p = 0
            else:
                b.append(index)
                e.append(q)
                q = 0

        # construct the new arrays for k+1 by concatenating intermediates
        ppa = a + b
        div = d + e

    return ppa, div

def block_width(blocks, X, order_mat, order_r_mat, divergence_mat, divergence_r_mat):

    for bi in blocks:
        indexsBi = bi[0]
        iniBi = bi[1]
        finBi = bi[2]

        amp_s=-1
        amp_d=-1

        for j in indexsBi:
            indexBi = get_diver_pos(j, iniBi, order_mat)
            #print(iniBi-1)
            #print(indexBi)
            print(divergence_mat[indexBi, iniBi])
            divergenza = divergence_mat[indexBi, iniBi]
            if divergenza < amp_s:
                amp_s = divergenza
            elif amp_s == -1:
                amp_s = divergenza

            indexBi_r = get_rev_diver_pos(j, finBi, order_r_mat)
            divergenza = divergence_r_mat[indexBi_r, finBi]
            if divergenza < amp_d:
                amp_d = divergenza
            elif amp_d == -1:
                amp_d = divergenza

        seq = X[indexsBi[0]]
        seqS = str(seq[iniBi - 1 - amp_s:iniBi-amp_s+1]) + "x" + str(seq[iniBi:finBi+1])
        print(amp_s)
        print("Sinistra:")
        print(seqS)
        H1Blocks.append([indexsBi, iniBi-1-amp_s, finBi, seqS])

        seq = X[indexsBi[0]]
        seqD = str(seq[iniBi:finBi+1]) + "x" + str(seq[finBi+2:finBi+2+amp_d])
        print(amp_d)
        print("Destra:")
        print(seqD)
        H1Blocks.append([indexsBi, iniBi, finBi+1+amp_d, seqD])

        print("H1 Computated")
        matprint(H1Blocks)


def matprint(mat):
    print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in mat]))

#ritorna la posizione di index nell'ordinamento a pos_diver
def get_diver_pos(ind, pos_diver, order_mat):
    for element in range(len(order_mat[:, pos_diver])):
        if order_mat[element, pos_diver] == ind:
            return element

#ritorna la posizione di index nell'ordinamento a pos_diver
def get_rev_diver_pos(ind, pos_diver, order_r_mat):
    for element in range(len(order_r_mat[:, pos_diver])):
        if order_r_mat[element, pos_diver] == ind:
            return element

if __name__ == '__main__':
    print_hi()
