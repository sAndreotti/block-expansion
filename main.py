import math

import numpy as np
import argparse

# Argomenti
parser = argparse.ArgumentParser()
parser.add_argument("--panel", help="Raw panel", default=".")

args = parser.parse_args()

H1Blocks = []


def print_hi():
    print('Hi')
    blocks = []
    block = []
    block.append([1, 2, 3, 11, 12, 13])  # indici blocco
    block.extend([4, 5])  # inizio e fine blocco
    block.append('10')  # sequenza
    blocks.append(block);
    print(blocks)

    X = read_matrix(args.panel)

    order_mat = np.empty((len(X), len(X[0]))).astype('int32')
    divergence_mat = np.empty((len(X), len(X[0]))).astype('int32')

    order_r_mat = np.empty((len(X), len(X[0]))).astype('int32')
    divergence_r_mat = np.empty((len(X), len(X[0]))).astype('int32')

    print()
    print("Matrix")
    matprint(X)
    print()

    index_mat = []
    diver_mat = []
    index_r_mat = []
    diver_r_mat = []

    for element in range(len(X[0]) + 1):
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

    X_r = np.fliplr(X)
    for element in range(len(X_r[0]) + 1):
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

    # print("Index")
    # matprint(order_mat)
    # print()
    # print()
    # print("Divergence")
    # matprint(divergence_mat)

    #print("Index Reverse")
    #matprint(order_r_mat)
    #print()
    #print()
    #print("Divergence Reverse")
    #matprint(divergence_r_mat)

    block_width(blocks, X, order_mat, order_r_mat, divergence_mat, divergence_r_mat)
    block_height(blocks, order_mat, order_r_mat, X)

def read_matrix(path):
    X = np.loadtxt(path, dtype='i', delimiter='\t')
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
        indexs_bi = bi[0]
        ini_bi = bi[1]
        fin_bi = bi[2]
        seq_bi = bi[3]

        amp_s = len(X[0])
        amp_d = len(X[0])
        ex_ind = indexs_bi[0]
        ex_div = get_diver_pos(indexs_bi[0], ini_bi, order_mat) + 1

        for j in indexs_bi:
            # Espansione verso sinistra
            if ini_bi != 1 and ini_bi != 0:
                index_bi = get_diver_pos(j, ini_bi, order_mat)

                if index_bi + 1 > len(X) - 1:
                    divergenza = ini_bi - divergence_mat[index_bi, ini_bi]
                else:
                    divergenza = ini_bi - divergence_mat[index_bi + 1, ini_bi]

                # Se l'indice che analizzo è staccato dagli altri
                if index_bi + 1 != ex_div and index_bi - 1 != ex_div:
                    # Calcolo divergenza tra le due stringhe, quando leggo una divergenza = 0
                    for i in reversed(range(ini_bi - 1)):
                        ppa, div = prefix_divergence_arrays([X[ex_ind], X[j]], i)
                        if div[1] == 0:
                            divergenza = ini_bi - i - 1
                            break;

                if divergenza < amp_s:
                    amp_s = divergenza

            # Espansione verso destra
            if fin_bi != len(X[0] - 1) and fin_bi != len(X[0] - 2):
                index_bi_r = get_diver_pos(j, fin_bi, order_r_mat)

                if index_bi_r + 1 > len(X) - 1:
                    divergenza = fin_bi - divergence_r_mat[index_bi_r, fin_bi]
                else:
                    divergenza = fin_bi - divergence_r_mat[index_bi_r + 1, fin_bi]

                # Se l'indice che analizzo è staccato dagli altri
                if index_bi_r + 1 != ex_div and index_bi_r - 1 != ex_div:
                    ex_rev = X[ex_ind]
                    j_rev = X[j]
                    # Calcolo divergenza tra le due stringhe, quando leggo una divergenza = 0
                    for i in range(len(X[0]) - fin_bi + 1):
                        ppa, div = prefix_divergence_arrays([ex_rev[::-1], j_rev[::-1]], i)
                        if div[1] == 0:
                            divergenza = fin_bi - i - 1
                            break;

                divergenza = abs(divergenza)
                if divergenza < amp_d:
                    amp_d = divergenza

            ex_ind = j
            ex_div = index_bi

        if amp_s != len(X[0]):
            seq = X[indexs_bi[0]]
            seqS = str(seq[ini_bi - 1 - amp_s:ini_bi - amp_s]) + "x" + str(seq[ini_bi:fin_bi + 1])
            H1Blocks.append([indexs_bi, ini_bi - 1 - amp_s, fin_bi, seqS])

        if amp_d != len(X[0]):
            seq = X[indexs_bi[0]]
            seqD = str(seq[ini_bi:fin_bi + 1]) + "x" + str(seq[fin_bi + 2:fin_bi + 2 + amp_d])
            H1Blocks.append([indexs_bi, ini_bi, fin_bi + 1 + amp_d, seqD])

        if ini_bi == 1:
            H1Blocks.append([indexs_bi, 0, fin_bi, "x" + seq_bi])
        if fin_bi == len(X[0] - 1):
            H1Blocks.append([indexs_bi, ini_bi, len(X[0]), seq_bi + "x"])

        #print()
        #print("H1 max site:")
        #matprint(H1Blocks)

def block_height(blocks, order_mat, order_r_mat, X):
    for bi in blocks:
        indexs_bi = bi[0].copy()
        ini_bi = bi[1]
        fin_bi = bi[2]
        seq_bi = bi[3]

        for sito in range(fin_bi - ini_bi + 1):
            # aggiungo l'offset per essere nel sito giusto
            sito = sito + ini_bi
            add_aplo = []

            for aplo in range(len(X)):
                left_ind = False
                right_ind = False

                if aplo not in indexs_bi:
                    if sito != ini_bi and sito != fin_bi:
                        if X[aplo][ini_bi:sito] == X[indexs_bi[0]][ini_bi:sito]:
                            # la parte sx è giusta
                            left_ind = True
                        if X[aplo][sito+1:fin_bi+1] == X[indexs_bi[0]][sito+1:fin_bi+1]:
                            # la parte dx è giusta
                            right_ind = True
                    elif sito == ini_bi:
                        if X[aplo][sito+1:fin_bi+1] == X[indexs_bi[0]][sito+1:fin_bi+1]:
                            # la parte dx è giusta
                            right_ind = True
                            left_ind = True
                    elif sito == fin_bi:
                        if X[aplo][ini_bi:sito] == X[indexs_bi[0]][ini_bi:sito]:
                            # la parte sx è giusta
                            left_ind = True
                            right_ind = True

                    if left_ind and right_ind:
                        add_aplo.append(aplo)

            if len(add_aplo) != 0:
                seq_h = seq_bi[:sito-ini_bi] + "x" + seq_bi[sito+1-ini_bi:]
                h1_indexs = indexs_bi + add_aplo
                H1Blocks.append([h1_indexs, ini_bi, fin_bi, seq_h])

    print()
    print("H1 Blocks:")
    matprint(H1Blocks)


def matprint(mat):
    print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in mat]))


# ritorna la posizione di index nell'ordinamento a pos_diver
def get_diver_pos(indice_da_trovare, sito, matrice_ordinamento):
    for element in range(len(matrice_ordinamento[:, sito])):
        if matrice_ordinamento[element, sito] == indice_da_trovare:
            return element


if __name__ == '__main__':
    print_hi()
