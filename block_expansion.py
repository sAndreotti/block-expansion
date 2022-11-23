import numpy as np
import argparse

# Argomenti
parser = argparse.ArgumentParser()
parser.add_argument("--panel", help="Raw panel", default=".")

args = parser.parse_args()

H1Blocks = []
blocks = []

def inizio():
    block = [[5, 6, 7, 16, 17, 18]]
    block.extend([6, 10])  # inizio e fine blocco
    block.append('00000')  # sequenza
    blocks.append(block)

    block2 = [[1, 2, 3, 11, 12, 13]]
    block2.extend([4, 5])
    block2.append('10')
    blocks.append(block2)


# creo con tutte le sequenze del blocco una matrice in cui la prima riga è anche l'ultima, così posso vedere tutte
# le divergenze di ogni riga con quella sopra, mi basta fare un for per ogni riga, guardo la divergenza a dx e a sx
# e prendo per ognuna la minore
def larghezza():
    for bi in blocks:
        indexs_bi = bi[0].copy()
        ini_bi = bi[1]
        fin_bi = bi[2]

        X = read_matrix(args.panel)
        aplo_mat = np.empty((len(indexs_bi) + 1, len(X[0]))).astype('int32')
        ind_mat = 0
        for aplo in range(len(X)):
            if aplo in indexs_bi:
                aplo_mat[ind_mat, :] = X[aplo]
                ind_mat = ind_mat + 1
        aplo_mat[ind_mat, :] = aplo_mat[0, :]
        rev_aplo_mat = np.fliplr(aplo_mat)

        # Ho creato la matrice
        mat_print(aplo_mat)

        # calcolo indici e divergenze
        order_mat, divergence_mat = pbwt(aplo_mat)
        rev_order_mat, rev_div_mat = pbwt(rev_aplo_mat)
        sx_min = len(aplo_mat[0])
        dx_min = len(aplo_mat[0])

        for element in range(len(aplo_mat)):
            if element != 0:
                # prendo la divergenza verso dx (numero caratteri diversi)
                order_element = index_at(element, ini_bi - 1, order_mat)
                # prendo la divergenza del mio indice al sito inizio-1
                # trasformo la divergenza in "n caratteri uguali"
                sx_div = ini_bi - 1 - divergence_mat[order_element, ini_bi - 1]
                # se il n caratteri uguali è < di min
                if sx_div < sx_min:
                    sx_min = sx_div

                # faccio la stessa cosa a sx
                rev_order_element = index_at(element, fin_bi + 1, rev_order_mat)
                dx_div = fin_bi + 1 - rev_div_mat[rev_order_element, fin_bi + 1]
                if dx_div < dx_min:
                    dx_min = dx_div

        # ho calcolato l'aumento della larghezza
        seq_sx = aplo_mat[0]
        seq_sx = str(seq_sx[ini_bi - 1 - sx_min:ini_bi - 1]) + "x" + str(seq_sx[ini_bi:fin_bi + 1])
        H1Blocks.append([indexs_bi, ini_bi - 1 - sx_min, fin_bi, seq_sx])

        seq_dx = aplo_mat[0]
        seq_dx = str(seq_dx[ini_bi: fin_bi + 1]) + "x" + str(seq_dx[fin_bi + 2: fin_bi + 2 + dx_min])
        H1Blocks.append([indexs_bi, ini_bi, fin_bi + 1 + dx_min, seq_dx])

    mat_print(H1Blocks)


# creo una matrice in cui la prima riga del blocco è la prima e l'ultima, le sequenze sono limitate all'intervallo
# del blocco, per ogni sequenza se div_sx + div_dx = 0 allora posso essere aggiunte (aggiunta la seq prima della
# riga in cui sono zero)
def altezza():
    for bi in blocks:
        indexs_bi = bi[0].copy()
        ini_bi = bi[1]
        fin_bi = bi[2]

        X = read_matrix(args.panel)

        # prendo la sequenza di riferimento
        seq_rif = X[indexs_bi[0]][ini_bi:fin_bi + 1]

        for sito in range(fin_bi + 1 - ini_bi):
            # suppongo sito come sito di errore, quindi ammetto la differenza di caratteri
            aplo_add = []

            for element in range(len(X)):
                if element not in indexs_bi:
                    # creo la matrice seq-blocco, seq-esame
                    aplo_mat = np.empty((2, fin_bi + 1 - ini_bi)).astype('int32')

                    aplo_mat[0, :] = seq_rif
                    aplo_mat[1, :] = X[element][ini_bi:fin_bi + 1]
                    rev_aplo_mat = np.fliplr(aplo_mat)

                    order_mat, divergence_mat = pbwt(aplo_mat)
                    rev_order_mat, rev_div_mat = pbwt(rev_aplo_mat)

                    # prendo sempre le divergenze nella riga 1 poichè sono rispetto alla precedente
                    sx_div = divergence_mat[1, sito]
                    dx_div = rev_div_mat[1, (fin_bi - ini_bi) - sito]

                    if sx_div + dx_div == 0:
                        aplo_add.append(element)

            if len(aplo_add) != 0:
                seq = aplo_mat[0]
                seq = str(seq[:sito]) + "x" + str(seq[sito + 1:])
                h1_indexs = indexs_bi + aplo_add
                H1Blocks.append([h1_indexs, ini_bi, fin_bi, seq])

    mat_print(H1Blocks)


def read_matrix(path):
    X = np.loadtxt(path, dtype='i', delimiter='\t')
    return X


def mat_print(mat):
    print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in mat]))


def pbwt(matrice):
    index_mat = []
    diver_mat = []
    order_mat = np.empty((len(matrice), len(matrice[0]))).astype('int32')
    divergence_mat = np.empty((len(matrice), len(matrice[0]))).astype('int32')

    for element in range(len(matrice[0]) + 1):
        ppa, div = prefix_divergence_arrays(matrice, element)
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

    return order_mat, divergence_mat


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


def index_at(elemento, posizione, matrice_ordinamento):
    for aplo in matrice_ordinamento[:, posizione]:
        if matrice_ordinamento[aplo, posizione] == matrice_ordinamento[elemento, posizione]:
            return aplo


def index_from(elemento, list_index):
    return list_index[elemento]


if __name__ == '__main__':
    inizio()
    larghezza()
    altezza()
