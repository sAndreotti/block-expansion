import argparse
import sys
from os import path

import numpy as np

# Argomenti
parser = argparse.ArgumentParser()
parser.add_argument("--panel", help="Pannello in formato 0/1", default=".")
parser.add_argument("--blocks", help="File contenente i blocchi in input", default=".")

args = parser.parse_args()
blocks = []


def parser():
    if args.panel == '.':
        sys.exit("Pannello non inserito")

    if args.blocks == '.':
        sys.exit("Blocchi non inseriti")

    if not path.isfile(args.panel):
        sys.exit("Pannello non valido")

    if not path.isfile(args.blocks):
        sys.exit("Blocchi non validi")

    f = open(args.blocks, "r")
    lines = f.readlines()

    if len(lines) == 0:
        sys.exit("Blocchi non presenti")

    # grandezza (start-fine, prima-linea:numero sequenze), [indici]
    # indici; inizio; fine; sequenza

    # X = np.loadtxt(args.panel, dtype='i', converters=float)
    p = open(args.panel, "r")
    lines_panel = p.readlines()

    count = -1

    for line in lines:
        count += 1

        if count % 10 and count < 100:
            start = line.index("[")
            end = line.index("]")
            indici = line[start + 1:end]

            start = line.index("(")
            end = line.index("-")
            inizio = line[start + 1:end]

            start = end
            end = line.index(",")
            fine = line[start + 1:end]

            start = line.index(",")
            end = line.index(":")
            panel_aplo = line[start + 1:end]
            seq = lines_panel[int(panel_aplo)]
            sequenza = seq[int(inizio):int(fine)+1]

            seq_str = ""
            for element in sequenza:
                seq_str = seq_str + str(element)

            if inizio != fine:
                print(f"{indici};{inizio};{fine};{seq_str}")


def bmtab():
    bm = open(args.panel, "r")
    lines = bm.readlines()

    n_line = []
    w_lines = []

    for line in range(len(lines)):
        e_list = list(lines[line])
        for i in range(len(e_list)):
            n_line.append(e_list[i])
            if i != len(e_list)-1:
                n_line.append('\t')

        # n_line.append('\n')
        if n_line[-2] == '\t' and n_line[-1] != '0' and n_line[-1] != '1':
            del n_line[-2]

        w_lines.append(''.join(map(str, n_line)))
        n_line = []

    bm.close()

    open(args.panel, "w").close()

    bm = open(args.panel, "w")
    [bm.write(element) for element in w_lines]

    bm.close()


if __name__ == '__main__':
    # parser()
    bmtab()
