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

    X = np.loadtxt(args.panel, dtype='i', delimiter='\t')

    for line in lines:
        start = line.index("[")
        end = line.index("]")
        indici = line[start+1:end]

        start = line.index("(")
        end = line.index("-")
        inizio = line[start+1:end]

        start = end
        end = line.index(",")
        fine = line[start+1:end]

        start = line.index(",")
        end = line.index(":")
        panel_aplo = line[start+1:end]
        seq = X[int(panel_aplo)]
        sequenza = seq[int(inizio):int(fine)]

        seq_str = ""
        for element in sequenza:
            seq_str = seq_str + str(element)

        print(f"{indici};{inizio};{fine};{seq_str}")


if __name__ == '__main__':
    parser()