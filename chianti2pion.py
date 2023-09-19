# Mpv10 format convert pypion toolkit.

# This script convert the original mellema cooling table to
# mpv10 format table

# Author: Arun Mathew
# Created: 21-10-2022

# Usage:

import numpy as np
import argparse


def mpv10_format(table,desti_dir, filename):

    print("Reading source ... ")
    Nrows = len(table)
    Ncols = len(table[0])
    print("Source table contain " + str(Nrows) + " rows and " + str(Ncols) + " columns")
    print("Converting to MPv10 format ... ")

    outfile = open(desti_dir + filename + ".cpp", "w")

    outfile.write('chinati_rate_'+ filename+ ' = \n {')

    for i in range(0, Nrows):
        line = "{"
        for j in range(1, Ncols):
            line = line + '{:.3f}'.format(table[i, j])
            if j != (Ncols - 1): line = line + ", "
            if j == (Ncols - 1) and i != Nrows - 1: line = line + "}, "
            if j == (Ncols - 1) and i == Nrows - 1: line = line + "}"

        outfile.write(line)
    outfile.write("}")
    outfile.close()




