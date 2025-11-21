#!/usr/bin/env python
import sys
import pandas as pd
import argparse
from pyace.preparedata import WEIGHTS_ENERGY_COLUMN, WEIGHTS_FORCES_COLUMN, normalize_energy_forces_weights
import numpy as np

parser = argparse.ArgumentParser(description="Add weights to pckl.gzip for PACE optimization.")
parser.add_argument('input',            help='Input pckl.gzip file')
parser.add_argument('--wE', type=float, help='weight for energy')

args = parser.parse_args()
filename = args.input
wE = args.wE
df=pd.read_pickle(filename, compression="gzip")
print (df.shape)
print (df.head())
outfile = filename.replace(".pckl.gzip", "_weights.pckl.gzip")

df[WEIGHTS_ENERGY_COLUMN] = wE/df["natoms"]

def generate_force_weights(row):
    n = int(row["natoms"])
    return np.ones(n)

df[WEIGHTS_FORCES_COLUMN] =  df.apply(generate_force_weights, axis=1)

#normalize_energy_forces_weights(df)
df.to_pickle(outfile, compression="gzip", protocol=4)

print (df.head())

