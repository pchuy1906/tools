#!/usr/bin/env python
import sys
import pandas as pd
import argparse
from pyace.preparedata import WEIGHTS_ENERGY_COLUMN, WEIGHTS_FORCES_COLUMN, normalize_energy_forces_weights
import numpy as np

parser = argparse.ArgumentParser(description="Add weights to pckl.gzip for PACE optimization.")
parser.add_argument('--file_pckl_1',      default="file1.pckl.gzip", help='Input pckl.gzip file')
parser.add_argument('--file_pckl_2',      default="file2.pckl.gzip", help='Input pckl.gzip file')
parser.add_argument('--wE_1', type=float, default=1.0,               help='weight for energy')
parser.add_argument('--wE_2', type=float, default=1.0,               help='weight for energy')

args = parser.parse_args()
file_pckl_1 = args.file_pckl_1
file_pckl_2 = args.file_pckl_2
wE_1 = args.wE_1
wE_2 = args.wE_2

def generate_force_weights(row):
    n = int(row["natoms"])
    return np.ones(n)

df1 = pd.read_pickle(file_pckl_1, compression="gzip")
print (df1.shape)
print (df1.head())
df1[WEIGHTS_ENERGY_COLUMN] = wE_1/df1["natoms"]
df1[WEIGHTS_FORCES_COLUMN] =  df1.apply(generate_force_weights, axis=1)
print (df1.shape)
print (df1.head())

df2 = pd.read_pickle(file_pckl_2, compression="gzip")
print (df2.shape)
print (df2.head())
df2[WEIGHTS_ENERGY_COLUMN] = wE_2/df2["natoms"]
df2[WEIGHTS_FORCES_COLUMN] =  df2.apply(generate_force_weights, axis=1)
print (df2.shape)
print (df2.head())

df = pd.concat([df1, df2], ignore_index=True)

df.to_pickle("all_data.pckl.gzip", compression="gzip", protocol=4)
print (df.head())

