
import sys
import argparse
import gzip
import pandas as pd
import numpy as np
import datetime

#######################################
# Functions
#######################################

#######################################
# Parse and set up argument parameters
#######################################

parser = argparse.ArgumentParser()
parser.add_argument("--dose", "-d", help="Mach input dose file", required=True)
parser.add_argument("--info", "-i", help="Mach input INFO file.", required=True)
parser.add_argument("--rsq", "-r", default = 0.3, help="R squared filter (default = 0.3).", required=True)
parser.add_argument("--out", "-o", help="Output file prefix. Output is gzipped.", required=True)

args = parser.parse_args()
args.rsq = float(args.rsq)

info = pd.read_table(args.info, header = 0, compression = "gzip")

# ix = (info['Rsq'] >= args.rsq) | (info["Genotyped"] == "Genotyped")
ix = (info['Rsq'] >= args.rsq)
info = info.loc[ix]

info_out_fn = args.out + ".info.gz"
info.to_csv(info_out_fn, sep = "\t", header = True, index = False, compression = "gzip")

dosefs = gzip.open(args.dose, 'rt')
dose_out_fn = args.out + ".dose.gz"
dose_out_fs = gzip.open(dose_out_fn, 'wt')
for line in dosefs:
    linel = line.strip().split('\t')
    id = linel.pop(0)
    typ = linel.pop(0)
    print(id)
    linel = pd.Series(linel)
    linel = linel[ix]
    lineol = [id, typ] + linel.tolist()
    lineo = '\t'.join(lineol)
    dose_out_fs.write(lineo + '\n')

dosefs.close()
dose_out_fs.close()

