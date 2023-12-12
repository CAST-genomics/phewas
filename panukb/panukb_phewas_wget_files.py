#!/usr/bin/env python3
"""
Perform phewas on all Pan-UKB phenotypes

Example:
./panukb_phewas_gwet_files.py
"""

import argparse
import os
import pandas as pd
import sys

# TODO set arg for output location
# TODO don't spit all the output to the terminal
def RunCmd(cmd):
	os.system(cmd)

def MSG(msg):
	sys.stderr.write("[panukb_phewas_wget_files.py]:"+msg.strip()+"\n")

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--manifest", help="Path to manifest file", type=str, default="panukb-manifest.csv")
    args = parser.parse_args()
    
    manifest = pd.read_csv(args.manifest)

    for index, row in manifest.iterrows():
            # Download header line and tabix index
            RunCmd("wget {}".format(row['aws_link_tabix']))

if __name__ == "__main__":
	main()    
    